"""
2I-network simulation: face-sharing and axis-sharing growth of icosahedral
networks in R^3, with observables (Laplacian spectral gap, Hausdorff dimension,
heat-kernel return) and a random-regular-graph null.

Pre-registered per 2I-network/PRE_REGISTRATION.md.
"""

from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass

import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import eigsh, ArpackNoConvergence
from scipy.sparse.csgraph import shortest_path, connected_components

PHI = (1.0 + np.sqrt(5.0)) / 2.0
SEED = 20260525


# ---------------------------------------------------------------------------
# Standard icosahedron geometry
# ---------------------------------------------------------------------------

def standard_icosahedron() -> np.ndarray:
    raw = np.array([
        [0,  1,  PHI], [0,  1, -PHI], [0, -1,  PHI], [0, -1, -PHI],
        [1,  PHI, 0], [1, -PHI, 0], [-1,  PHI, 0], [-1, -PHI, 0],
        [PHI, 0,  1], [PHI, 0, -1], [-PHI, 0,  1], [-PHI, 0, -1],
    ], dtype=np.float64)
    return raw / np.sqrt(1.0 + PHI * PHI)


V_STD = standard_icosahedron()
EDGE_LEN = 2.0 / np.sqrt(1.0 + PHI * PHI)


def compute_edges_faces():
    edges: set[tuple[int, int]] = set()
    for i in range(12):
        for j in range(i + 1, 12):
            d = np.linalg.norm(V_STD[i] - V_STD[j])
            if abs(d - EDGE_LEN) < 1e-6:
                edges.add((i, j))
    faces: list[tuple[int, int, int]] = []
    for i in range(12):
        for j in range(i + 1, 12):
            if (i, j) not in edges:
                continue
            for k in range(j + 1, 12):
                if (i, k) in edges and (j, k) in edges:
                    faces.append((i, j, k))
    # antipodes for axes
    antipodes = {}
    for i in range(12):
        for j in range(12):
            if np.linalg.norm(V_STD[i] + V_STD[j]) < 1e-9:
                antipodes[i] = j
                break
    axes: list[tuple[int, int]] = []
    seen = set()
    for i, j in antipodes.items():
        if i not in seen:
            axes.append((i, j))
            seen.add(i)
            seen.add(j)
    return list(edges), faces, axes


EDGES_STD, FACES_STD, AXES_STD = compute_edges_faces()
assert len(EDGES_STD) == 30, len(EDGES_STD)
assert len(FACES_STD) == 20, len(FACES_STD)
assert len(AXES_STD) == 6, len(AXES_STD)


def face_normal_and_centroid(face_verts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c = face_verts.mean(axis=0)
    v1 = face_verts[1] - face_verts[0]
    v2 = face_verts[2] - face_verts[0]
    n = np.cross(v1, v2)
    n = n / np.linalg.norm(n)
    return n, c


def reflect_points(pts: np.ndarray, plane_pt: np.ndarray, plane_n: np.ndarray) -> np.ndarray:
    return pts - 2.0 * np.outer(np.dot(pts - plane_pt, plane_n), plane_n)


# ---------------------------------------------------------------------------
# Vertex pool with spatial hashing
# ---------------------------------------------------------------------------

class VertexPool:
    def __init__(self, eps: float = 1e-6):
        self.eps = eps
        self.cell = max(eps * 10.0, 1e-5)
        self.cells: dict[tuple[int, int, int], list[int]] = {}
        self.positions: list[np.ndarray] = []

    def add_or_get(self, pos: np.ndarray) -> int:
        key = (int(np.floor(pos[0] / self.cell)),
               int(np.floor(pos[1] / self.cell)),
               int(np.floor(pos[2] / self.cell)))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    k = (key[0] + dx, key[1] + dy, key[2] + dz)
                    for idx in self.cells.get(k, ()):
                        if np.linalg.norm(self.positions[idx] - pos) < self.eps:
                            return idx
        new_idx = len(self.positions)
        self.positions.append(pos.copy())
        self.cells.setdefault(key, []).append(new_idx)
        return new_idx


# ---------------------------------------------------------------------------
# Face-sharing growth
# ---------------------------------------------------------------------------

def build_face_shared(N_target: int) -> tuple[list[list[int]], VertexPool, set[tuple[int, int]]]:
    """
    Returns:
        icos_idx: list of 12-tuples of global vertex indices, one per icosahedron.
        pool: VertexPool with all positions.
        edges_graph: set of (i, j) with i < j giving the network edges
                      (two icosahedra share at least one face = 3 vertices).
    """
    pool = VertexPool(eps=1e-6)
    icos_local_verts: list[np.ndarray] = []
    icos_idx: list[list[int]] = []

    # Initial icosahedron
    init_pos = V_STD.copy()
    init_idx = [pool.add_or_get(p) for p in init_pos]
    icos_local_verts.append(init_pos)
    icos_idx.append(init_idx)

    face_open: list[list[bool]] = [[True] * 20]
    # Map sorted-triple-of-global-indices -> first icosahedron to claim it
    face_map: dict[tuple[int, int, int], tuple[int, int]] = {}
    for f, (a, b, c) in enumerate(FACES_STD):
        tri = tuple(sorted((init_idx[a], init_idx[b], init_idx[c])))
        face_map[tri] = (0, f)

    queue: deque[tuple[int, int]] = deque((0, f) for f in range(20))
    edges_graph: set[tuple[int, int]] = set()

    while len(icos_idx) < N_target and queue:
        ico_i, f_i = queue.popleft()
        if not face_open[ico_i][f_i]:
            continue

        parent_verts = icos_local_verts[ico_i]
        face_v_pos = parent_verts[list(FACES_STD[f_i])]
        n, c = face_normal_and_centroid(face_v_pos)
        new_pos = reflect_points(parent_verts, c, n)

        new_glob_idx = [pool.add_or_get(p) for p in new_pos]
        # Sanity: the three face vertices must coincide with parent's
        for a in FACES_STD[f_i]:
            assert new_glob_idx[a] == icos_idx[ico_i][a], (
                f"Reflection failed to preserve face vertex {a}")

        new_ico = len(icos_idx)
        icos_local_verts.append(new_pos)
        icos_idx.append(new_glob_idx)
        face_open.append([True] * 20)

        # Process the 20 faces of new icosahedron
        for f_new in range(20):
            a, b, c2 = FACES_STD[f_new]
            tri = tuple(sorted((new_glob_idx[a], new_glob_idx[b], new_glob_idx[c2])))
            if tri in face_map:
                other_ico, other_face = face_map[tri]
                # Close both faces
                face_open[new_ico][f_new] = False
                face_open[other_ico][other_face] = False
                if other_ico != new_ico:
                    e = (min(other_ico, new_ico), max(other_ico, new_ico))
                    edges_graph.add(e)
            else:
                face_map[tri] = (new_ico, f_new)
                queue.append((new_ico, f_new))

    return icos_idx, pool, edges_graph


# ---------------------------------------------------------------------------
# Axis-sharing growth (secondary construction)
# ---------------------------------------------------------------------------

def build_vertex_shared(N_target: int) -> tuple[list[list[int]], VertexPool, set[tuple[int, int]]]:
    """
    Vertex-sharing (tertiary, diagnostic): each new icosahedron shares ONE
    vertex with its parent. Canonical placement: new icosahedron is parent
    translated by 2 * direction_k, where direction_k is the unit vector from
    parent center to its k-th vertex. The new icosahedron's antipode-of-k
    vertex coincides with parent's k-th vertex.
    """
    pool = VertexPool(eps=1e-6)
    icos_local_verts: list[np.ndarray] = []
    icos_idx: list[list[int]] = []
    icos_centers: list[np.ndarray] = []

    init_pos = V_STD.copy()
    init_idx = [pool.add_or_get(p) for p in init_pos]
    icos_local_verts.append(init_pos)
    icos_idx.append(init_idx)
    icos_centers.append(np.zeros(3))

    vert_open: list[list[bool]] = [[True] * 12]
    # antipode of vertex k in V_STD
    antipode = {}
    for i in range(12):
        for j in range(12):
            if np.linalg.norm(V_STD[i] + V_STD[j]) < 1e-9:
                antipode[i] = j
                break
    vert_map: dict[int, list[tuple[int, int]]] = {}
    for k in range(12):
        vert_map.setdefault(init_idx[k], []).append((0, k))

    queue: deque[tuple[int, int]] = deque((0, k) for k in range(12))
    edges_graph: set[tuple[int, int]] = set()

    while len(icos_idx) < N_target and queue:
        ico_i, k_i = queue.popleft()
        if not vert_open[ico_i][k_i]:
            continue

        parent_verts = icos_local_verts[ico_i]
        c_parent = icos_centers[ico_i]
        v_shared = parent_verts[k_i]
        direction_k = (v_shared - c_parent) / np.linalg.norm(v_shared - c_parent)
        # New center: along +direction_k beyond v_shared by one radius
        c_new = v_shared + direction_k * 1.0  # circumradius = 1
        new_pos = c_new + (parent_verts - c_parent)  # same orientation, translated

        new_glob_idx = [pool.add_or_get(p) for p in new_pos]
        # The new icosahedron's antipode-of-k vertex should equal parent's k
        ant_k = antipode[k_i]
        assert new_glob_idx[ant_k] == icos_idx[ico_i][k_i], (
            f"vertex-share placement broke: {new_glob_idx[ant_k]} vs {icos_idx[ico_i][k_i]}")

        new_ico = len(icos_idx)
        icos_local_verts.append(new_pos)
        icos_idx.append(new_glob_idx)
        icos_centers.append(c_new)
        vert_open.append([True] * 12)

        # Update vert_map and detect any vertex coincidence with other icosahedra
        for k_new in range(12):
            g = new_glob_idx[k_new]
            others = vert_map.get(g, [])
            if others:
                for other_ico, other_k in others:
                    if other_ico == new_ico:
                        continue
                    e = (min(other_ico, new_ico), max(other_ico, new_ico))
                    edges_graph.add(e)
            vert_map.setdefault(g, []).append((new_ico, k_new))
            # mark the matching pair as closed (the shared vertex with parent)
            if k_new == ant_k:
                vert_open[new_ico][k_new] = False
                vert_open[ico_i][k_i] = False
            else:
                queue.append((new_ico, k_new))

    return icos_idx, pool, edges_graph


def build_axis_shared(N_target: int) -> tuple[list[list[int]], VertexPool, set[tuple[int, int]]]:
    """
    Axis-sharing: each new icosahedron attaches to a parent by sharing one
    5-fold axis (2 antipodal vertices). The new icosahedron is obtained by
    rotating the parent by 180 degrees about the shared axis (which flips
    the icosahedron through the axis). This gives a deterministic placement.
    """
    pool = VertexPool(eps=1e-6)
    icos_local_verts: list[np.ndarray] = []
    icos_idx: list[list[int]] = []

    init_pos = V_STD.copy()
    init_idx = [pool.add_or_get(p) for p in init_pos]
    icos_local_verts.append(init_pos)
    icos_idx.append(init_idx)

    axis_open: list[list[bool]] = [[True] * 6]
    # Map sorted-pair-of-global-indices -> (ico, axis_idx)
    axis_map: dict[tuple[int, int], tuple[int, int]] = {}
    for a_idx, (i, j) in enumerate(AXES_STD):
        pair = tuple(sorted((init_idx[i], init_idx[j])))
        axis_map[pair] = (0, a_idx)

    queue: deque[tuple[int, int]] = deque((0, a) for a in range(6))
    edges_graph: set[tuple[int, int]] = set()

    while len(icos_idx) < N_target and queue:
        ico_i, a_i = queue.popleft()
        if not axis_open[ico_i][a_i]:
            continue

        parent_verts = icos_local_verts[ico_i]
        vi, vj = AXES_STD[a_i]
        # Axis direction
        p0 = parent_verts[vi]
        p1 = parent_verts[vj]
        axis_pt = 0.5 * (p0 + p1)
        axis_dir = p1 - p0
        axis_dir = axis_dir / np.linalg.norm(axis_dir)

        # Rotation by 180 degrees about the axis through axis_pt with direction axis_dir
        # For each point v: v_new = 2 * (proj on axis) - v + 2 * axis_pt
        # i.e., reflection through the axis line
        rel = parent_verts - axis_pt
        proj_lengths = rel @ axis_dir
        proj = np.outer(proj_lengths, axis_dir)
        rotated = 2.0 * proj - rel + axis_pt

        new_glob_idx = [pool.add_or_get(p) for p in rotated]
        # Sanity: axis vertices preserved
        assert new_glob_idx[vi] == icos_idx[ico_i][vi]
        assert new_glob_idx[vj] == icos_idx[ico_i][vj]

        new_ico = len(icos_idx)
        icos_local_verts.append(rotated)
        icos_idx.append(new_glob_idx)
        axis_open.append([True] * 6)

        for a_new in range(6):
            vi2, vj2 = AXES_STD[a_new]
            pair = tuple(sorted((new_glob_idx[vi2], new_glob_idx[vj2])))
            if pair in axis_map:
                other_ico, other_axis = axis_map[pair]
                axis_open[new_ico][a_new] = False
                axis_open[other_ico][other_axis] = False
                if other_ico != new_ico:
                    e = (min(other_ico, new_ico), max(other_ico, new_ico))
                    edges_graph.add(e)
            else:
                axis_map[pair] = (new_ico, a_new)
                queue.append((new_ico, a_new))

    return icos_idx, pool, edges_graph


# ---------------------------------------------------------------------------
# Network analytics
# ---------------------------------------------------------------------------

@dataclass
class NetworkResult:
    N: int
    n_edges: int
    avg_degree: float
    n_components: int
    largest_component_size: int
    lambda1_raw: float
    lambda1_per_deg: float
    lambda1_x_N23: float
    dim_d: float
    dim_residual: float
    heat_alpha: float | None
    eigenvalues: np.ndarray  # smallest 8


def network_to_graph(icos_idx: list[list[int]], edges: set[tuple[int, int]]) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(range(len(icos_idx)))
    G.add_edges_from(edges)
    return G


def compute_laplacian_eigs(G: nx.Graph, k: int = 8) -> np.ndarray:
    n = G.number_of_nodes()
    if n <= 12:
        L = nx.laplacian_matrix(G).toarray().astype(np.float64)
        evals = np.sort(np.linalg.eigvalsh(L))
        return evals[: min(k, n)]
    L = nx.laplacian_matrix(G).astype(np.float64)
    try:
        # sigma=0 shift-invert for smallest eigenvalues; pad k if requested > n-2
        k_req = min(k, n - 2)
        evals = eigsh(L, k=k_req, sigma=0.0, which='LM', return_eigenvectors=False,
                      tol=1e-8, maxiter=5000)
        return np.sort(evals)
    except (ArpackNoConvergence, ValueError, RuntimeError):
        # Fall back: dense
        evals = np.sort(np.linalg.eigvalsh(L.toarray()))
        return evals[: min(k, n)]


def spectral_gap(evals: np.ndarray) -> float:
    # smallest non-zero eigenvalue (>= 1e-9)
    nz = evals[evals > 1e-9]
    return float(nz[0]) if len(nz) > 0 else float('nan')


def hausdorff_dimension_fit(G: nx.Graph, n_sources: int = 200, seed: int = SEED) -> tuple[float, float]:
    """Return (d, residual_RMS) from log-log linear regression of N(R) vs R."""
    n = G.number_of_nodes()
    if n < 6:
        return float('nan'), float('nan')
    rng = np.random.default_rng(seed)
    sources = rng.choice(n, size=min(n_sources, n), replace=False)
    # Compute mean N(R) over sources
    max_dist = 0
    nrs: list[np.ndarray] = []
    for s in sources:
        d = nx.single_source_shortest_path_length(G, int(s))
        if not d:
            continue
        vec = np.zeros(max(max(d.values()), max_dist) + 1, dtype=np.int64)
        for _, v in d.items():
            if v < len(vec):
                vec[v] += 1
        # cumulative
        cum = np.cumsum(vec)
        nrs.append(cum)
        max_dist = max(max_dist, len(cum) - 1)
    if not nrs:
        return float('nan'), float('nan')
    # Pad to same length
    M = max(len(v) for v in nrs)
    A = np.zeros((len(nrs), M), dtype=np.float64)
    for i, v in enumerate(nrs):
        A[i, : len(v)] = v
        A[i, len(v):] = v[-1] if len(v) > 0 else 0
    mean_NR = A.mean(axis=0)
    # R values: 1..M-1 (skip R=0 since N(0)=1)
    R_vals = np.arange(1, M)
    N_vals = mean_NR[1:M]
    valid = N_vals > 0
    R_vals = R_vals[valid]
    N_vals = N_vals[valid]
    if len(R_vals) < 4:
        return float('nan'), float('nan')
    # Middle 60% of range to avoid small-R and finite-size effects
    lo = int(0.2 * len(R_vals))
    hi = int(0.8 * len(R_vals))
    if hi - lo < 3:
        lo, hi = 0, len(R_vals)
    Rfit = R_vals[lo:hi]
    Nfit = N_vals[lo:hi]
    logR = np.log(Rfit)
    logN = np.log(Nfit)
    if np.std(logR) < 1e-9:
        return float('nan'), float('nan')
    slope, intercept = np.polyfit(logR, logN, 1)
    pred = slope * logR + intercept
    rms = float(np.sqrt(np.mean((logN - pred) ** 2)))
    return float(slope), rms


def heat_kernel_alpha(G: nx.Graph, eigs_all: np.ndarray | None = None,
                      t_range: tuple[float, float] = (1.0, 50.0)) -> float | None:
    """
    Compute p(t) = (1/N) sum_k exp(-t * lambda_k) and fit alpha from p(t) ~ t^(-alpha).
    Uses all eigenvalues if available; else estimates from a sparse subset.
    """
    n = G.number_of_nodes()
    if n < 30:
        return None
    try:
        if eigs_all is None:
            L = nx.laplacian_matrix(G).astype(np.float64)
            if n <= 1500:
                eigs_all = np.linalg.eigvalsh(L.toarray())
            else:
                # too expensive; return None
                return None
        ts = np.geomspace(t_range[0], t_range[1], 30)
        p_ts = np.array([np.mean(np.exp(-t * eigs_all)) for t in ts])
        # Subtract floor (long-time = 1/N from zero-mode)
        floor = 1.0 / n
        p_signal = p_ts - floor
        valid = p_signal > 1e-12
        if valid.sum() < 5:
            return None
        slope, _ = np.polyfit(np.log(ts[valid]), np.log(p_signal[valid]), 1)
        return float(-slope)
    except Exception:
        return None


def analyse_graph(G: nx.Graph, N_label: int, want_heat: bool = True) -> NetworkResult:
    n = G.number_of_nodes()
    m = G.number_of_edges()
    avg_d = (2.0 * m) / n if n > 0 else 0.0
    A = nx.to_scipy_sparse_array(G).astype(np.float64)
    n_components, _ = connected_components(A, directed=False)
    sizes = [len(c) for c in nx.connected_components(G)]
    largest = max(sizes) if sizes else 0

    # Work on largest connected component for spectral
    if n_components > 1:
        big = max(nx.connected_components(G), key=len)
        Gbig = G.subgraph(big).copy()
        # relabel
        Gbig = nx.convert_node_labels_to_integers(Gbig)
    else:
        Gbig = G

    eigs = compute_laplacian_eigs(Gbig, k=8)
    lam1 = spectral_gap(eigs)
    d, res = hausdorff_dimension_fit(Gbig)
    alpha = heat_kernel_alpha(Gbig) if want_heat else None

    return NetworkResult(
        N=n, n_edges=m, avg_degree=avg_d,
        n_components=n_components, largest_component_size=largest,
        lambda1_raw=lam1,
        lambda1_per_deg=lam1 / avg_d if avg_d > 0 else float('nan'),
        lambda1_x_N23=lam1 * (n ** (2.0 / 3.0)),
        dim_d=d, dim_residual=res,
        heat_alpha=alpha,
        eigenvalues=eigs,
    )


# ---------------------------------------------------------------------------
# Null model: random regular graph
# ---------------------------------------------------------------------------

def null_distribution(N: int, target_avg_deg: float, n_samples: int = 20,
                      seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate n_samples random regular graphs with N nodes and degree
    closest to target_avg_deg (even N*d required by networkx).
    Returns (lambda1 samples, d samples).
    """
    d = int(round(target_avg_deg))
    if (N * d) % 2 != 0:
        d += 1
    if d < 2:
        d = 2
    if d >= N:
        d = N - 1 if N % 2 == 0 else N - 2

    lam1s = []
    ds = []
    rng = np.random.default_rng(seed)
    for s in range(n_samples):
        sub_seed = int(rng.integers(0, 2**31 - 1))
        try:
            G = nx.random_regular_graph(d, N, seed=sub_seed)
        except nx.NetworkXError:
            continue
        if not nx.is_connected(G):
            comps = nx.connected_components(G)
            big = max(comps, key=len)
            G = G.subgraph(big).copy()
            G = nx.convert_node_labels_to_integers(G)
        eigs = compute_laplacian_eigs(G, k=4)
        lam1s.append(spectral_gap(eigs))
        dim, _ = hausdorff_dimension_fit(G, n_sources=min(100, G.number_of_nodes()))
        ds.append(dim)
    return np.array(lam1s), np.array(ds)


# ---------------------------------------------------------------------------
# Sweep driver
# ---------------------------------------------------------------------------

def run_sweep(N_targets: list[int], construction: str = 'face',
              run_null: bool = True, verbose: bool = True) -> dict:
    """Run the full sweep for a single construction. Returns dict of results."""
    results: dict = {}
    if construction == 'face':
        builder = build_face_shared
    elif construction == 'axis':
        builder = build_axis_shared
    elif construction == 'vertex':
        builder = build_vertex_shared
    else:
        raise ValueError(construction)

    t0 = time.time()
    icos_idx, pool, edges = builder(max(N_targets))
    t_build = time.time() - t0
    if verbose:
        print(f"[{construction}] built up to N={len(icos_idx)} in {t_build:.1f}s, "
              f"{len(pool.positions)} vertices, {len(edges)} edges")

    G_full = network_to_graph(icos_idx, edges)

    for N in N_targets:
        if N > len(icos_idx):
            if verbose:
                print(f"  N={N}: only got {len(icos_idx)} icosahedra; skipping")
            continue
        # Take subgraph on first N nodes (BFS-ordered by construction)
        Gn = G_full.subgraph(range(N)).copy()
        Gn = nx.convert_node_labels_to_integers(Gn)
        t_a = time.time()
        res = analyse_graph(Gn, N, want_heat=(N <= 1500))
        t_anal = time.time() - t_a
        results[N] = {
            'signal': res,
            't_analyse': t_anal,
        }
        if verbose:
            print(f"  N={N}: edges={res.n_edges} avg_deg={res.avg_degree:.2f} "
                  f"comp={res.n_components} biggest={res.largest_component_size} "
                  f"lam1={res.lambda1_raw:.4f} lam1*N^2/3={res.lambda1_x_N23:.2f} "
                  f"d={res.dim_d:.3f} ({t_anal:.1f}s)")

        if run_null:
            t_n = time.time()
            lam_null, d_null = null_distribution(N, res.avg_degree, n_samples=20)
            t_null = time.time() - t_n
            results[N]['null_lam1'] = lam_null
            results[N]['null_d'] = d_null
            z = (res.lambda1_raw - lam_null.mean()) / (lam_null.std() + 1e-12)
            results[N]['z_lambda1'] = z
            if verbose:
                print(f"  N={N} null: lam1={lam_null.mean():.4f}+-{lam_null.std():.4f} "
                      f"d={np.nanmean(d_null):.3f}+-{np.nanstd(d_null):.3f} z={z:.2f} ({t_null:.1f}s)")

    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--construction', choices=['face', 'axis', 'vertex'], default='face')
    parser.add_argument('--N', type=int, nargs='+', default=[12, 60, 120])
    parser.add_argument('--no-null', action='store_true')
    args = parser.parse_args()
    run_sweep(args.N, construction=args.construction, run_null=not args.no_null)
