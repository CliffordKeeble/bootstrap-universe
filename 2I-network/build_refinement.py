"""
v0.2 refinement implementations: iterated edge subdivision (Scheme A) and
barycentric cell subdivision (Scheme B) of the 600-cell graph.

Per pre-registration PRE_REGISTRATION_v0_2.md.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass

import numpy as np
import networkx as nx
from scipy.sparse.linalg import eigsh, ArpackNoConvergence
from scipy.sparse.csgraph import connected_components

from reference_600cell import build_600cell_graph, build_600cell_vertices, PHI


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

V_S3_2I = 2.0 * np.pi * np.pi / 120.0           # vol(S^3/2I), unit S^3 convention
TARGET_LAMBDA1_X_N23 = 168.0 * V_S3_2I ** (2.0 / 3.0)  # ≈ 50.436
SEED = 20260526


# ---------------------------------------------------------------------------
# Scheme A: iterated edge subdivision
# ---------------------------------------------------------------------------

def edge_subdivide(G: nx.Graph) -> nx.Graph:
    """Return new graph: each edge (u, v) -> (u, w), (w, v), with w fresh."""
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    # New vertex IDs start after max existing
    next_id = max(G.nodes()) + 1 if G.nodes() else 0
    for (u, v) in G.edges():
        w = next_id
        next_id += 1
        H.add_node(w)
        H.add_edge(u, w)
        H.add_edge(w, v)
    return H


def scheme_A_levels(max_level: int = 4) -> list[nx.Graph]:
    G0 = build_600cell_graph()
    G0 = nx.convert_node_labels_to_integers(G0)
    graphs = [G0]
    for lvl in range(max_level):
        Gnext = edge_subdivide(graphs[-1])
        graphs.append(Gnext)
    return graphs


# ---------------------------------------------------------------------------
# Scheme B: barycentric cell subdivision of the 600-cell
# ---------------------------------------------------------------------------

def find_600cell_tetrahedra() -> list[tuple[int, int, int, int]]:
    """
    Identify the 600 tetrahedral cells of the 600-cell.

    A tetrahedral cell is a 4-vertex subset where every pair of vertices is
    connected by a 600-cell edge AND the four vertices share a common centroid
    that is the cell center (at distance r_cell from origin on R^4).

    The 600-cell has 600 tetrahedral cells; each cell is a 4-clique in the
    edge graph. So we enumerate all 4-cliques.
    """
    V = build_600cell_vertices()
    G = build_600cell_graph()
    G = nx.convert_node_labels_to_integers(G)
    # Find all 4-cliques. The 600-cell has exactly 600 of them.
    cliques4 = [tuple(sorted(c)) for c in nx.find_cliques(G) if len(c) == 4]
    # Some 4-cliques might be subsumed by larger cliques (find_cliques returns
    # maximal cliques only). For the 600-cell, maximal cliques are 4-cliques
    # (no 5-cliques since each vertex link is an icosahedron with no K_5).
    cliques4 = list(set(cliques4))
    assert len(cliques4) == 600, f"Expected 600 cells, got {len(cliques4)}"
    return cliques4


def scheme_B_level1() -> nx.Graph:
    """Barycentric cell subdivision: add a vertex per tetrahedral cell."""
    G = build_600cell_graph()
    G = nx.convert_node_labels_to_integers(G)
    cells = find_600cell_tetrahedra()
    H = G.copy()
    next_id = max(H.nodes()) + 1
    for cell in cells:
        c = next_id
        next_id += 1
        H.add_node(c)
        for v in cell:
            H.add_edge(c, v)
    return H


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------

@dataclass
class LevelResult:
    scheme: str
    level: int
    N: int
    E: int
    avg_degree: float
    n_components: int
    lambda1_raw: float
    lambda1_x_N23: float
    deviation_pct: float
    eigenvalues: np.ndarray


def smallest_eigs(G: nx.Graph, k: int = 6) -> np.ndarray:
    n = G.number_of_nodes()
    L = nx.laplacian_matrix(G).astype(np.float64)
    if n <= 1500:
        evals = np.sort(np.linalg.eigvalsh(L.toarray()))
        return evals[: min(k, n)]
    k_req = min(k, n - 2)
    try:
        evals = eigsh(L, k=k_req, sigma=0.0, which='LM',
                      return_eigenvectors=False, tol=1e-9, maxiter=10000)
        return np.sort(evals)
    except (ArpackNoConvergence, RuntimeError, ValueError) as e:
        # Try without shift-invert
        evals = eigsh(L, k=k_req, which='SM',
                      return_eigenvectors=False, tol=1e-9, maxiter=20000)
        return np.sort(evals)


def gap(evals: np.ndarray) -> float:
    nz = evals[evals > 1e-9]
    return float(nz[0]) if len(nz) > 0 else float('nan')


def analyse_level(G: nx.Graph, scheme: str, level: int) -> LevelResult:
    n = G.number_of_nodes()
    m = G.number_of_edges()
    avg_d = 2.0 * m / n
    # Connectivity
    A = nx.to_scipy_sparse_array(G).astype(np.float64)
    n_comp, _ = connected_components(A, directed=False)
    if n_comp > 1:
        big = max(nx.connected_components(G), key=len)
        G_use = G.subgraph(big).copy()
        G_use = nx.convert_node_labels_to_integers(G_use)
    else:
        G_use = G
    eigs = smallest_eigs(G_use, k=6)
    lam1 = gap(eigs)
    lam1_x = lam1 * (n ** (2.0 / 3.0))
    dev = (lam1_x - TARGET_LAMBDA1_X_N23) / TARGET_LAMBDA1_X_N23 * 100.0
    return LevelResult(
        scheme=scheme, level=level, N=n, E=m, avg_degree=avg_d,
        n_components=n_comp,
        lambda1_raw=lam1,
        lambda1_x_N23=lam1_x,
        deviation_pct=dev,
        eigenvalues=eigs,
    )


# ---------------------------------------------------------------------------
# Null
# ---------------------------------------------------------------------------

def null_distribution(G: nx.Graph, n_samples: int = 20, seed: int = SEED) -> np.ndarray:
    """
    Generate n_samples random graphs matched on (N, E, degree sequence).
    Return array of lambda1*N^(2/3) values from the largest component.
    """
    n = G.number_of_nodes()
    degree_seq = [d for _, d in G.degree()]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_samples):
        sub_seed = int(rng.integers(0, 2**31 - 1))
        # Try expected_degree_graph first
        try:
            Gn = nx.expected_degree_graph(degree_seq, seed=sub_seed, selfloops=False)
            # Remove self-loops just in case
            Gn.remove_edges_from(nx.selfloop_edges(Gn))
        except Exception:
            # Fallback: configuration model
            Gn = nx.configuration_model(degree_seq, seed=sub_seed)
            Gn = nx.Graph(Gn)
            Gn.remove_edges_from(nx.selfloop_edges(Gn))
        # Largest component
        if not nx.is_connected(Gn):
            big = max(nx.connected_components(Gn), key=len)
            Gn = Gn.subgraph(big).copy()
            Gn = nx.convert_node_labels_to_integers(Gn)
        eigs = smallest_eigs(Gn, k=4)
        lam = gap(eigs)
        out.append(lam * (n ** (2.0 / 3.0)))
    return np.array(out)


if __name__ == '__main__':
    print(f"Target lambda_1 * N^(2/3) = 168 * (2pi^2/120)^(2/3) = {TARGET_LAMBDA1_X_N23:.4f}")
