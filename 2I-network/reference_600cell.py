"""
Reference computation: build the 600-cell explicitly as a 12-regular graph
on 120 vertices and compute its Laplacian spectrum. This is the known-answer
target the brief mentions: 'compute or look up the 600-cell's spectrum'.

The 120 vertices of the 600-cell on S^3 (unit 3-sphere in R^4):
  - 8 permutations of (+/-1, 0, 0, 0)
  - 16 sign-choices of (+/-1/2, +/-1/2, +/-1/2, +/-1/2)
  - 96 even permutations of (0, +/-1/2, +/-phi/2, +/-1/(2 phi))

Edges of the 600-cell connect vertices at chord distance 2 sin(pi/10), i.e.,
at unit-sphere dot product cos(pi/5) = phi/2.
"""

from __future__ import annotations

import itertools
import numpy as np
import networkx as nx
from scipy.sparse.linalg import eigsh

PHI = (1.0 + np.sqrt(5.0)) / 2.0


def even_permutations(t: tuple) -> list[tuple]:
    perms = list(itertools.permutations(t))
    # An even permutation has signature +1; built-in permutations of length 4
    # we'd need parity. Use insertion-sort-style sign or itertools.permutations
    # filtered by parity.
    def parity(p):
        s = 0
        a = list(p)
        for i in range(len(a)):
            for j in range(i + 1, len(a)):
                if a[i] > a[j]:
                    s += 1
        return s % 2

    # to distinguish identical entries, use indices
    n = len(t)
    base_indices = list(range(n))
    even = set()
    for p in itertools.permutations(base_indices):
        if parity(p) == 0:
            even.add(tuple(t[i] for i in p))
    return list(even)


def build_600cell_vertices() -> np.ndarray:
    verts = set()
    # 8 permutations of (1,0,0,0) with signs
    for i in range(4):
        for s in (+1, -1):
            v = [0.0] * 4
            v[i] = float(s)
            verts.add(tuple(v))
    # 16 sign-choices of (1/2, 1/2, 1/2, 1/2)
    for signs in itertools.product([-0.5, 0.5], repeat=4):
        verts.add(tuple(signs))
    # 96 even permutations of (0, +/-1/2, +/-phi/2, +/-1/(2 phi))
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            for s3 in (+1, -1):
                tup = (0.0, s1 * 0.5, s2 * PHI / 2.0, s3 / (2.0 * PHI))
                for p in even_permutations(tup):
                    verts.add(p)
    V = np.array(sorted(verts), dtype=np.float64)
    assert V.shape == (120, 4), V.shape
    # Verify unit sphere
    norms = np.linalg.norm(V, axis=1)
    assert np.allclose(norms, 1.0), norms
    return V


def build_600cell_graph() -> nx.Graph:
    V = build_600cell_vertices()
    # Edges: pairs at dot product = phi/2 (i.e., chord distance 2 sin(pi/10))
    target = PHI / 2.0
    dots = V @ V.T
    G = nx.Graph()
    G.add_nodes_from(range(120))
    for i in range(120):
        for j in range(i + 1, 120):
            if abs(dots[i, j] - target) < 1e-9:
                G.add_edge(i, j)
    return G


def compute_spectrum(G: nx.Graph) -> np.ndarray:
    L = nx.laplacian_matrix(G).astype(np.float64).toarray()
    return np.sort(np.linalg.eigvalsh(L))


if __name__ == '__main__':
    G = build_600cell_graph()
    print(f"600-cell graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"avg degree {2*G.number_of_edges()/G.number_of_nodes():.2f}")
    assert G.number_of_nodes() == 120
    # Expected: 720 edges (each vertex has degree 12)
    assert G.number_of_edges() == 720, G.number_of_edges()
    assert nx.is_connected(G)
    evals = compute_spectrum(G)
    print(f"Smallest 10 Laplacian eigenvalues: {evals[:10]}")
    print(f"Spectral gap (smallest nonzero): lambda_1 = {evals[1]:.6f}")
    print(f"Normalized lambda_1 / <deg> = {evals[1] / 12.0:.6f}")
    print(f"lambda_1 * N^(2/3) = {evals[1] * (120 ** (2.0/3.0)):.4f}")
    # Save eigenvalues
    np.savetxt('reference_600cell_eigenvalues.csv', evals, delimiter=',', fmt='%.10f')
    print("Saved eigenvalues to reference_600cell_eigenvalues.csv")
