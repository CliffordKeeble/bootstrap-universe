"""
Sub A v2: properly count W as a subspace of M_4(Q(sqrt 5)).

Treat each 4x4 matrix as a 16-vector. The anticommutator with G0 is a
linear endomorphism on M_4. W = ker(M -> {M, G0}) is its kernel.
True dimension of W is at most 16 (typically 8 for G0 with G0^2 = scalar I).
"""

from __future__ import annotations

import itertools
import sympy as sp
from sympy import sqrt, Matrix, eye, zeros, kronecker_product, simplify

phi = (1 + sqrt(5)) / 2
root5 = 2 * phi - 1
I2 = eye(2)
I4 = eye(4)
K = Matrix([[0, 1], [1, 0]])
J = Matrix([[0, -1], [1, 0]])
Z = K * J
Gamma_seed = Matrix([[0, 1], [root5, 0]])
Gamma_adj = Matrix([[0, -1], [root5, 0]])

G0 = kronecker_product(Z, Gamma_seed)


def matrix_to_vec(M):
    """Flatten 4x4 to 16-vector (column major)."""
    return Matrix([M[i, j] for i in range(4) for j in range(4)])


def vec_to_matrix(v):
    """Inflate 16-vector to 4x4."""
    return Matrix(4, 4, list(v))


def anticomm_endomorphism(G0):
    """Build the 16x16 matrix representing M -> {M, G0} on Mat(4)."""
    # Express in standard basis E_{ij} for 4x4 matrices.
    # E_{ij} has 1 at (i, j), 0 elsewhere. Index = 4*i + j.
    cols = []
    for i in range(4):
        for j in range(4):
            E = sp.zeros(4, 4)
            E[i, j] = 1
            AC = simplify(E * G0 + G0 * E)
            cols.append(matrix_to_vec(AC))
    A = Matrix.hstack(*cols)
    return A


def find_W_basis():
    """Find an M_4-basis of W = {M : {M, G0} = 0}."""
    A = anticomm_endomorphism(G0)
    null = A.nullspace()
    print(f"W has true dim {len(null)} over Q(sqrt 5)")
    W_matrices = [vec_to_matrix(v) for v in null]
    return W_matrices


def main_search(K_bound: int = 1):
    print("=" * 70)
    print(f"Sub A v2 - Method 1 - search bound K = {K_bound}")
    print("=" * 70)
    W = find_W_basis()
    print(f"\nW basis matrices ({len(W)} elements):")
    for i, M in enumerate(W):
        # Simplify and print top-left block for sanity
        print(f"  w[{i}]: shape {M.shape}, top-left = {M[:2, :2].tolist()}")

    # Verify each w in W satisfies {w, G0} = 0
    for i, w in enumerate(W):
        if not simplify(w * G0 + G0 * w) == sp.zeros(4, 4):
            print(f"  WARNING: w[{i}] fails anticommutator")
            return None

    n_W = len(W)
    print(f"\nSearch space size: 9^{n_W} = {9**n_W} for K=1, a+b*phi with a,b in {{-1,0,1}}")
    if 9**n_W > 10**8:
        print(f"  Brute force at K={K_bound} infeasible.")
        if K_bound == 1:
            print(f"  Falling back to single-coefficient search: M = c * w_i for some i, c in Z[phi] bounded.")
            # Test each basis element scaled
            found_pos = []
            found_neg = []
            for i, w in enumerate(W):
                # Try c = 1, phi, 1+phi, etc. up to K=1
                for a in range(-K_bound, K_bound + 1):
                    for b in range(-K_bound, K_bound + 1):
                        if a == 0 and b == 0:
                            continue
                        c = a + b * phi
                        M = simplify(c * w)
                        Msq = simplify(M * M)
                        if Msq == root5 * I4:
                            found_pos.append((M, f"({a}+{b}*phi)*w[{i}]"))
                        elif Msq == -root5 * I4:
                            found_neg.append((M, f"({a}+{b}*phi)*w[{i}]"))
            print(f"  Single-element candidates: +sqrt(5)*I count = {len(found_pos)}, "
                  f"-sqrt(5)*I count = {len(found_neg)}")
            for M, lbl in found_pos[:5]:
                print(f"    +sqrt(5): {lbl}")
            for M, lbl in found_neg[:5]:
                print(f"    -sqrt(5): {lbl}")
            return W, found_pos, found_neg

    # Full brute force if space is small enough
    coeff_pairs = [(a, b) for a in range(-K_bound, K_bound + 1) for b in range(-K_bound, K_bound + 1)]
    print(f"  Brute force search: {len(coeff_pairs)**n_W} combinations")
    cands_pos = []
    cands_neg = []
    for combo in itertools.product(coeff_pairs, repeat=n_W):
        if all(p == (0, 0) for p in combo):
            continue
        M = sp.zeros(4, 4)
        for w, (a, b) in zip(W, combo):
            c = a + b * phi
            if c != 0:
                M = M + c * w
        M = simplify(M)
        Msq = simplify(M * M)
        if Msq == root5 * I4:
            cands_pos.append((M, str(combo)))
        elif Msq == -root5 * I4:
            cands_neg.append((M, str(combo)))
    return W, cands_pos, cands_neg


if __name__ == '__main__':
    import sys
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    result = main_search(K_bound=K)
    if result is not None:
        W, pos, neg = result
        print(f"\nTotal +sqrt(5)*I candidates: {len(pos)}")
        print(f"Total -sqrt(5)*I candidates: {len(neg)}")
