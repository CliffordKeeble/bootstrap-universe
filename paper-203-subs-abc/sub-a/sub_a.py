"""
Sub A: signature enumeration scope extension.

Method 1: bounded-coefficient search in the ℤ[φ]-linear span anticommuting
with G0 = Z⊗Γ_seed.

Step 1: Compute W = {M in V : {M, G0} = 0} via linear algebra over Q(sqrt 5).
Step 2: For each candidate M in W with ℤ[φ]-coefficients bounded by K,
  check M² = +/- sqrt(5)*I_4.
Step 3: Among candidates with M² = −√5·I_4 (spacelike), find 3-cliques
  of mutually anticommuting matrices.
"""

from __future__ import annotations

import itertools
import sympy as sp
from sympy import sqrt, Matrix, eye, zeros, kronecker_product, simplify

# Foundation (matches paper-203-algebra/subs.py)
phi = (1 + sqrt(5)) / 2
root5 = 2 * phi - 1
I2 = eye(2)
I4 = eye(4)
K = Matrix([[0, 1], [1, 0]])
J = Matrix([[0, -1], [1, 0]])
Z = K * J  # diag(1, -1)
Gamma_seed = Matrix([[0, 1], [root5, 0]])
Gamma_adj = Matrix([[0, -1], [root5, 0]])

# Gamma_0 = Z x Gamma_seed (the unique +sqrt(5)*I literal candidate from Sub 3)
G0 = kronecker_product(Z, Gamma_seed)

# Basis of 2x2 matrices over ℤ[φ]
BASIS_2X2 = {
    'I': I2, 'K': K, 'J': J, 'Z': Z,
    'Gseed': Gamma_seed, 'Gadj': Gamma_adj,
}

# Build the 36 tensor-product basis for V
V_BASIS = []
V_NAMES = []
for nA, A in BASIS_2X2.items():
    for nB, B in BASIS_2X2.items():
        V_BASIS.append(kronecker_product(A, B))
        V_NAMES.append(f"{nA}x{nB}")


def anticomm(M, N):
    return simplify(M * N + N * M)


def is_zero(M):
    return simplify(M) == sp.zeros(*M.shape)


def step1_compute_W() -> tuple[list[Matrix], list[str]]:
    """
    Find the subspace W of V consisting of matrices anticommuting with G0.

    For each basis element b_i in V, compute {b_i, G0}. The result is a 4x4
    matrix in the 36-dim ambient. Stack these as columns of a big linear
    system, find the kernel.

    Returns: basis of W as a list of Matrix.
    """
    print("Step 1: computing W = {M : {M, G0} = 0} via linear algebra...")
    # Express each {b_i, G0} as a 16-vector (4x4 entries flattened).
    # The kernel of the resulting 16x36 linear map (over Q(sqrt 5)) is W.
    rows = []
    for b in V_BASIS:
        ac = simplify(b * G0 + G0 * b)
        flat = [ac[i, j] for i in range(4) for j in range(4)]
        rows.append(flat)
    # Build a 36x16 matrix where rows are flattened anticommutators
    A_mat = Matrix(rows)  # 36 x 16
    # We want coefficients c_i such that sum_i c_i * {b_i, G0} = 0.
    # That's the (left) nullspace of A_mat, equivalently the nullspace
    # of A_mat^T (16 x 36) acting on c (36-vector).
    null = A_mat.T.nullspace()
    print(f"  W has dimension {len(null)} over Q(sqrt 5)")
    # Each null vector is a 36-dim coefficient vector — express the corresponding
    # element of V as a sum.
    W_basis = []
    W_basis_names = []
    for vec in null:
        # vec is 36-dim; construct M = sum_i vec[i] * V_BASIS[i]
        M = sp.zeros(4, 4)
        for i in range(36):
            M = M + simplify(vec[i]) * V_BASIS[i]
        M = simplify(M)
        # Check {M, G0} = 0
        if not is_zero(anticomm(M, G0)):
            print(f"  WARNING: null vec failed anticommutator check, skipping")
            continue
        W_basis.append(M)
        # Identify name as the dominant basis element if simple
        terms = []
        for i in range(36):
            v_i = simplify(vec[i])
            if v_i != 0:
                terms.append(f"{V_NAMES[i]}({v_i})")
        W_basis_names.append(" + ".join(terms[:3]) + ("..." if len(terms) > 3 else ""))
    return W_basis, W_basis_names


def matrices_with_pm_root5_squares_in_W(W_basis: list[Matrix], K: int = 1) -> list[tuple[Matrix, str, int]]:
    """
    Enumerate matrices in W of the form M = sum_i c_i * w_i with
    c_i in {-K, ..., K} * {1, phi}, check M² = +/- sqrt(5)*I_4.

    Returns list of (M, label, sign) where sign = +1 if M² = +√5·I, -1 if -√5·I.
    """
    n_W = len(W_basis)
    # Coefficient values: {(a, b) : a in -K..K, b in -K..K}
    coeff_range = list(range(-K, K + 1))
    coeff_pairs = [(a, b) for a in coeff_range for b in coeff_range]
    print(f"\nStep 2/3: enumerating M = sum c_i * w_i with c_i = a + b*phi, "
          f"a, b in {{{-K}, ..., {K}}}; |W| = {n_W}, total combos = "
          f"{len(coeff_pairs)**n_W}")

    if len(coeff_pairs) ** n_W > 10 ** 8:
        print(f"  WARNING: search space too large ({len(coeff_pairs)**n_W}); "
              f"may need to reduce K or n_W")
        return []

    out = []
    seen_keys = set()  # avoid duplicates up to overall sign
    count = 0
    for combo in itertools.product(coeff_pairs, repeat=n_W):
        # Skip the all-zero combo
        if all(p == (0, 0) for p in combo):
            continue
        count += 1
        # Build M
        M = sp.zeros(4, 4)
        coeffs = []
        for w, (a, b) in zip(W_basis, combo):
            c = a + b * phi
            coeffs.append((a, b))
            if c != 0:
                M = M + c * w
        M = simplify(M)
        # Deduplicate up to sign
        # canonical key: tuple of coefficients (a,b) and (-a,-b) → pick lex-smaller
        key_pos = tuple(coeffs)
        key_neg = tuple((-a, -b) for (a, b) in coeffs)
        canonical = min(key_pos, key_neg)
        if canonical in seen_keys:
            continue
        seen_keys.add(canonical)
        # Compute M^2
        Msq = simplify(M * M)
        # Check if M^2 is a scalar multiple of I_4 with that scalar = ±√5
        if Msq == root5 * I4:
            out.append((M, str(combo), +1))
        elif Msq == -root5 * I4:
            out.append((M, str(combo), -1))
    print(f"  enumerated {count} non-zero combos; found {len(out)} matrices with M² = +/- sqrt(5)*I_4")
    return out


def find_clique_with_signature(candidates: list[tuple[Matrix, str, int]],
                                 target_signature: tuple[int, int]) -> list | None:
    """
    target_signature = (p, q): p elements with sign +1, q with sign -1.
    For (1, 3) Minkowski: target = (1, 3). But G0 is fixed (+1); we need
    (0, 3) in the remaining 3 gammas.
    For (G1, G2, G3) all spacelike (sign = -1), signature (1, 3) overall.

    Returns list of three matrices forming a clique, or None.
    """
    # Filter to sign = -1 (spacelike) candidates
    spacelike = [(M, lbl) for M, lbl, sgn in candidates if sgn == -1]
    print(f"\nLooking for 3-clique of mutually anticommuting spacelike "
          f"matrices among {len(spacelike)} candidates...")

    if len(spacelike) < 3:
        return None

    # Build anticommutation graph
    n = len(spacelike)
    ac_graph = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if is_zero(anticomm(spacelike[i][0], spacelike[j][0])):
                ac_graph[i][j] = True
                ac_graph[j][i] = True
    # Find 3-cliques
    for combo in itertools.combinations(range(n), 3):
        if all(ac_graph[a][b] for a, b in itertools.combinations(combo, 2)):
            return [spacelike[i] for i in combo]
    return None


def main(K_bound: int = 1):
    print("=" * 70)
    print(f"Sub A — Method 1 — search bound K = {K_bound}")
    print("=" * 70)

    W, W_names = step1_compute_W()
    print(f"\nW basis (dim {len(W)}):")
    for i, (w, name) in enumerate(zip(W, W_names)):
        print(f"  w[{i}] ~ {name}")

    cands = matrices_with_pm_root5_squares_in_W(W, K=K_bound)
    if not cands:
        print(f"\nNo candidates found at K={K_bound}.")
        print("Verdict: NOT CONSTRUCTIBLE at K = {K_bound} (no +/- sqrt(5)*I candidates).")
        return None

    print(f"\nSignature distribution of M² values:")
    pos = sum(1 for _, _, s in cands if s == +1)
    neg = sum(1 for _, _, s in cands if s == -1)
    print(f"  +√5·I: {pos}")
    print(f"  -√5·I: {neg}")

    clique = find_clique_with_signature(cands, target_signature=(0, 3))
    if clique is None:
        print(f"\nVerdict at K = {K_bound}: NOT CONSTRUCTIBLE — no 3-clique of "
              f"mutually anticommuting spacelike matrices found.")
        return cands, None
    else:
        print(f"\nVerdict: CONSTRUCTIBLE at K = {K_bound}!")
        for i, (M, lbl) in enumerate(clique):
            print(f"  G^{i+1}: {lbl[:80]}...")
            print(f"  M^2 = {simplify(M*M).tolist()}")
        return cands, clique


if __name__ == '__main__':
    import sys
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    main(K_bound=K)
