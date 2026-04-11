# golden_dirac.py
# Golden Dirac operator — gamma matrices over Z[phi]
# 2I Universe Programme / Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 11 2026

import sympy as sp
from sympy import sqrt, Rational, Matrix, eye, zeros, kronecker_product

# ── Constants ─────────────────────────────────────────────────────────────────

phi = (1 + sqrt(5)) / 2      # golden ratio
psi = 1 - phi                 # golden conjugate = -1/phi
root5 = 2*phi - 1             # √5 ∈ Z[phi] — native, no import

# ── Seed matrices (2×2 over Z[phi]) ──────────────────────────────────────────
#
# Γ_seed  = [[0,  1], [√5,  0]]    Γ_seed²  = +√5·I₂
# Γ_adj   = [[0, -1], [√5,  0]]    Γ_adj²   = −√5·I₂
#
# {Γ_seed, Γ_adj} = 0              — they anticommute
#
# Γ_adj is Γ_seed with a single sign flip. The pair forms a 2D
# Clifford algebra over Z[phi] with norm √5.

Gamma_seed = Matrix([
    [0,       1],
    [root5,   0]
])

Gamma_adj = Matrix([
    [0,      -1],
    [root5,   0]
])

# ── Building blocks (2×2 drivers) ────────────────────────────────────────────
#
# K = σ₁  (swap),   K² = +I
# J ≈ iσ₂ (rotation), J² = -I
#
# {K, J} = 0 — they anticommute

K = Matrix([[0, 1], [1,  0]])    # K² =  I
J = Matrix([[0,-1], [1,  0]])    # J² = -I

# ── Three golden gamma matrices (4×4 over Z[phi]) ────────────────────────────
#
# G1 = K ⊗ Γ_seed    →  (G1)² = K²⊗Γ_seed² = +√5·I₄
# G2 = J ⊗ Γ_seed    →  (G2)² = J²⊗Γ_seed² = −√5·I₄
# G3 = I ⊗ Γ_adj     →  (G3)² = I²⊗Γ_adj²  = −√5·I₄
#
# Anticommutation:
#   {G1,G2} = {K,J}⊗(Γ_seed·Γ_seed) = 0⊗(...) = 0  ✓
#   {G1,G3} = (K·I+I·K)⊗? = K⊗{Γ_seed,Γ_adj} = K⊗0 = 0  ✓
#   {G2,G3} = (J·I+I·J)⊗? = J⊗{Γ_seed,Γ_adj} = J⊗0 = 0  ✓
#
# Signature: (+,−,−) — Minkowski, from the algebra, not imposed.
# No i anywhere. Every entry in Z[phi]. Discrete from the start.

G1 = kronecker_product(K, Gamma_seed)
G2 = kronecker_product(J, Gamma_seed)
G3 = kronecker_product(eye(2), Gamma_adj)

# ── Verification ──────────────────────────────────────────────────────────────

def anticommutator(A, B):
    return sp.simplify(A*B + B*A)

def verify_all(G1, G2, G3, root5):
    I4 = eye(4)
    results = {}

    results['(G1)² = +√5·I₄'] = sp.simplify(G1*G1 - root5*I4) == zeros(4)
    results['(G2)² = −√5·I₄'] = sp.simplify(G2*G2 + root5*I4) == zeros(4)
    results['{G1,G2} = 0']     = anticommutator(G1,G2) == zeros(4)
    results['(G3)² = −√5·I₄'] = sp.simplify(G3*G3 + root5*I4) == zeros(4)
    results['{G1,G3} = 0']     = anticommutator(G1,G3) == zeros(4)
    results['{G2,G3} = 0']     = anticommutator(G2,G3) == zeros(4)

    return results

print("=" * 60)
print("GOLDEN DIRAC — GAMMA MATRIX VERIFICATION")
print("=" * 60)

results = verify_all(G1, G2, G3, root5)

all_pass = True
for relation, passed in results.items():
    status = "PASS" if passed else "FAIL"
    if not passed:
        all_pass = False
    print(f"  {status}   {relation}")

print()

if all_pass:
    print("All relations satisfied. Three golden gammas confirmed over Z[phi].")
    print()
    print("-- G1 (K x Gamma_seed) --")
    sp.pprint(G1)
    print()
    print("-- G2 (J x Gamma_seed) --")
    sp.pprint(G2)
    print()
    print("-- G3 (I x Gamma_adj) --")
    sp.pprint(G3)
    print()
    print("Full golden Clifford relation:")
    print("  {G^u, G^v} = 2 * g^{uv} * sqrt(5) * I_4")
    print("  Signature: (+,-,-) from the algebra, not imposed.")
    print("  No i anywhere. Entries in Z[phi]. Discrete from the start.")
    print()
    print("Construction:")
    print("  Gamma_seed = [[0, 1], [sqrt(5), 0]]     seed^2 = +sqrt(5)*I")
    print("  Gamma_adj  = [[0,-1], [sqrt(5), 0]]     adj^2  = -sqrt(5)*I")
    print("  {seed, adj} = 0")
    print()
    print("  G1 = K x seed    K = [[0,1],[1,0]]      K^2 = +I")
    print("  G2 = J x seed    J = [[0,-1],[1,0]]     J^2 = -I")
    print("  G3 = I x adj     I = identity")
else:
    print("VERIFICATION FAILED.")
    for relation, passed in results.items():
        if not passed:
            print(f"  Failed: {relation}")

# ── Eigenvalue hunt ───────────────────────────────────────────────────────────
# Compute eigenvalues of G1 + G2 + G3 as proxy for Dirac operator
# These point toward the minimum Planck eigenvalue

if all_pass:
    print()
    print("-- Eigenvalue hunt (proxy Dirac operator G1+G2+G3) --")
    D_proxy = G1 + G2 + G3
    eigenvals = D_proxy.eigenvals()
    print("Eigenvalues:")
    for val, mult in eigenvals.items():
        print(f"  {sp.simplify(val)}  (multiplicity {mult})")
    print()
    print("Minimum nonzero |eigenvalue| is the Planck candidate.")
