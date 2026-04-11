# eigenvalue_analysis.py
# Golden Dirac — eigenvalue and chirality analysis
# 2I Universe Programme / Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 11 2026

import sympy as sp
from sympy import sqrt, Rational, Matrix, eye, zeros, kronecker_product, Abs, Min

# ── Constants ─────────────────────────────────────────────────────────────────

phi = (1 + sqrt(5)) / 2
psi = 1 - phi
root5 = 2*phi - 1

# ── Seed matrices ─────────────────────────────────────────────────────────────

Gamma_seed = Matrix([
    [0,       1],
    [root5,   0]
])

Gamma_adj = Matrix([
    [0,      -1],
    [root5,   0]
])

# ── Building blocks ───────────────────────────────────────────────────────────

K = Matrix([[0, 1], [1,  0]])
J = Matrix([[0,-1], [1,  0]])

# ── Gamma matrices ────────────────────────────────────────────────────────────

G1 = kronecker_product(K, Gamma_seed)
G2 = kronecker_product(J, Gamma_seed)
G3 = kronecker_product(eye(2), Gamma_adj)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Individual gamma matrix eigenvalues
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("1. INDIVIDUAL GAMMA MATRIX EIGENVALUES")
print("=" * 60)

for name, M in [("G1", G1), ("G2", G2), ("G3", G3)]:
    evals = M.eigenvals()
    print(f"\n  {name} eigenvalues:")
    for val, mult in evals.items():
        v = sp.simplify(val)
        try:
            nv = complex(v.evalf())
            if abs(nv.imag) < 1e-15:
                print(f"    {v}  (mult {mult})  = {nv.real:.8f}")
            else:
                print(f"    {v}  (mult {mult})  = {nv.real:.8f} + {nv.imag:.8f}i")
        except:
            print(f"    {v}  (mult {mult})")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Proxy Dirac operator spectrum
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("2. PROXY DIRAC OPERATOR D = G1 + G2 + G3")
print("=" * 60)

D_proxy = G1 + G2 + G3

print("\n  D_proxy matrix:")
sp.pprint(D_proxy, wrap_line=False)

evals_D = D_proxy.eigenvals()
print("\n  Eigenvalues:")
for val, mult in evals_D.items():
    v = sp.simplify(val)
    try:
        nv = complex(v.evalf())
        if abs(nv.imag) < 1e-15:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f}")
        else:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f} + {nv.imag:.8f}i")
    except:
        print(f"    {v}  (mult {mult})")

# Minimum nonzero |eigenvalue|
print("\n  Minimum nonzero |eigenvalue|:")
nonzero = [sp.Abs(v) for v in evals_D.keys()
           if sp.simplify(v) != 0]
if nonzero:
    # Evaluate numerically to find minimum
    num_vals = [(v, float(v.evalf())) for v in nonzero]
    min_val = min(num_vals, key=lambda x: x[1])
    print(f"    |lambda_min| = {sp.simplify(min_val[0])}")
    print(f"    Numerical    = {min_val[1]:.8f}")
    print(f"    5^(1/4)      = {float(sp.Rational(5)**sp.Rational(1,4)):.8f}")
    print(f"    Match: {abs(min_val[1] - float(sp.Rational(5)**sp.Rational(1,4))) < 1e-10}")
else:
    print("    All eigenvalues are zero!")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Chirality operator
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("3. CHIRALITY OPERATOR G5 = G1*G2*G3")
print("=" * 60)

Gamma5 = sp.simplify(G1 * G2 * G3)

print("\n  G5 matrix:")
sp.pprint(Gamma5, wrap_line=False)

print("\n  G5 eigenvalues:")
evals_G5 = Gamma5.eigenvals()
for val, mult in evals_G5.items():
    v = sp.simplify(val)
    try:
        nv = complex(v.evalf())
        if abs(nv.imag) < 1e-15:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f}")
        else:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f} + {nv.imag:.8f}i")
    except:
        print(f"    {v}  (mult {mult})")

G5_sq = sp.simplify(Gamma5 * Gamma5)
print("\n  (G5)^2:")
sp.pprint(G5_sq, wrap_line=False)

# Check if (G5)^2 = +/- sqrt(5) * I
I4 = eye(4)
if G5_sq == root5 * I4:
    print("\n  (G5)^2 = +sqrt(5) * I_4   *** GOLDEN CHIRALITY ***")
elif G5_sq == -root5 * I4:
    print("\n  (G5)^2 = -sqrt(5) * I_4   *** GOLDEN CHIRALITY ***")
else:
    print(f"\n  (G5)^2 is not +/-sqrt(5)*I_4")
    # Check what scalar multiple of I it might be
    if G5_sq[0,1] == 0 and G5_sq[0,2] == 0 and G5_sq[0,3] == 0:
        diag = sp.simplify(G5_sq[0,0])
        print(f"  Diagonal entries = {diag}")
        print(f"  Numerical = {float(diag.evalf()):.8f}")
        print(f"  sqrt(5)   = {float(root5.evalf()):.8f}")
        print(f"  5         = 5.0")
        print(f"  5*sqrt(5) = {float((5*root5).evalf()):.8f}")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Seed vs Adjoint chirality
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("4. SEED vs ADJOINT CHIRALITY")
print("=" * 60)

commutator     = sp.simplify(Gamma_seed * Gamma_adj - Gamma_adj * Gamma_seed)
anticommutator = sp.simplify(Gamma_seed * Gamma_adj + Gamma_adj * Gamma_seed)

print("\n  [Gamma_seed, Gamma_adj] (commutator):")
sp.pprint(commutator, wrap_line=False)

print("\n  {Gamma_seed, Gamma_adj} (anticommutator):")
sp.pprint(anticommutator, wrap_line=False)

print("\n  Gamma_seed eigenvalues:")
evals_seed = Gamma_seed.eigenvals()
for val, mult in evals_seed.items():
    v = sp.simplify(val)
    try:
        nv = complex(v.evalf())
        if abs(nv.imag) < 1e-15:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f}")
        else:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f} + {nv.imag:.8f}i")
    except:
        print(f"    {v}  (mult {mult})")

print("\n  Gamma_adj eigenvalues:")
evals_adj = Gamma_adj.eigenvals()
for val, mult in evals_adj.items():
    v = sp.simplify(val)
    try:
        nv = complex(v.evalf())
        if abs(nv.imag) < 1e-15:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f}")
        else:
            print(f"    {v}  (mult {mult})  = {nv.real:.8f} + {nv.imag:.8f}i")
    except:
        print(f"    {v}  (mult {mult})")

# Check the big result
seed_evals_list = list(evals_seed.keys())
adj_evals_list  = list(evals_adj.keys())

print("\n  === THE BIG QUESTION ===")
print(f"  Seed eigenvalues:  {[sp.simplify(v) for v in seed_evals_list]}")
print(f"  Adj eigenvalues:   {[sp.simplify(v) for v in adj_evals_list]}")

# Check if seed eigenvalues are real +/- 5^(1/4)
# and adj eigenvalues are imaginary +/- i*5^(1/4)
seed_numerical = [complex(v.evalf()) for v in seed_evals_list]
adj_numerical  = [complex(v.evalf()) for v in adj_evals_list]

print(f"\n  Seed numerical: {seed_numerical}")
print(f"  Adj numerical:  {adj_numerical}")

seed_real = all(abs(v.imag) < 1e-10 for v in seed_numerical)
adj_imag  = all(abs(v.real) < 1e-10 for v in adj_numerical)

print(f"\n  Seed eigenvalues purely real?      {seed_real}")
print(f"  Adj eigenvalues purely imaginary?  {adj_imag}")

if seed_real and adj_imag:
    print("\n  *** i EMERGES FROM THE GOLDEN SPLIT ***")
    print("  Gamma_seed has real eigenvalues +/- 5^(1/4)")
    print("  Gamma_adj has imaginary eigenvalues +/- i*5^(1/4)")
    print("  Complex numbers are a CONSEQUENCE, not an assumption.")
elif seed_real and not adj_imag:
    print("\n  Seed is real, adj is NOT purely imaginary.")
    print("  Chirality exists but i does not emerge cleanly.")
elif not seed_real:
    print("\n  Seed eigenvalues are not purely real.")
    print("  The chiral structure may be more subtle.")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Product Gamma_seed * Gamma_adj
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("5. PRODUCT Gamma_seed * Gamma_adj")
print("=" * 60)

product = sp.simplify(Gamma_seed * Gamma_adj)
print("\n  Gamma_seed * Gamma_adj:")
sp.pprint(product, wrap_line=False)
print(f"\n  = {sp.simplify(product[0,0])} * I_2" if product[0,1]==0 and product[1,0]==0 else "")
print(f"  Eigenvalues: {sp.simplify(product.eigenvals())}")

product_rev = sp.simplify(Gamma_adj * Gamma_seed)
print("\n  Gamma_adj * Gamma_seed:")
sp.pprint(product_rev, wrap_line=False)
print(f"\n  = {sp.simplify(product_rev[0,0])} * I_2" if product_rev[0,1]==0 and product_rev[1,0]==0 else "")

print()
print("=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
