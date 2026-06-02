"""
Paper 157 rigour tests — PRE-LOCK SCOPING (not the locked test).

Purpose: establish the representation-theory ground truth that the two
locks (Test A operator, Test B correction form) must be chosen against.
This script computes only standard, exact facts about A5 and the
augmentation ideal of its regular representation, plus the precise size
and sign of the 0.02% residual. It declares NO pass/null verdict — that
waits for the frozen lock per Pattern 75.

Run: python scoping_rep_theory.py
"""

import sympy as sp

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2          # golden ratio
psi = (1 - sqrt5) / 2          # conjugate

# ---------------------------------------------------------------------------
# A5 character table.
# Conjugacy classes (size):  1A(1)  2A(15)  3A(20)  5A(12)  5B(12)
# Irreps (dim):              1, 3, 3', 4, 5
# 5A character values for the two 3-dims involve phi, psi (golden).
# ---------------------------------------------------------------------------
class_names = ["1A", "2A", "3A", "5A", "5B"]
class_sizes = [1, 15, 20, 12, 12]

# rows = irreps, cols = classes (in class_names order)
irreps = {
    "1":  (1, [1, 1, 1, 1, 1]),
    "3":  (3, [3, -1, 0, phi, psi]),
    "3'": (3, [3, -1, 0, psi, phi]),
    "4":  (4, [4, 0, 1, -1, -1]),
    "5":  (5, [5, 1, -1, 0, 0]),
}

G = sum(class_sizes)
print(f"|A5| = {G}")
print(f"sum of squared irrep dims = {sum(d*d for d,_ in irreps.values())}  (must equal |A5|)")
assert sum(d*d for d, _ in irreps.values()) == G

# --- 1. Orthonormality of characters (verifies the table is correct) ------
def inner(chiA, chiB):
    return sp.Rational(1, G) * sum(
        sz * a * sp.conjugate(b) for sz, a, b in zip(class_sizes, chiA, chiB)
    )

print("\n--- character orthonormality (should be I) ---")
keys = list(irreps)
ok = True
for i in keys:
    row = []
    for j in keys:
        val = sp.simplify(inner(irreps[i][1], irreps[j][1]))
        row.append(val)
        expected = 1 if i == j else 0
        if val != expected:
            ok = False
    print(f"  <{i:>2}| {[str(x) for x in row]}")
print(f"orthonormal: {ok}")
assert ok

# --- 2. Regular-rep augmentation ideal decomposition ----------------------
# Regular rep contains each irrep V_i with multiplicity = dim(V_i).
# Augmentation ideal = regular rep minus the single trivial copy.
print("\n--- augmentation ideal of regular rep (dim |A5|-1 = 59) ---")
aug_total = 0
blocks = []
for name, (d, _) in irreps.items():
    mult = d                      # multiplicity in regular rep = dim
    if name == "1":
        mult -= 1                 # remove the one trivial copy
    block_dim = mult * d
    aug_total += block_dim
    if mult > 0:
        blocks.append((name, d, mult, block_dim))
    print(f"  irrep {name:>2}: dim {d}, multiplicity {mult:>1}, block dim {block_dim:>2}")
print(f"augmentation ideal dim = {aug_total}  (must equal 59)")
assert aug_total == 59
print(f"non-trivial isotypic blocks: {len(blocks)}  "
      f"-> dims { {b[0]: b[3] for b in blocks} }")

# --- 3. Are the blocks INEQUIVALENT? (Schur: invariant op is block-scalar)-
# An A5-invariant operator acts as a scalar on each *isotypic* block.
# The four non-trivial blocks {3, 3', 4, 5} are pairwise inequivalent
# irreps, so the four block-scalars are independent. Equipartition (one
# common eigenvalue across all 59 dims) therefore requires the operator
# to be scalar on the WHOLE space -- i.e. it carries NO spectral gap.
print("\n--- class-sum (central) operator eigenvalues per irrep ---")
print("    omega_C(V) = |C| * chi_V(C) / dim(V);  central ops act as this scalar")
print(f"    {'class':>6} | " + " | ".join(f"{b[0]:>8}" for b in blocks))
for ci, cname in enumerate(class_names):
    if cname == "1A":
        continue   # identity class sum is trivial (= |C|=1 -> scalar 1/dim, no gap)
    row = []
    for name, d, mult, bd in blocks:
        chi = irreps[name][1][ci]
        omega = sp.nsimplify(class_sizes[ci] * chi / d)
        row.append(omega)
    vals = [sp.simplify(x) for x in row]
    spread = "UNIFORM" if len(set(vals)) == 1 else "NON-uniform"
    print(f"    {cname:>6} | " + " | ".join(f"{str(v):>8}" for v in vals)
          + f"   <- {spread}")

print("""
Structural conclusion (Schur's lemma):
  The augmentation ideal splits into 4 pairwise-inequivalent isotypic
  blocks. Every A5-invariant operator is a scalar on each block, and the
  four scalars are independent. A uniform spectrum across all 59 modes
  (equipartition) requires all four scalars equal -> the operator is a
  multiple of the identity -> it has NO gap. A non-trivial spectral GAP
  is, by construction, a distinction between modes; equipartition is the
  ABSENCE of distinction. The two pull opposite ways.
  => 'spectral gap forces equipartition' is in structural tension; this
     is the landscape the Test A operator-lock must be chosen against.
""")

# --- 4. Test B scoping: precise residual size and SIGN --------------------
print("--- Test B scoping: the 0.02% residual ---")
L = 2 * sp.log(phi) / sqrt5
a_inv = 59 / L
codata = sp.Float("137.035999177", 15)
L_f = L.evalf(15)
a_inv_f = a_inv.evalf(15)
gap_abs = a_inv_f - codata
gap_frac = gap_abs / a_inv_f
print(f"  L(1,chi5)        = {L_f}")
print(f"  alpha^-1 = 59/L  = {a_inv_f}")
print(f"  CODATA 2022      = {codata}")
print(f"  residual (abs)   = {gap_abs}   ({'OVERSHOOT' if gap_abs>0 else 'undershoot'})")
print(f"  residual (frac)  = {gap_frac.evalf(4)}  -> need a NEGATIVE correction of this size")
print(f"  effective count  = alpha^-1_CODATA * L = {(codata*L_f).evalf(8)}  (vs 59)")

print("\n  scale menu of FORCED augmentation quantities (no fitted coeff):")
for label, expr in [
    ("1/119",          sp.Rational(1, 119)),
    ("1/239",          sp.Rational(1, 239)),
    ("59/119",         sp.Rational(59, 119)),
    ("59/239",         sp.Rational(59, 239)),
    ("1/119^2",        sp.Rational(1, 119**2)),
    ("1/239^2",        sp.Rational(1, 239**2)),
    ("1/(59*119)",     sp.Rational(1, 59*119)),
    ("1/(59*239)",     sp.Rational(1, 59*239)),
    ("L/239",          L / 239),
]:
    print(f"    {label:>12} = {expr.evalf(5)}")
print(f"\n  target fractional correction = {gap_frac.evalf(4)}")
print("  (NOTE: scale comparison only. No 'best form' is selected here -")
print("   the single correction form is locked in the pre-reg, then run once.)")
