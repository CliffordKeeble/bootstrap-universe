"""
Paper 157 rigour test v0.2 — the §2 FORCING GATE.

The gate is the test (brief v0.2 §2, §6). It asks, for each candidate
weighting rule W1/W2/W3: can the four block weights {w3, w3', w4, w5} be
stated from STRUCTURE ALONE (character table / dimensions / a forced
operator spectrum) WITHOUT reference to the target 58.9815, AND is the
rule forced (not four free knobs)?

If NO candidate passes -> HALT, "NULL BY CONSTRUCTION" (brief §2, §6,
the expected safe outcome). The gate is evaluated on FORCEDNESS and
SCALE; it does NOT compute N_eff against the ±0.003 window for any
candidate (that is §3, reached only if a candidate passes the gate).
Showing that ALL candidates fail is documenting the null, not a
multiple-comparison search for a winner (brief §3.6).

Run: python gate_evaluation.py
"""

import sympy as sp

sqrt5 = sp.sqrt(5)
phi = (1 + sqrt5) / 2
psi = (1 - sqrt5) / 2

# A5 character table (verified orthonormal in scoping_rep_theory.py)
# classes:        1A(1) 2A(15) 3A(20) 5A(12) 5B(12)
class_order = [1, 2, 3, 5, 5]
class_size = [1, 15, 20, 12, 12]
chi = {
    "3":  [3, -1, 0, phi, psi],
    "3'": [3, -1, 0, psi, phi],
    "4":  [4, 0, 1, -1, -1],
    "5":  [5, 1, -1, 0, 0],
}
block_dim = {"3": 9, "3'": 9, "4": 16, "5": 25}   # multiplicity*dim in reg rep
irrep_dim = {"3": 3, "3'": 3, "4": 4, "5": 5}

# Target context (NOT used inside any weighting rule; shown only to state
# the scale the gate must reach).
import mpmath as mp
mp.mp.dps = 30
L = 2 * mp.log((1 + mp.sqrt(5)) / 2) / mp.sqrt(5)
codata = mp.mpf("137.035999177")
N_uniform = 59
N_target = codata * L                      # = 58.9815
deficit = N_uniform - N_target             # = 0.0185
frac = deficit / N_uniform                 # = 3.13e-4
print("=== scale the gate must reach (context only) ===")
print(f"  uniform count      N = {N_uniform}")
print(f"  CODATA-match count N = {mp.nstr(N_target, 8)}")
print(f"  deficit            = {mp.nstr(deficit, 6)}   (fractional {mp.nstr(frac, 4)})")
print(f"  => a passing forced weighting must perturb uniform weights")
print(f"     (all w_i = 1) by ~3e-4, i.e. land at 1 +/- ~0.001.\n")

# ---------------------------------------------------------------------------
# CANDIDATE W1 — chi5 character weighting
# Map chi5 (Dirichlet char mod 5) onto blocks via element order:
#   chi5(ord): ord1->chi5(1)=+1, ord2->chi5(2)=-1, ord3->chi5(3)=-1,
#              ord5->chi5(5)=0.
# A natural forced contraction: w(block) ~ sum_C |C| chi_block(C) chi5(ord C).
# ---------------------------------------------------------------------------
chi5_of = {1: 1, 2: -1, 3: -1, 5: 0}       # Legendre symbol mod 5
print("=== CANDIDATE W1: chi5 character weighting ===")
w1 = {}
for b, row in chi.items():
    s = sum(class_size[i] * row[i] * chi5_of[class_order[i]] for i in range(5))
    w1[b] = sp.simplify(s)
print(f"  raw chi5-contraction per block: {dict((b, str(v)) for b, v in w1.items())}")
print("  -> values are O(10), one negative (block 4 = -16). To become")
print("     multiplicative weights ~1 they need a free additive shift +")
print("     scale: TWO free knobs not fixed by structure. And the spread")
print("     is O(1), not O(1e-3).")
print("  GATE: FAIL. Not forced (shift+scale free); wrong scale (O(1)).\n")

# ---------------------------------------------------------------------------
# CANDIDATE W2 — dimension / "second augmentation" weighting
# Idea: a forced quantity equal to the deficit 0.0185, or removing one
# more forced direction.
# ---------------------------------------------------------------------------
print("=== CANDIDATE W2: dimension / second-augmentation weighting ===")
print("  Removing one more whole mode: N_eff = 59 - 1 = 58  (residual ~1.7%,")
print("  overshoots the gap by ~55x and wrong side of 58.98). Integer mode")
print("  removal cannot reach 58.98.")
print("  Forced fractional deficit must equal 0.0185. Nearest clean forms:")
for label, val in [
    ("1/59",  sp.Rational(1, 59)),
    ("1/54",  sp.Rational(1, 54)),
    ("L/2",   None),
]:
    if val is not None:
        print(f"    {label:>6} = {float(val):.5f}")
print(f"    deficit  = {mp.nstr(deficit,5)}")
print("  1/59 = 0.01695 (off 9%); 1/54 = 0.01852 fits but 54 is not a")
print("  forced A5/augmentation quantity (not 59, 60, 11, 119, 239, ...).")
print("  Any '54' is reverse-engineered from the target = peeking.")
print("  GATE: FAIL. No forced dimensional quantity equals 0.0185 without")
print("        reference to 58.98; integer mode removal lands on 58.\n")

# ---------------------------------------------------------------------------
# CANDIDATE W3 — spectral-gap weighting (eigenvalue ratios)
# Weights = forced eigenvalue ratios of a NAMED operator on the 4 blocks.
# A class-sum (central) operator is the most natural icosahedral choice;
# show its per-block eigenvalues for the three rotation-axis classes.
# ---------------------------------------------------------------------------
print("=== CANDIDATE W3: spectral-gap (eigenvalue-ratio) weighting ===")
print("  Class-sum eigenvalue omega_C(V) = |C| chi_V(C)/dim(V) per block:")
axis_class = {"2A (15, edge 2-fold)": 1, "3A (20, face 3-fold)": 2,
              "5A (12, vertex 5-fold)": 3}
for cname, ci in axis_class.items():
    vals = {b: sp.simplify(sp.nsimplify(class_size[ci] * chi[b][ci] / irrep_dim[b]))
            for b in chi}
    print(f"    {cname:<24}: " + ", ".join(f"{b}={vals[b]}" for b in chi))
print("  Eigenvalues are O(1) integers/algebraics, include 0 and negatives,")
print("  so they are not directly multiplicative weights; forming weights")
print("  needs a free shift+scale (unforced), and the generating-set choice")
print("  (which class?) is itself unforced -> different classes give")
print("  different weights. Spread is O(1), not O(1e-3).")
print("  GATE: FAIL. Operator/generating-set choice free; wrong scale.\n")

# ---------------------------------------------------------------------------
print("=" * 64)
print("GATE VERDICT: NULL BY CONSTRUCTION.")
print("No candidate (W1/W2/W3) yields a forced, peek-free weighting rule")
print("at the required ~3e-4 scale. Every rep-theoretic weighting from")
print("{characters, dimensions, operator spectra} produces O(1) deviations")
print("from uniform; the residual needs O(1e-3). Bridging 3 orders of")
print("magnitude requires a fitted small parameter = the four-knobs trap")
print("the brief (§1, §4) forbids.")
print("=> Per brief §2/§6: HALT. The 0.031% residual is the honest")
print("   precision ceiling of a convergence result, not a forced block-")
print("   weight correction. 157 stays OBSERVED at corrected 0.031%;")
print("   convergence argument (Pattern 9) stands; no numerology acquired.")
print("=" * 64)
