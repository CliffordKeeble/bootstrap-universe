#!/usr/bin/env python3
"""
Paper 48 §5 Audit — reproduction script.

Pre-registered Pattern 75 audit on the m_e c^2 icosahedral relation
    M = m_e c^2 / eV  ~=  2 (alpha^-1 - 1) (alpha^-1)^2 / (V - chi)

Brief: paper_48_s5_audit_brief_v0_1.md (committed before this script ran).
Single run. Locked specifications. Seed 48. No post-hoc tuning.

Part A  — arithmetic verification (A1) and algebraic-equivalence check (A2).
Part B  — Pattern 75 null distribution: search-space match-density (B1),
          with optional tighter audit (B2) iff B1 returns >= 1% (high branch).

Run:  python3 paper_48_s5_audit.py
Outputs: prints to stdout and writes paper_48_s5_audit_results.md
"""

import math
import os

import numpy as np

try:
    from scipy.stats import beta as _beta
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

from mpmath import mp, mpf, fabs, nstr

mp.dps = 40  # plenty for 12 significant figures in Part A

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(HERE, "paper_48_s5_audit_results.md")

# --------------------------------------------------------------------------
# Locked constants (CODATA 2018, per brief)
# --------------------------------------------------------------------------
ALPHA_INV = mpf("137.035999084")          # CODATA 2018 inverse fine-structure
V, CHI = 12, 2                             # icosahedral primitives
VMC = V - CHI                              # = 10
M = mpf("510998.95000")                    # m_e c^2 in eV (0.51099895000 MeV)
RINF_HC = mpf("13.605693122994")           # Rydberg energy R_inf h c in eV


def sig(x, n=12):
    """Format an mpf to n significant figures."""
    return nstr(x, n)


# ==========================================================================
# PART A1 — arithmetic verification
# ==========================================================================
def part_a1():
    rhs = 2 * (ALPHA_INV - 1) * ALPHA_INV**2 / VMC
    abs_res = fabs(M - rhs)
    rel_res = (M - rhs) / M            # signed
    rel_abs = fabs(rel_res)
    ppm = rel_res * mpf(10) ** 6
    ppb = rel_res * mpf(10) ** 9
    ppt = rel_res * mpf(10) ** 12
    pct = rel_res * 100
    return {
        "rhs": rhs,
        "M": M,
        "abs_res": abs_res,
        "rel_res": rel_res,
        "rel_abs": rel_abs,
        "ppm": ppm,
        "ppb": ppb,
        "ppt": ppt,
        "pct": pct,
    }


# ==========================================================================
# PART A2 — algebraic-equivalence check
# ==========================================================================
def part_a2(a1):
    e0_v3 = (ALPHA_INV - 1) / VMC                 # §3.2 fossil, dimensionless number
    res_32 = fabs(e0_v3 - RINF_HC) / RINF_HC      # §3.2 residual
    res_5 = a1["rel_abs"]                          # §5 residual (from A1)
    # If the algebraic equivalence holds, res_32 and res_5 should be equal up to
    # CODATA internal consistency. Report their difference in ppb of res_5.
    diff = fabs(res_32 - res_5)
    diff_ppb_of_res = (diff / res_5) * mpf(10) ** 9 if res_5 != 0 else mpf("nan")
    # Cross-check the exact physics identity M = 2 * R_inf_hc * (alpha^-1)^2.
    m_from_rinf = 2 * RINF_HC * ALPHA_INV**2
    codata_consistency_ppm = fabs(M - m_from_rinf) / M * mpf(10) ** 6
    equivalence_holds = diff_ppb_of_res < 10  # brief: agree within ~10 ppb
    return {
        "e0_v3": e0_v3,
        "res_32": res_32,
        "res_5": res_5,
        "diff": diff,
        "diff_ppb_of_res": diff_ppb_of_res,
        "m_from_rinf": m_from_rinf,
        "codata_consistency_ppm": codata_consistency_ppm,
        "equivalence_holds": equivalence_holds,
    }


# ==========================================================================
# PART B — Pattern 75 null distribution
# ==========================================================================
# Locked search space (Paper 127 §2 template):
#   S = {e, pi, sqrt2, sqrt3, sqrt5} U {1..30}            |S| = 35
#   1..4 terms; each term a product of 1..3 distinct primitives from S,
#   each raised to an integer exponent in {-5..-1, 1..5}; terms signed +/-.
ALLOWED_EXP = np.array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5], dtype=np.int64)

PRIMS = np.array(
    [math.e, math.pi, math.sqrt(2.0), math.sqrt(3.0), math.sqrt(5.0)]
    + [float(i) for i in range(1, 31)],
    dtype=np.float64,
)
assert PRIMS.size == 35, PRIMS.size

M_F = float(M)  # 510998.95


def sample_total(rng, B):
    """Vectorised draw of B expression values from the locked search space."""
    n_terms = rng.integers(1, 5, B)               # 1..4
    total = np.zeros(B, dtype=np.float64)
    for slot in range(4):
        active = slot < n_terms
        k = rng.integers(1, 4, B)                  # 1..3 primitives in this term
        # three distinct primitive indices in [0, 35)
        i1 = rng.integers(0, 35, B)
        i2 = rng.integers(0, 34, B)
        i2 += (i2 >= i1)
        i3 = rng.integers(0, 33, B)
        lo_i = np.minimum(i1, i2)
        hi_i = np.maximum(i1, i2)
        i3 += (i3 >= lo_i)
        i3 += (i3 >= hi_i)
        e1 = rng.choice(ALLOWED_EXP, B)
        e2 = rng.choice(ALLOWED_EXP, B)
        e3 = rng.choice(ALLOWED_EXP, B)
        # primitives beyond k contribute factor 1 (exponent 0)
        e2 = np.where(k >= 2, e2, 0)
        e3 = np.where(k >= 3, e3, 0)
        val = PRIMS[i1] ** e1 * PRIMS[i2] ** e2 * PRIMS[i3] ** e3
        sign = rng.choice(np.array([-1.0, 1.0]), B)
        term = sign * val
        total += np.where(active, term, 0.0)
    return total


def run_null(seed, N, tol, chunk=200_000):
    rng = np.random.default_rng(seed)
    lo, hi = M_F - tol, M_F + tol
    matches = 0
    drawn = 0
    while drawn < N:
        B = min(chunk, N - drawn)
        tot = sample_total(rng, B)
        matches += int(np.count_nonzero((tot >= lo) & (tot <= hi)))
        drawn += B
    return matches, drawn, lo, hi


def upper_bound_95(matches, N):
    """One-sided 95% upper bound on the match fraction."""
    if matches == 0:
        return 3.0 / N, "rule of three (3/N)"
    if _HAVE_SCIPY:
        # Clopper-Pearson upper limit at 95% (one-sided)
        ub = _beta.ppf(0.95, matches + 1, N - matches)
        return float(ub), "Clopper-Pearson 95% one-sided"
    # Normal approximation fallback
    p = matches / N
    ub = p + 1.645 * math.sqrt(p * (1 - p) / N)
    return ub, "normal approx (scipy unavailable)"


# ==========================================================================
# Driver + results document
# ==========================================================================
def main():
    a1 = part_a1()
    a2 = part_a2(a1)

    print("=" * 70)
    print("PART A1 — arithmetic verification")
    print("=" * 70)
    print(f"  RHS  = {sig(a1['rhs'])} eV")
    print(f"  M    = {sig(a1['M'])} eV")
    print(f"  |M-RHS| = {sig(a1['abs_res'])} eV")
    print(f"  rel residual = {sig(a1['ppm'],6)} ppm")
    print(f"               = {sig(a1['ppb'],6)} ppb")
    print(f"               = {sig(a1['ppt'],6)} ppt")
    print(f"               = {sig(a1['pct'],2)} %")

    print("\n" + "=" * 70)
    print("PART A2 — algebraic-equivalence check")
    print("=" * 70)
    print(f"  E0_v3 = (alpha^-1 - 1)/(V-chi) = {sig(a2['e0_v3'])}")
    print(f"  §3.2 residual |E0_v3 - Rinf_hc|/Rinf_hc = {sig(a2['res_32']*10**6,6)} ppm")
    print(f"  §5   residual |M - RHS|/M               = {sig(a2['res_5']*10**6,6)} ppm")
    print(f"  difference between residuals = {sig(a2['diff_ppb_of_res'],4)} ppb (of residual)")
    print(f"  equivalence holds (<10 ppb)? {a2['equivalence_holds']}")
    print(f"  [cross-check] M vs 2*Rinf_hc*(alpha^-1)^2: {sig(a2['codata_consistency_ppm'],3)} ppm")

    # ---- Part B1 ----
    SEED = 48
    N_B1 = 1_000_000
    TOL_B1 = 51.1
    print("\n" + "=" * 70)
    print("PART B1 — Pattern 75 null (N=1e6, tol=100 ppm=51.1 eV, seed=48)")
    print("=" * 70)
    m_b1, n_b1, lo_b1, hi_b1 = run_null(SEED, N_B1, TOL_B1)
    frac_b1 = m_b1 / n_b1
    ub_b1, ub_b1_method = upper_bound_95(m_b1, n_b1)
    print(f"  window = [{lo_b1:.2f}, {hi_b1:.2f}] eV")
    print(f"  N drawn = {n_b1}")
    print(f"  matches = {m_b1}")
    print(f"  match fraction = {frac_b1:.3e}")
    print(f"  95% upper bound = {ub_b1:.3e}  [{ub_b1_method}]")
    print(f"  pre-registered threshold = 1.0e-2")
    b1_sparse = frac_b1 < 0.01
    print(f"  match-rate < 1% threshold? {b1_sparse}")

    # ---- Part B2 (only if B1 is in the high branch, >= 1%) ----
    b2 = None
    if frac_b1 >= 0.01:
        N_B2 = 10_000_000
        TOL_B2 = 5.11
        print("\n" + "=" * 70)
        print("PART B2 — tighter audit (N=1e7, tol=10 ppm=5.11 eV, seed=48)")
        print("=" * 70)
        m_b2, n_b2, lo_b2, hi_b2 = run_null(SEED, N_B2, TOL_B2)
        frac_b2 = m_b2 / n_b2
        ub_b2, ub_b2_method = upper_bound_95(m_b2, n_b2)
        print(f"  matches = {m_b2}, fraction = {frac_b2:.3e}, 95% UB = {ub_b2:.3e}")
        b2 = dict(N=n_b2, matches=m_b2, frac=frac_b2, ub=ub_b2,
                  ub_method=ub_b2_method, tol=TOL_B2, lo=lo_b2, hi=hi_b2)
    else:
        print("\n[Part B2 skipped — B1 cleanly < 1% (not the high/inconclusive branch).]")

    write_results(a1, a2,
                  dict(seed=SEED, N=n_b1, matches=m_b1, frac=frac_b1, ub=ub_b1,
                       ub_method=ub_b1_method, tol=TOL_B1, lo=lo_b1, hi=hi_b1,
                       sparse=b1_sparse),
                  b2)
    print(f"\nResults written to {RESULTS_PATH}")


def write_results(a1, a2, b1, b2):
    L = []
    w = L.append
    w("# Paper 48 §5 Audit — Results")
    w("")
    w("Pre-registered Pattern 75 audit on the m_e c² icosahedral relation.")
    w("Brief: `paper_48_s5_audit_brief_v0_1.md` (committed before this script ran).")
    w("Generated by `paper_48_s5_audit.py`. Single run, seed 48, no post-hoc tuning.")
    w("")
    w("CODATA 2018 inputs: α⁻¹ = 137.035999084, V−χ = 10, "
      "M = m_e c²/eV = 510998.95000 eV, R∞hc = 13.605693122994 eV.")
    w("")
    w("---")
    w("")
    w("## Part A1 — arithmetic verification")
    w("")
    w("Relation: `RHS = 2 (α⁻¹ − 1) (α⁻¹)² / (V − χ)`")
    w("")
    w("| Quantity | Value |")
    w("|---|---|")
    w(f"| RHS (12 s.f.) | {sig(a1['rhs'])} eV |")
    w(f"| M = m_e c²/eV | {sig(a1['M'])} eV |")
    w(f"| Absolute residual \\|M − RHS\\| | {sig(a1['abs_res'])} eV |")
    w(f"| Relative residual (M − RHS)/M | {sig(a1['ppm'],6)} ppm |")
    w(f"| | {sig(a1['ppb'],6)} ppb |")
    w(f"| | {sig(a1['ppt'],6)} ppt |")
    w(f"| | {sig(a1['pct'],2)} % |")
    w("")
    w(f"**The §5 relation matches M to {sig(a1['ppm'],4)} ppm "
      f"({sig(a1['pct'],2)} %), not 0.0008%.** "
      "The v4.1 figure of 0.0008% (≈8 ppm) is refuted: it is low by ~20×.")
    w("")
    w("## Part A2 — algebraic-equivalence check")
    w("")
    w("§3.2 fossil: `E₀_v3 = (α⁻¹ − 1)/(V − χ)` compared to R∞hc (eV).")
    w("§5 relation residual from A1. The two are algebraically linked by "
      "`R∞hc = m_e c² α² / 2`, i.e. `M = 2 R∞hc (α⁻¹)²` and "
      "`RHS = 2 E₀_v3 (α⁻¹)²`, so the factor 2(α⁻¹)² cancels and the two "
      "relative residuals must be equal.")
    w("")
    w("| Quantity | Value |")
    w("|---|---|")
    w(f"| E₀_v3 | {sig(a2['e0_v3'])} |")
    w(f"| §3.2 residual \\|E₀_v3 − R∞hc\\|/R∞hc | {sig(a2['res_32']*10**6,6)} ppm |")
    w(f"| §5 residual \\|M − RHS\\|/M | {sig(a2['res_5']*10**6,6)} ppm |")
    w(f"| Difference between residuals | {sig(a2['diff_ppb_of_res'],4)} ppb (of residual) |")
    w(f"| Cross-check: M vs 2 R∞hc (α⁻¹)² | {sig(a2['codata_consistency_ppm'],3)} ppm |")
    w("")
    verdict = ("**CONFIRMED**" if a2["equivalence_holds"] else "**REFUTED**")
    w(f"Algebraic-equivalence claim: {verdict} "
      f"(residuals agree to {sig(a2['diff_ppb_of_res'],3)} ppb; brief threshold ~10 ppb). "
      "The §5 relation is the §3.2 fossil observation scaled by the exact "
      "Rydberg–electron-mass identity; it carries the same precision.")
    w("")
    w("## Part B1 — Pattern 75 null distribution")
    w("")
    w("Locked search space (Paper 127 §2 template): "
      "S = {e, π, √2, √3, √5} ∪ {1..30}, |S| = 35; 1–4 terms; each term a "
      "product of 1–3 distinct primitives raised to integer exponents in "
      "{−5..−1, 1..5}; signed ±. α⁻¹ excluded (treated as external constant).")
    w("")
    w("| Parameter | Value |")
    w("|---|---|")
    w(f"| Target M | {M_F:.2f} eV |")
    w(f"| Tolerance | {b1['tol']:.1f} eV (100 ppm) |")
    w(f"| Window | [{b1['lo']:.2f}, {b1['hi']:.2f}] eV |")
    w(f"| Seed | {b1['seed']} |")
    w(f"| N samples drawn | {b1['N']:,} |")
    w(f"| Matches within tolerance | {b1['matches']:,} |")
    w(f"| Match fraction | {b1['frac']:.3e} |")
    w(f"| 95% upper bound | {b1['ub']:.3e} ({b1['ub_method']}) |")
    w(f"| Pre-registered threshold | 1.0×10⁻² (1%) |")
    w("")
    if b1["sparse"]:
        w(f"**Verdict: match-rate ≪ 1% (fraction {b1['frac']:.2e}, 95% UB "
          f"{b1['ub']:.2e}). The §5 OBSERVED status earns non-trivial "
          "standing** — M lies in a sparse region of the icosahedral-primitive "
          "search space; the construction is not generic coincidence.")
    else:
        w(f"**Verdict: match-rate ≥ 1% (fraction {b1['frac']:.2e}). The §5 "
          "observation is consistent with generic match-density in this "
          "primitive space; OBSERVED status reflects coincidence rather than "
          "structural connection.** (See Part B2 for tighter discrimination.)")
    w("")
    if b2 is not None:
        w("## Part B2 — tighter audit (10 ppm)")
        w("")
        w("| Parameter | Value |")
        w("|---|---|")
        w(f"| Tolerance | {b2['tol']:.2f} eV (10 ppm) |")
        w(f"| Window | [{b2['lo']:.2f}, {b2['hi']:.2f}] eV |")
        w(f"| N samples drawn | {b2['N']:,} |")
        w(f"| Matches | {b2['matches']:,} |")
        w(f"| Match fraction | {b2['frac']:.3e} |")
        w(f"| 95% upper bound | {b2['ub']:.3e} ({b2['ub_method']}) |")
        w("")
    else:
        w("## Part B2 — not run")
        w("")
        w("B1 returned a clean result below the 1% threshold (not the "
          "high/inconclusive branch), so the tighter 10 ppm audit was not "
          "required per the brief.")
        w("")
    w("---")
    w("")
    w("🐕☕⬡")
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
