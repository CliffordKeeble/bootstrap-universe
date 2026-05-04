# zeros_data.py - shared zero-data loader for the Riemann x digit-gap test
# Bootstrap Universe Programme - Brief 4 May 2026 (CinC), riemann_digit_gap_test
#
# Provides cached lists of the first ~100 positive imaginary parts of zeros
# of zeta(s) and L(s, chi_5) for use by the three tests in this folder.
# Computes them once on first run and caches as CSV; subsequent runs load
# instantly. CSVs are committed to the repo so the tests are reproducible
# without waiting for a fresh zero-finder.
#
# Both data sources are CANONICAL and well-tested:
#   - zeta zeros   : mpmath.zetazero(k) at dps=30; first 100 zeros take ~30 s
#                    to compute. Cross-checked against Odlyzko's tables to
#                    ~12 decimal digits in the LMFDB-anchored Layer 3a work.
#   - L(chi_5) zero: extends golden-chebyshev/lchi5_zeros.py: first 25 are
#                    LMFDB-anchored (4 decimals) + Lambda-bisection refined
#                    to dps=30; zeros 26+ found by Lambda-sign-change sweep
#                    on the critical line, refined by the same bisection.
#
# Selberg-null: explicit functions Delta_selberg(T), N_zeta_selberg(T),
# N_chi5_selberg(T). The constants are documented inline.

import csv
import os
import sys
from pathlib import Path

import mpmath
from mpmath import mp, mpf, mpc, fabs

# Allow imports from the sibling golden-chebyshev/ directory for the
# already-validated Lambda_chi5 + LMFDB seed list.
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "golden-chebyshev"))
from lchi5_zeros import (
    LMFDB_ZEROS_25,
    Lambda_chi5,
    refine_zero as refine_lchi5_zero,
    CHI5_VALUES,
)

CACHE_DIR = Path(__file__).resolve().parent / "zeros_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
ZETA_CSV   = CACHE_DIR / "riemann_zeros_dps30_n100.csv"
LCHI5_CSV  = CACHE_DIR / "lchi5_zeros_dps30_n150.csv"

DEFAULT_DPS = 30


# ----------------------------------------------------------------------
# Selberg analytic null - written out explicitly per Pattern 14 / 28.
# ----------------------------------------------------------------------
#
# Riemann-von Mangoldt for zeta(s):
#     N_zeta(T) = (T/2pi) log(T/2pi e) + 7/8 + S(T) + O(1/T)
# where S(T) = (1/pi) arg zeta(1/2 + iT) is the bounded oscillatory term.
# The "smooth" Selberg main term (= the analytic null we subtract):
#     N_zeta_smooth(T) = (T/2pi) log(T/2pi e) + 7/8
#
# For a primitive Dirichlet character chi mod q, the corresponding formula
# is (e.g. Iwaniec-Kowalski Thm 5.8):
#     N_chi(T) = (T/2pi) log(qT/2pi e) + delta_a + S_chi(T) + O(1/T)
# where delta_a depends on parity:
#   - even chi (a=0, gamma factor Gamma_R(s) = pi^(-s/2) Gamma(s/2)):
#       delta_a = -1/8
#   - odd chi  (a=1, gamma factor Gamma_R(s+1) = pi^(-(s+1)/2) Gamma((s+1)/2)):
#       delta_a = +3/8
#
# chi_5 (Legendre symbol n -> (n/5)): chi_5(-1) = (-1/5) = (-1)^((5-1)/2)
# = (-1)^2 = +1, so chi_5 is EVEN -> delta_a = -1/8.
# Conductor q = 5.
#
# Difference of smooth parts:
#     Delta_smooth(T) = N_zeta_smooth(T) - N_chi5_smooth(T)
#                     = (T/2pi) [log(T/2pi e) - log(5T/2pi e)] + (7/8 - (-1/8))
#                     = (T/2pi) log(1/5) + 1
#                     = -(T/2pi) log(5) + 1
#
# Slope: -log(5)/(2pi) ~ -0.2562; constant: +1.

def N_zeta_smooth(T):
    """Smooth Selberg main term for the zero-counting function of zeta(s)."""
    T = mpf(T)
    return T / (2 * mpmath.pi) * mpmath.log(T / (2 * mpmath.pi * mpmath.e)) + mpf(7) / 8


def N_chi5_smooth(T):
    """Smooth Selberg main term for L(s, chi_5).

    chi_5 is primitive of conductor 5 and EVEN (parity delta_a = -1/8).
    """
    T = mpf(T)
    q = mpf(5)
    return T / (2 * mpmath.pi) * mpmath.log(q * T / (2 * mpmath.pi * mpmath.e)) - mpf(1) / 8


def Delta_smooth(T):
    """N_zeta_smooth(T) - N_chi5_smooth(T). The Selberg null for Test 1."""
    T = mpf(T)
    return - T / (2 * mpmath.pi) * mpmath.log(5) + mpf(1)


# ----------------------------------------------------------------------
# Zeta zeros - cached from mpmath.zetazero at dps=30.
# ----------------------------------------------------------------------

def _compute_zeta_zeros(n, dps):
    print("  [zeros_data] computing first %d zeta zeros at dps=%d..."
          % (n, dps), flush=True)
    with mp.workdps(dps + 5):
        out = []
        for k in range(1, n + 1):
            rho = mpmath.zetazero(k)
            out.append(rho.imag)
            if k % 25 == 0:
                print("    ... %d / %d" % (k, n), flush=True)
    return out


def load_zeta_zeros(n=100, dps=DEFAULT_DPS, force_recompute=False):
    """Return list of first n imaginary parts of nontrivial zeta zeros."""
    if not force_recompute and ZETA_CSV.exists():
        zeros = []
        with open(ZETA_CSV, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                zeros.append(mpf(row["gamma"]))
        if len(zeros) >= n:
            return zeros[:n]
    zeros = _compute_zeta_zeros(n, dps)
    with open(ZETA_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "gamma"])
        for k, g in enumerate(zeros, start=1):
            w.writerow([k, mpmath.nstr(g, dps + 5, strip_zeros=False)])
    return zeros


# ----------------------------------------------------------------------
# L(chi_5) zeros - first 25 LMFDB-seeded, 26+ found by Lambda-sign-change.
# ----------------------------------------------------------------------

def _sweep_lambda_chi5_for_zeros(t_min, t_max, step, dps):
    """Find candidate zero brackets by Lambda(t) sign changes on a grid.

    Returns list of (a, b) bracket pairs, where Lambda(a) * Lambda(b) < 0.
    Step should be substantially below the local zero spacing
    2*pi / log(5*t / 2*pi) - at t ~ 200, spacing ~ 1.26, so step 0.4 is safe.
    """
    print("  [zeros_data] sweeping Lambda_chi5 on [%g, %g] step %g (dps=%d)..."
          % (t_min, t_max, step, dps), flush=True)
    brackets = []
    t = mpf(t_min)
    f_prev = Lambda_chi5(t, dps)
    while t < t_max:
        t_next = t + step
        if t_next > t_max:
            t_next = mpf(t_max)
        f_next = Lambda_chi5(t_next, dps)
        if f_prev * f_next < 0:
            brackets.append((t, t_next))
        t = t_next
        f_prev = f_next
    print("    found %d sign-change brackets" % len(brackets), flush=True)
    return brackets


def _compute_lchi5_zeros(n, dps):
    """Compute first n positive zeros of L(s, chi_5) at dps precision."""
    out = []
    # Stage 1: first 25 from LMFDB seeds, refined by bisection.
    print("  [zeros_data] refining first 25 LMFDB-anchored zeros...", flush=True)
    for seed in LMFDB_ZEROS_25:
        z = refine_lchi5_zero(seed, dps=dps)
        out.append(z)
    if n <= 25:
        return out[:n]
    # Stage 2: extend by Lambda sign-change sweep beyond LMFDB_ZEROS_25[-1].
    t_start = float(LMFDB_ZEROS_25[-1]) + 0.5
    # Generous upper bound: zeros at average spacing ~1.5, so first ~80 zeros
    # are within roughly t = 200. Sweep slightly past 200 for a safety margin.
    t_end = 220.0
    step = 0.4
    brackets = _sweep_lambda_chi5_for_zeros(t_start, t_end, step, dps)
    print("  [zeros_data] refining %d swept brackets..." % len(brackets),
          flush=True)
    for i, (a, b) in enumerate(brackets, start=26):
        # Use the bracket midpoint as seed; bisection-with-secant.
        mid = (a + b) / 2
        z = refine_lchi5_zero(mid, dps=dps, delta=(b - a) / 2 + mpf("0.001"))
        out.append(z)
        if (i - 25) % 25 == 0:
            print("    ... refined zero %d (gamma~%.2f)" % (i, float(z)),
                  flush=True)
        if len(out) >= n:
            break
    return out[:n]


def load_lchi5_zeros(n=150, dps=DEFAULT_DPS, force_recompute=False):
    """Return list of first n positive imaginary parts of L(s, chi_5) zeros."""
    if not force_recompute and LCHI5_CSV.exists():
        zeros = []
        with open(LCHI5_CSV, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                zeros.append(mpf(row["gamma"]))
        if len(zeros) >= n:
            return zeros[:n]
    zeros = _compute_lchi5_zeros(n, dps)
    with open(LCHI5_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "gamma"])
        for k, g in enumerate(zeros, start=1):
            w.writerow([k, mpmath.nstr(g, dps + 5, strip_zeros=False)])
    return zeros


# ----------------------------------------------------------------------
# Sanity self-test (run directly).
# ----------------------------------------------------------------------

def _self_test():
    mp.dps = DEFAULT_DPS + 5
    print("Loading zeta zeros (n=100)...")
    zz = load_zeta_zeros(n=100)
    print("  first 3:", [mpmath.nstr(x, 12) for x in zz[:3]])
    print("  count <= 140:", sum(1 for g in zz if g <= 140),
          "  (brief expects ~47)")
    print("  count <= 184:", sum(1 for g in zz if g <= 184),
          "  (brief expects ~71)")

    print()
    print("Loading L(chi_5) zeros (n=150)...")
    lz = load_lchi5_zeros(n=150)
    print("  first 3:", [mpmath.nstr(x, 12) for x in lz[:3]])
    print("  count <= 140:", sum(1 for g in lz if g <= 140))
    print("  count <= 184:", sum(1 for g in lz if g <= 184))
    print("  max gamma:", mpmath.nstr(lz[-1], 8))

    print()
    print("Selberg null sanity check at T=140:")
    print("  N_zeta_smooth(140)  =", mpmath.nstr(N_zeta_smooth(140), 6))
    print("  N_chi5_smooth(140)  =", mpmath.nstr(N_chi5_smooth(140), 6))
    print("  Delta_smooth(140)   =", mpmath.nstr(Delta_smooth(140), 6))
    print("  -(T/2pi) log 5 + 1  =", mpmath.nstr(
        -mpf(140) / (2 * mpmath.pi) * mpmath.log(5) + 1, 6))
    print()
    print("  Selberg null at T=184:")
    print("  Delta_smooth(184)   =", mpmath.nstr(Delta_smooth(184), 6))

    # ------------------------------------------------------------------
    # Calibration audit (Pattern 14: the calibration is itself a lemma).
    # Per CinC 4 May 2026: residuals at non-boundary heights should sit
    # near zero with typical magnitude ~ sqrt(log T) / pi ~ 0.7. If all
    # five non-boundary residuals are within +/- 1.5 with no obvious
    # systematic drift, the linear-plus-constant null is calibrated and
    # the windowed boundary test is meaningful. If any is outside, the
    # constants need revisiting before the test runs.
    # ------------------------------------------------------------------
    print()
    print("Calibration audit at non-boundary heights:")
    print("  expected typical |R| ~ sqrt(log T) / pi ~ 0.7 at this scale")
    print("  pass condition: all |R| < 1.5  (no systematic drift)")
    print()
    audit_T = [60, 80, 100, 120, 160]
    audit_rows = []
    for T in audit_T:
        n_z = sum(1 for g in zz if g <= T)
        n_l = sum(1 for g in lz if g <= T)
        delta_emp = n_z - n_l
        delta_pred = float(Delta_smooth(T))
        R = delta_emp - delta_pred
        expected = float(mpmath.sqrt(mpmath.log(T)) / mpmath.pi)
        audit_rows.append((T, n_z, n_l, delta_emp, delta_pred, R, expected))

    print("  | T   | N_zeta | N_chi5 | Delta_emp | Delta_pred | R     | expected ~ |")
    print("  |-----|--------|--------|-----------|------------|-------|------------|")
    for T, nz, nl, de, dp, R, exp in audit_rows:
        print("  | %3d | %6d | %6d | %9d | %10.4f | %+5.3f | %10.3f |"
              % (T, nz, nl, de, dp, R, exp))

    audit_pass = all(abs(R) < 1.5 for _, _, _, _, _, R, _ in audit_rows)
    print()
    R_values = [R for _, _, _, _, _, R, _ in audit_rows]
    R_mean = sum(R_values) / len(R_values)
    print("  R mean (drift indicator): %+.4f  (small = no drift)" % R_mean)
    print("  R range: [%+.3f, %+.3f]" % (min(R_values), max(R_values)))
    if audit_pass:
        print("  CALIBRATION: PASS - linear-plus-constant null is calibrated.")
        print("  -> windowed boundary test is meaningful; safe to proceed.")
    else:
        print("  CALIBRATION: FAIL - residual exceeds +/-1.5 at a non-boundary height.")
        print("  -> Selberg constants need revisiting before the boundary test runs.")


if __name__ == "__main__":
    _self_test()
