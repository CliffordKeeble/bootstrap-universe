# phase_lock_null.py - Critical-line null test for the golden phase lock
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a follow-up Task 1
#
# Per Mr Code Brief Addendum 3 (3 May 2026, CinC), Task 1.
#
# Purpose
# -------
# Confirm that the golden phase lock
#   arg(P(chi_2)/P(chi_3)) in {-arctan(1/phi), pi - arctan(1/phi)}
# is a property of the critical line Re(s) = 1/2, NOT a property specific
# to zeta zeros.
#
# Mr Adversary's algebraic completion of Paper 125 Theorem 1 shows
#   P(chi_k)(s) = 120^s * L(s, conj(chi_k))
# so
#   P(chi_2)(s) / P(chi_3)(s) = L(s, chi_3) / L(s, chi_2),
# whose phase is locked by the L-function functional equation on Re(s) = 1/2,
# independent of whether s sits at a zeta zero or anywhere else on the line.
# Tests at random critical-line points should therefore reproduce the same
# locked phases (at the same precision) as tests at zeta zeros. This script
# provides the numerical confirmation.
#
# Method (reproducibility-critical; carry into Paper 196 methods section)
# -----------------------------------------------------------------------
# 1. Pre-compute the first 100 zeta-zero imaginary parts (gamma_k) at dps=50.
# 2. Sample t ~ Uniform(T_MIN, T_MAX) using random.seed(SEED). Reject any
#    candidate t with |t - gamma_k| < REJECT_RADIUS for k <= REJECT_WINDOW;
#    resample until N_RANDOM accepted. Rejection count is logged to console
#    and CSV header.
# 3. At each accepted t, compute s = 0.5 + i*t and the phase residual via
#    character_projections.P, sharing locked_targets()/classify_phase()
#    with phase_lock_test.py.
#
# Pass condition: all N_RANDOM points show residual < 1e-35.
#
# Expected outcome: residuals indistinguishable from the zero-based test
# (dps=50 roundoff floor, ~10^-48). If so, this is confirmation, not a
# finding. If the null distribution differs from the zero-based distribution
# in any informative way, that IS a finding and we stop and report.

import csv
import os
import random
import sys

import mpmath as mp
from mpmath import mp as _mp_ctx, mpc, mpf, fabs

from character_projections import P
from phase_lock_test import locked_targets, classify_phase

# Reproducibility constants. DO NOT alter without recording the change in
# the Paper 196 methods section: re-running this script with these constants
# will reproduce the exact 50 random t-values used here.
SEED = 42
N_RANDOM = 50
T_MIN = 14.0           # ~ gamma_1 = 14.13
T_MAX = 144.0          # ~ gamma_50 = 143.7
REJECT_RADIUS = 0.1    # rejection-zone half-width around each gamma_k
REJECT_WINDOW = 100    # check candidates against gamma_k for k = 1..100


def precompute_zeros(n):
    """Return list of imaginary parts of first n zeta zeros at current dps."""
    return [_mp_ctx.zetazero(k).imag for k in range(1, n + 1)]


def too_close_to_zero(t, gammas, radius):
    """True if t is within `radius` of any gamma in gammas (mpf comparison)."""
    t_mp = mpf(t)
    for g in gammas:
        if fabs(t_mp - g) < radius:
            return True
    return False


def sample_t_avoiding_zeros(n, t_min, t_max, gammas, radius, rng):
    """Draw n t-values uniform on (t_min, t_max), rejecting those near zeros.

    Returns (accepted_list, rejected_count).
    """
    accepted = []
    n_rejected = 0
    while len(accepted) < n:
        t = rng.uniform(t_min, t_max)
        if too_close_to_zero(t, gammas, radius):
            n_rejected += 1
        else:
            accepted.append(t)
    return accepted, n_rejected


def load_zero_residuals(csv_path):
    """Load effective residuals from PR #2's phase_lock_results.csv."""
    if not os.path.exists(csv_path):
        return None
    residuals = []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            r50 = mpf(row['residual_dps50'])
            r80_raw = row.get('residual_dps80', '').strip()
            r80 = mpf(r80_raw) if r80_raw else None
            rerun = row['dps80_rerun_triggered'] == 'true'
            residuals.append(r80 if rerun and r80 is not None else r50)
    return residuals


def main():
    _mp_ctx.dps = 50
    tol = mpf(10) ** -35

    print("phase_lock_null - critical-line null test for the golden phase lock")
    print("=" * 72)
    print("dps = %d, target tol = 1e-35, N_RANDOM = %d"
          % (_mp_ctx.dps, N_RANDOM))
    print("t-range: [%.1f, %.1f]  (covers first 50 zeta zeros)"
          % (T_MIN, T_MAX))
    print("Seed = %d, reject |t - gamma_k| < %s for k <= %d"
          % (SEED, REJECT_RADIUS, REJECT_WINDOW))
    print()

    # Step 1: pre-compute zeros for the rejection check.
    print("Pre-computing first %d zeta zeros for rejection check..."
          % REJECT_WINDOW)
    gammas = precompute_zeros(REJECT_WINDOW)
    print("  gamma_1   = %s" % mp.nstr(gammas[0], 12))
    print("  gamma_50  = %s" % mp.nstr(gammas[49], 12))
    print("  gamma_100 = %s" % mp.nstr(gammas[-1], 12))
    print()

    # Step 2: sample t-values, rejecting those near zeros.
    rng = random.Random(SEED)
    t_values, n_rejected = sample_t_avoiding_zeros(
        N_RANDOM, T_MIN, T_MAX, gammas, mpf(str(REJECT_RADIUS)), rng)
    accept_rate = 100.0 * N_RANDOM / (N_RANDOM + n_rejected)
    print("Sampled %d t-values  |  rejected %d candidates  |  acceptance %.1f%%"
          % (N_RANDOM, n_rejected, accept_rate))
    print()

    # Step 3: compute phase lock at each random t.
    target_minus, target_plus = locked_targets()
    print("Per-point progress ('.' = below tol, 'F' = at/above tol):")
    results = []
    for i, t in enumerate(t_values, start=1):
        s = mpc(mpf("0.5"), mpf(t))
        p2 = P(s, 2)
        p3 = P(s, 3)
        ratio = p2 / p3
        phase = mp.arg(ratio)
        magnitude = fabs(ratio)
        direction, residual = classify_phase(phase, target_minus, target_plus)
        results.append({
            'index': i, 't': t, 'phase': phase, 'magnitude': magnitude,
            'direction': direction, 'residual': residual,
        })
        sys.stdout.write('.' if residual < tol else 'F')
        sys.stdout.flush()
    print()
    print()

    # ---- Summary ------------------------------------------------------------
    pass_count = sum(1 for r in results if r['residual'] < tol)
    max_res = max(r['residual'] for r in results)
    mean_res = sum(r['residual'] for r in results) / mpf(len(results))
    mag_min = min(r['magnitude'] for r in results)
    mag_max = max(r['magnitude'] for r in results)
    n_plus = sum(1 for r in results if r['direction'] == '+')
    n_minus = N_RANDOM - n_plus

    print("Null sweep summary (%d random critical-line points):" % N_RANDOM)
    print("  Pass count:       %d / %d" % (pass_count, N_RANDOM))
    print("  Max residual:     %s" % mp.nstr(max_res, 6))
    print("  Mean residual:    %s" % mp.nstr(mean_res, 6))
    print("  Magnitude range:  [%s, %s]"
          % (mp.nstr(mag_min, 4), mp.nstr(mag_max, 4)))
    print("  Direction split:  + : %d, - : %d" % (n_plus, n_minus))
    print()

    # ---- Side-by-side comparison vs zero-based residuals --------------------
    zero_csv = 'phase_lock_results.csv'
    zero_residuals = load_zero_residuals(zero_csv)
    if zero_residuals:
        zero_max = max(zero_residuals)
        zero_mean = sum(zero_residuals) / mpf(len(zero_residuals))
        print("Side-by-side comparison vs zero-based test (%d zeros from %s):"
              % (len(zero_residuals), zero_csv))
        print("                 |   zero-based   |   null (random)")
        print("  ---------------+----------------+----------------")
        print("  Max residual   | %14s | %14s"
              % (mp.nstr(zero_max, 4), mp.nstr(max_res, 4)))
        print("  Mean residual  | %14s | %14s"
              % (mp.nstr(zero_mean, 4), mp.nstr(mean_res, 4)))
        if zero_max > 0:
            print("  Ratio max(null)/max(zero) = %s"
                  % mp.nstr(max_res / zero_max, 4))
        floor_threshold = mpf(10) ** -40
        if max_res < floor_threshold and zero_max < floor_threshold:
            print()
            print("  Both at dps=50 roundoff floor (< 1e-40).")
            print("  Lock confirmed as a critical-line property; zeros are")
            print("  not the locus of the phenomenon.")
        else:
            print()
            print("  *** Distributions NOT both at floor - investigate. ***")
    else:
        print("(zero-based comparison skipped - %s not found)" % zero_csv)
    print()

    # ---- CSV output ---------------------------------------------------------
    csv_path = 'phase_lock_null_results.csv'
    with open(csv_path, 'w', newline='') as f:
        f.write("# Critical-line null test for the golden phase lock\n")
        f.write("# Seed: %d  N: %d  t-range: [%.1f, %.1f]  dps: %d\n"
                % (SEED, N_RANDOM, T_MIN, T_MAX, _mp_ctx.dps))
        f.write("# Rejection: |t - gamma_k| >= %s for k <= %d  "
                "(rejected %d candidates, acceptance %.1f%%)\n"
                % (REJECT_RADIUS, REJECT_WINDOW, n_rejected, accept_rate))
        w = csv.writer(f)
        w.writerow(['index', 't', 'phase_radians', 'magnitude_ratio',
                    'direction', 'residual'])
        for r in results:
            w.writerow([
                r['index'],
                "%.17g" % r['t'],
                mp.nstr(r['phase'], 40),
                mp.nstr(r['magnitude'], 20),
                r['direction'],
                mp.nstr(r['residual'], 6),
            ])
    print("CSV written: %s  (%d rows + 3 header lines)"
          % (csv_path, N_RANDOM))
    print()

    print("=" * 72)
    overall_pass = (pass_count == N_RANDOM)
    print("Task 1 (null test): %s"
          % ("PASS - lock confirmed as critical-line property"
             if overall_pass else "FAIL - investigate before continuing"))


if __name__ == "__main__":
    main()
