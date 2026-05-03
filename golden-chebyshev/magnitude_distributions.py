# magnitude_distributions.py - Magnitude distribution comparison
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a follow-up Task 3
#
# Per Mr Code Brief Addendum 3 (3 May 2026, CinC), Task 3.
#
# Purpose
# -------
# Mr Adversary's question: is anything *magnitude*-related zero-specific, or
# are the magnitudes |P(chi_2)|/|P(chi_3)| at zeros just samples from the
# same distribution as at random critical-line points? Tasks 1 and 2 confirm
# the *phase* lock is a critical-line property; Task 3 asks whether the
# magnitudes carry any zero-specific signal.
#
# Method
# ------
# 1. Compute |P(chi_2)|/|P(chi_3)| at first 200 zeta zeros (extending the
#    PR #2 sweep from 50 to 200). Output: phase_lock_results_n200.csv
#    (separate from PR #2's phase_lock_results.csv, which is preserved).
# 2. Compute the same at 200 random critical-line points sampled from
#    [14, 400] (extended t-range to cover gamma_1..gamma_200 ~ [14, 397]),
#    seed=42, same rejection rule with REJECT_WINDOW=250 to cover the
#    extended range. Output: phase_lock_null_results_n200.csv.
# 3. Compare the two magnitude distributions with a two-sample Kolmogorov-
#    Smirnov test. Plot overlaid log-x histograms.
# 4. Watch for the rho_3 pattern: at N=50, rho_3 had the largest |P(chi_2)|/
#    |P(chi_3)| (9.22) AND the smallest off-line residual (0.005). Does
#    extreme-magnitude correlate with anything in the larger sample? Report
#    top-5 extremes (largest and smallest) with gamma values for inspection.
#
# Outcome interpretation (per brief)
# ----------------------------------
# - Distributions coincide (K-S p > 0.1): zeros are not magnitude-special;
#   Paper 125's "handoff" language is decorative; Paper 196 should adjust.
# - Distributions differ (K-S p < 0.05): zeros ARE magnitude-special; there
#   is a zero-specific phenomenon worth a paper.
#
# Exploratory. No pass/fail. Report findings honestly either way.
#
# Performance note: 200 zetazero calls + 400 P_2_and_3 calls at dps=50 takes
# several minutes. Optimization: P_2_and_3 caches the four S(s, r) values
# and combines them for both characters in one pass, halving the Hurwitz
# zeta work per s.

import csv
import os
import random
import sys
import time

import numpy as np

import mpmath as mp
from mpmath import mp as _mp_ctx, mpc, mpf, fabs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hurwitz_120_mod5 import S
from character_projections import CHI_MOD5
from phase_lock_test import locked_targets, classify_phase

# Reproducibility constants for the N=200 extension. Distinct from the N=50
# constants in phase_lock_null.py (different t-range, different reject
# window) because the t-range was extended to match gamma_1..gamma_200.
SEED = 42
N = 200
T_MIN = 14.0
T_MAX = 400.0           # gamma_200 ~ 396.4
REJECT_RADIUS = 0.1
REJECT_WINDOW = 250     # cover all zeros with gamma <= T_MAX comfortably


def P_2_and_3(s):
    """Return (P(chi_2)(s), P(chi_3)(s)) sharing the four S(s, r) calls.

    Halves the Hurwitz zeta work vs. calling P(s, 2) and P(s, 3) separately.
    """
    s_vals = [None, None, None, None, None]  # 1-indexed
    for r in (1, 2, 3, 4):
        s_vals[r] = S(s, r)
    chi2 = CHI_MOD5[2]
    chi3 = CHI_MOD5[3]
    p2 = mpc(0)
    p3 = mpc(0)
    for r in (1, 2, 3, 4):
        p2 += complex(chi2[r]).conjugate() * s_vals[r]
        p3 += complex(chi3[r]).conjugate() * s_vals[r]
    return p2, p3


def compute_at(s):
    """Compute (phase, magnitude, direction, residual) at complex s."""
    p2, p3 = P_2_and_3(s)
    ratio = p2 / p3
    phase = mp.arg(ratio)
    target_minus, target_plus = locked_targets()
    direction, residual = classify_phase(phase, target_minus, target_plus)
    return phase, fabs(ratio), direction, residual


def precompute_zeros(n):
    """Return list of imaginary parts of first n zeta zeros at current dps."""
    out = []
    t0 = time.time()
    for k in range(1, n + 1):
        out.append(_mp_ctx.zetazero(k).imag)
        if k % 50 == 0:
            print("    zeros 1..%d computed (%.1fs elapsed)" % (k, time.time() - t0))
    return out


def too_close_to_zero(t, gammas, radius):
    t_mp = mpf(t)
    for g in gammas:
        if fabs(t_mp - g) < radius:
            return True
    return False


def sample_t_avoiding_zeros(n, t_min, t_max, gammas, radius, rng):
    accepted = []
    n_rejected = 0
    while len(accepted) < n:
        t = rng.uniform(t_min, t_max)
        if too_close_to_zero(t, gammas, radius):
            n_rejected += 1
        else:
            accepted.append(t)
    return accepted, n_rejected


def ks_2samp_basic(sample1, sample2):
    """Two-sample Kolmogorov-Smirnov test with asymptotic p-value.

    Standard textbook implementation - K-S statistic D is the supremum of
    |F1(x) - F2(x)| over the combined sample; p-value via the asymptotic
    Kolmogorov distribution (good for n >= 50 per side).
    """
    s1 = np.sort(np.asarray(sample1, dtype=float))
    s2 = np.sort(np.asarray(sample2, dtype=float))
    n1 = len(s1)
    n2 = len(s2)
    data_all = np.concatenate([s1, s2])
    cdf1 = np.searchsorted(s1, data_all, side='right') / n1
    cdf2 = np.searchsorted(s2, data_all, side='right') / n2
    D = np.max(np.abs(cdf1 - cdf2))
    n_eff = n1 * n2 / (n1 + n2)
    # Stephens (1970) modified statistic for asymptotic distribution.
    lam = (np.sqrt(n_eff) + 0.12 + 0.11 / np.sqrt(n_eff)) * D
    # Q(lambda) = 2 * sum_{k=1}^inf (-1)^(k-1) exp(-2 k^2 lambda^2)
    p = 0.0
    for k in range(1, 100):
        term = 2.0 * ((-1) ** (k - 1)) * np.exp(-2.0 * k * k * lam * lam)
        p += term
        if abs(term) < 1e-12:
            break
    p = max(0.0, min(1.0, p))
    return D, p


def main():
    _mp_ctx.dps = 50

    print("magnitude_distributions - N=%d zeros vs N=%d random points" % (N, N))
    print("=" * 76)
    print("dps=%d, T-range=[%.0f, %.0f], seed=%d, reject |t-gamma_k|<%s for k<=%d"
          % (_mp_ctx.dps, T_MIN, T_MAX, SEED, REJECT_RADIUS, REJECT_WINDOW))
    print()

    # ---- Step A: pre-compute first 250 zeros (covers both sweep + reject) -
    print("Step A: pre-computing first %d zeta zeros (used by both sweeps)..."
          % REJECT_WINDOW)
    gammas_all = precompute_zeros(REJECT_WINDOW)
    gammas_for_sweep = gammas_all[:N]   # first 200 for the zero sweep
    print("  gamma_1   = %s" % mp.nstr(gammas_all[0], 10))
    print("  gamma_200 = %s" % mp.nstr(gammas_all[N - 1], 10))
    print("  gamma_250 = %s" % mp.nstr(gammas_all[-1], 10))
    print()

    # ---- Step B: zero-based sweep at N=200 --------------------------------
    print("Step B: computing |P(chi_2)|/|P(chi_3)| at first %d zeta zeros..." % N)
    zero_results = []
    t0 = time.time()
    for k, gamma in enumerate(gammas_for_sweep, start=1):
        s = mpc(mpf("0.5"), gamma)
        phase, mag, direction, residual = compute_at(s)
        zero_results.append({
            'k': k, 'gamma': gamma, 'phase': phase, 'magnitude': mag,
            'direction': direction, 'residual': residual,
        })
        if k % 50 == 0:
            print("    zero %d / %d done (%.1fs elapsed)" % (k, N, time.time() - t0))
    print("    done.")
    print()

    # ---- Step C: random null sweep at N=200 -------------------------------
    print("Step C: sampling and computing at %d random critical-line points..." % N)
    rng = random.Random(SEED)
    t_values, n_rejected = sample_t_avoiding_zeros(
        N, T_MIN, T_MAX, gammas_all, mpf(str(REJECT_RADIUS)), rng)
    accept_rate = 100.0 * N / (N + n_rejected)
    print("  sampled %d t-values; rejected %d (acceptance %.1f%%)"
          % (N, n_rejected, accept_rate))

    null_results = []
    t0 = time.time()
    for i, t in enumerate(t_values, start=1):
        s = mpc(mpf("0.5"), mpf(t))
        phase, mag, direction, residual = compute_at(s)
        null_results.append({
            'index': i, 't': t, 'phase': phase, 'magnitude': mag,
            'direction': direction, 'residual': residual,
        })
        if i % 50 == 0:
            print("    point %d / %d done (%.1fs elapsed)" % (i, N, time.time() - t0))
    print("    done.")
    print()

    # ---- Step D: write the two extended CSVs ------------------------------
    zero_csv = 'phase_lock_results_n200.csv'
    with open(zero_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['k', 't_gamma', 'phase_radians', 'magnitude_ratio',
                    'direction', 'residual_dps50',
                    'dps80_rerun_triggered', 'residual_dps80'])
        for r in zero_results:
            w.writerow([
                r['k'], mp.nstr(r['gamma'], 40),
                mp.nstr(r['phase'], 40), mp.nstr(r['magnitude'], 20),
                r['direction'], mp.nstr(r['residual'], 6),
                'false', '',
            ])
    print("CSV written: %s  (%d rows)" % (zero_csv, N))

    null_csv = 'phase_lock_null_results_n200.csv'
    with open(null_csv, 'w', newline='') as f:
        f.write("# Critical-line null test, N=200 extension (Task 3)\n")
        f.write("# Seed: %d  N: %d  t-range: [%.1f, %.1f]  dps: %d\n"
                % (SEED, N, T_MIN, T_MAX, _mp_ctx.dps))
        f.write("# Rejection: |t-gamma_k| >= %s for k <= %d  "
                "(rejected %d, acceptance %.1f%%)\n"
                % (REJECT_RADIUS, REJECT_WINDOW, n_rejected, accept_rate))
        w = csv.writer(f)
        w.writerow(['index', 't', 'phase_radians', 'magnitude_ratio',
                    'direction', 'residual'])
        for r in null_results:
            w.writerow([
                r['index'], "%.17g" % r['t'],
                mp.nstr(r['phase'], 40), mp.nstr(r['magnitude'], 20),
                r['direction'], mp.nstr(r['residual'], 6),
            ])
    print("CSV written: %s  (%d rows + 3 header lines)" % (null_csv, N))
    print()

    # ---- Step E: extract magnitudes and run K-S ---------------------------
    zero_mags = np.array([float(r['magnitude']) for r in zero_results])
    null_mags = np.array([float(r['magnitude']) for r in null_results])

    log_zero = np.log10(zero_mags)
    log_null = np.log10(null_mags)

    D, p = ks_2samp_basic(log_zero, log_null)

    print("Magnitude distribution summary:")
    print("                     |   zero-based   |   null (random)")
    print("  -------------------+----------------+----------------")
    print("  N                  | %14d | %14d" % (N, N))
    print("  median magnitude   | %14.4f | %14.4f"
          % (np.median(zero_mags), np.median(null_mags)))
    print("  mean magnitude     | %14.4f | %14.4f"
          % (np.mean(zero_mags), np.mean(null_mags)))
    print("  min  magnitude     | %14.4g | %14.4g"
          % (zero_mags.min(), null_mags.min()))
    print("  max  magnitude     | %14.4g | %14.4g"
          % (zero_mags.max(), null_mags.max()))
    print("  median log10(mag)  | %14.4f | %14.4f"
          % (np.median(log_zero), np.median(log_null)))
    print("  std    log10(mag)  | %14.4f | %14.4f"
          % (np.std(log_zero), np.std(log_null)))
    print()
    print("Two-sample Kolmogorov-Smirnov on log10(magnitude):")
    print("  D = %.4f" % D)
    print("  p = %.4g" % p)
    if p > 0.1:
        print("  -> distributions are statistically indistinguishable (p > 0.1)")
        print("     Zeros are NOT magnitude-special. Paper 196 framing should")
        print("     drop zero-specific magnitude language.")
    elif p < 0.05:
        print("  -> distributions DIFFER significantly (p < 0.05)")
        print("     Zeros ARE magnitude-special. Worth flagging for CinC.")
    else:
        print("  -> intermediate p (0.05 <= p <= 0.1); inconclusive")
    print()

    # ---- Step F: rho_3 pattern check -- top extremes ----------------------
    print("Top extreme magnitudes in N=%d zero sweep (rho_3 pattern check):" % N)
    sorted_by_mag = sorted(zero_results, key=lambda r: float(r['magnitude']))
    print("  Smallest 5:")
    for r in sorted_by_mag[:5]:
        print("    rho_%-3d  gamma=%9.4f  |P(chi_2)|/|P(chi_3)| = %g"
              % (r['k'], float(r['gamma']), float(r['magnitude'])))
    print("  Largest 5:")
    for r in sorted_by_mag[-5:]:
        print("    rho_%-3d  gamma=%9.4f  |P(chi_2)|/|P(chi_3)| = %g"
              % (r['k'], float(r['gamma']), float(r['magnitude'])))
    print()

    # ---- Step G: combined CSV + plots -------------------------------------
    combined_csv = 'magnitude_distributions.csv'
    with open(combined_csv, 'w', newline='') as f:
        f.write("# Combined magnitude data, zeros vs random; N_each=%d\n" % N)
        f.write("# K-S statistic D=%.6f, p=%.6g (on log10(magnitude))\n" % (D, p))
        w = csv.writer(f)
        w.writerow(['source', 'index_or_k', 't_or_gamma', 'magnitude'])
        for r in zero_results:
            w.writerow(['zero', r['k'], mp.nstr(r['gamma'], 12),
                        mp.nstr(r['magnitude'], 12)])
        for r in null_results:
            w.writerow(['null', r['index'], "%.6f" % r['t'],
                        mp.nstr(r['magnitude'], 12)])
    print("Combined CSV written: %s  (%d rows + 2 header lines)"
          % (combined_csv, 2 * N))

    plot_path = 'magnitude_distributions.png'
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: overlaid log-x histograms
    ax = axes[0]
    log_min = min(log_zero.min(), log_null.min()) - 0.2
    log_max = max(log_zero.max(), log_null.max()) + 0.2
    bins = np.linspace(log_min, log_max, 30)
    ax.hist(log_zero, bins=bins, alpha=0.55, label='zeros (N=%d)' % N,
            color='steelblue', edgecolor='black', linewidth=0.5)
    ax.hist(log_null, bins=bins, alpha=0.55, label='random (N=%d)' % N,
            color='darkorange', edgecolor='black', linewidth=0.5)
    ax.set_xlabel('log10 |P(chi_2)| / |P(chi_3)|')
    ax.set_ylabel('count')
    ax.set_title('Magnitude distribution: zeros vs random critical-line points\n'
                 'K-S D = %.4f, p = %.3g' % (D, p))
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Right: empirical CDFs (where K-S literally measures the gap)
    ax = axes[1]
    sorted_z = np.sort(log_zero)
    sorted_n = np.sort(log_null)
    cdf_z = np.arange(1, len(sorted_z) + 1) / len(sorted_z)
    cdf_n = np.arange(1, len(sorted_n) + 1) / len(sorted_n)
    ax.plot(sorted_z, cdf_z, drawstyle='steps-post', label='zeros',
            color='steelblue', linewidth=1.3)
    ax.plot(sorted_n, cdf_n, drawstyle='steps-post', label='random',
            color='darkorange', linewidth=1.3)
    ax.set_xlabel('log10 |P(chi_2)| / |P(chi_3)|')
    ax.set_ylabel('empirical CDF')
    ax.set_title('Empirical CDFs (K-S D = supremum gap)')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=120)
    plt.close()
    print("Plot written: %s" % plot_path)
    print()

    print("=" * 76)
    print("Task 3 (magnitude distribution): exploratory, no pass/fail.")
    print("Verdict: K-S p = %.4g  ->  %s" %
          (p,
           "indistinguishable" if p > 0.1 else
           "different" if p < 0.05 else
           "inconclusive"))


if __name__ == "__main__":
    main()
