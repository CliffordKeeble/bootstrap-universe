# sign_distribution.py - Sign sequence analysis at first 200 zeta zeros
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a follow-up Task 4
#
# Per Mr Code Brief Addendum 3 (3 May 2026, CinC), Task 4.
#
# Purpose
# -------
# The +/- direction at each zero records which side of the locked line
# arg(P(chi_2)/P(chi_3)) lands on. After Tasks 1-3 we know:
#   - The phase lock is a critical-line property, not zero-specific (Task 1)
#   - It breaks symmetrically off-line (Task 2)
#   - The magnitudes at zeros are not zero-special (Task 3, K-S p=0.53)
# So Task 4 asks the last open zero-specific question: does the +/- sequence
# carry any structure - periodicity, autocorrelation, correlation with zero
# spacing? If yes: there's a zero-specific finding worth a paper. If no:
# Paper 196 records the rule-out and moves on.
#
# Pre-test observation: at N=50, the direction split was 21/29 in PR #2's
# zero-based test and 19/31 in Task 1's null. Both ~60/40 but NOT
# significant at N=50. At N=200, the same ratio WOULD be significant
# (binomial threshold for 5%% is roughly 87/113 vs 100/100). This script
# checks whether the imbalance survives.
#
# Tests
# -----
#   1. Direction count + binomial test for fairness (H0: p = 0.5)
#   2. Wald-Wolfowitz runs test (sequence-randomness given the marginal)
#   3. Autocorrelation at lags 1..20, with white-noise 95%% band
#   4. Pearson correlation of sign and sign-change with zero spacing
#
# Inputs:  phase_lock_results_n200.csv (produced by magnitude_distributions.py)
# Outputs: sign_analysis.csv, sign_runs.txt, sign_plots.png
#
# Exploratory; no pass/fail. Report what's there and what isn't.

import csv
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ZERO_CSV = 'phase_lock_results_n200.csv'
N_LAGS = 20


def two_sided_p_normal(z):
    """Two-sided p-value from a standard-normal z statistic."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def load_signs_and_gammas(csv_path):
    """Read phase_lock_results_n200.csv -> (gammas array, signs array as +/-1)."""
    gammas = []
    signs = []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            gammas.append(float(row['t_gamma']))
            d = row['direction']
            if d == '+':
                signs.append(+1)
            elif d == '-':
                signs.append(-1)
            else:
                raise ValueError("unexpected direction value: %r" % d)
    return np.array(gammas), np.array(signs)


def binomial_normal_test(n_plus, N):
    """Two-sided binomial test (H0: p=0.5) via normal approximation."""
    z = (n_plus - N / 2.0) / math.sqrt(N / 4.0)
    return z, two_sided_p_normal(z)


def runs_test(signs):
    """Wald-Wolfowitz runs test on a binary +/-1 sequence.

    Returns (observed_runs, mu_R, sigma_R, z, p_two_sided), or None if
    one class is empty (test undefined).
    """
    N = len(signs)
    n_plus = int(np.sum(signs == +1))
    n_minus = N - n_plus
    if n_plus == 0 or n_minus == 0:
        return None
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    mu_R = 2.0 * n_plus * n_minus / N + 1.0
    var_R = (2.0 * n_plus * n_minus * (2.0 * n_plus * n_minus - N)
             / (N * N * (N - 1)))
    sigma_R = math.sqrt(var_R)
    z = (runs - mu_R) / sigma_R
    return runs, mu_R, sigma_R, z, two_sided_p_normal(z)


def autocorrelation(signs, max_lag):
    """Autocorrelation r_h for h = 1..max_lag (Bartlett-style)."""
    s = signs.astype(float)
    s_centered = s - s.mean()
    var = float(np.sum(s_centered * s_centered))
    rs = []
    for h in range(1, max_lag + 1):
        rs.append(float(np.sum(s_centered[:-h] * s_centered[h:]) / var))
    return np.array(rs)


def main():
    print("sign_distribution - sign sequence analysis at first 200 zeta zeros")
    print("=" * 76)

    if not os.path.exists(ZERO_CSV):
        print("ERROR: %s not found. Run magnitude_distributions.py first."
              % ZERO_CSV)
        return 1

    gammas, signs = load_signs_and_gammas(ZERO_CSV)
    N = len(signs)
    n_plus = int(np.sum(signs == +1))
    n_minus = N - n_plus
    spacings = np.diff(gammas)  # length N-1
    print("Loaded %d signs from %s" % (N, ZERO_CSV))
    print("  +: %d  (%.1f%%)" % (n_plus, 100.0 * n_plus / N))
    print("  -: %d  (%.1f%%)" % (n_minus, 100.0 * n_minus / N))
    print("  spacing range: [%.4f, %.4f], mean = %.4f"
          % (spacings.min(), spacings.max(), spacings.mean()))
    print()

    # -----------------------------------------------------------------
    # Test 1: direction count via two-sided binomial
    # -----------------------------------------------------------------
    print("Test 1: marginal fairness (binomial, H0: p = 0.5)")
    z_bin, p_bin = binomial_normal_test(n_plus, N)
    print("  observed: %d / %d  (%.0f / %.0f split)"
          % (n_plus, N, 100.0 * n_plus / N, 100.0 * n_minus / N))
    print("  z = %+.4f, p = %.4f"
          % (z_bin, p_bin))
    if p_bin < 0.05:
        print("  -> deviation from 50/50 IS significant at the 5%% level")
    elif p_bin < 0.1:
        print("  -> deviation from 50/50 is borderline (5-10%% range)")
    else:
        print("  -> deviation from 50/50 is NOT significant (consistent with fair coin)")
    print()

    # -----------------------------------------------------------------
    # Test 2: Wald-Wolfowitz runs test
    # -----------------------------------------------------------------
    print("Test 2: Wald-Wolfowitz runs test (sequence randomness given marginal)")
    rt = runs_test(signs)
    runs, mu_R, sigma_R, z_r, p_r = rt
    print("  observed runs:  %d" % runs)
    print("  expected runs:  %.2f +/- %.2f  (under H0: random ordering)" % (mu_R, sigma_R))
    print("  z = %+.4f, p = %.4f" % (z_r, p_r))
    if p_r < 0.05:
        if runs < mu_R:
            print("  -> FEWER runs than expected: clustering of like signs (significant)")
        else:
            print("  -> MORE runs than expected: anti-clustering / oscillation (significant)")
    else:
        print("  -> consistent with random ordering (no runs structure)")
    print()

    # -----------------------------------------------------------------
    # Test 3: autocorrelation up to lag 20
    # -----------------------------------------------------------------
    print("Test 3: autocorrelation (lags 1..%d)" % N_LAGS)
    autocorrs = autocorrelation(signs, N_LAGS)
    ci_95 = 1.96 / math.sqrt(N)
    print("  95%% white-noise band: +/- %.4f" % ci_95)
    print("  lag |   r_lag    | outside band?")
    print("  ----+------------+--------------")
    n_outside = 0
    for h, r in enumerate(autocorrs, start=1):
        flag = "YES" if abs(r) > ci_95 else "no"
        if abs(r) > ci_95:
            n_outside += 1
        print("  %3d | %+10.4f | %s" % (h, r, flag))
    expected_outside = N_LAGS * 0.05
    print()
    if n_outside == 0:
        print("  -> 0 / %d lags exceed the band; no autocorrelation structure"
              % N_LAGS)
    else:
        print("  -> %d / %d lags exceed band (expected ~%.1f by chance)"
              % (n_outside, N_LAGS, expected_outside))
        if n_outside > 2 * expected_outside:
            print("     -> notably more than chance; possible autocorrelation structure")
        else:
            print("     -> within chance; no clear autocorrelation structure")
    print()

    # -----------------------------------------------------------------
    # Test 4: correlation with zero spacing
    # -----------------------------------------------------------------
    print("Test 4: correlation with zero spacing (gamma_{k+1} - gamma_k)")
    s_truncated = signs[:-1].astype(float)               # signs aligned with spacings
    sign_changes = (signs[1:] != signs[:-1]).astype(float)  # 1 = sign flips at k->k+1
    pearson_sign_spacing = float(np.corrcoef(s_truncated, spacings)[0, 1])
    pearson_change_spacing = float(np.corrcoef(sign_changes, spacings)[0, 1])
    spacing_at_plus = spacings[s_truncated == +1]
    spacing_at_minus = spacings[s_truncated == -1]
    spacing_at_change = spacings[sign_changes == 1]
    spacing_at_nochange = spacings[sign_changes == 0]
    ci_corr = 1.96 / math.sqrt(N - 1)
    print("  Pearson(sign[k], spacing[k])        = %+.4f   (95%% CI ~ +/- %.4f)"
          % (pearson_sign_spacing, ci_corr))
    print("  Pearson(sign-change[k], spacing[k]) = %+.4f"
          % pearson_change_spacing)
    print()
    print("  Conditional means:")
    print("    spacing | sign = +     : %.4f  (n = %d)"
          % (spacing_at_plus.mean(), len(spacing_at_plus)))
    print("    spacing | sign = -     : %.4f  (n = %d)"
          % (spacing_at_minus.mean(), len(spacing_at_minus)))
    print("    spacing | sign change  : %.4f  (n = %d)"
          % (spacing_at_change.mean(), len(spacing_at_change)))
    print("    spacing | same sign    : %.4f  (n = %d)"
          % (spacing_at_nochange.mean(), len(spacing_at_nochange)))
    if abs(pearson_sign_spacing) > ci_corr or abs(pearson_change_spacing) > ci_corr:
        print("  -> at least one correlation exceeds the 95%% null band")
    else:
        print("  -> no correlation exceeds the 95%% null band")
    print()

    # -----------------------------------------------------------------
    # Outputs
    # -----------------------------------------------------------------
    csv_path = 'sign_analysis.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['k', 'gamma', 'sign', 'spacing_to_next'])
        for k in range(N):
            sp = "%.6f" % spacings[k] if k < N - 1 else ''
            w.writerow([k + 1, "%.10f" % gammas[k],
                        '+' if signs[k] == 1 else '-', sp])
    print("CSV written: %s  (%d rows)" % (csv_path, N))

    txt_path = 'sign_runs.txt'
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Sign sequence analysis - Layer 3a follow-up Task 4\n")
        f.write("=" * 60 + "\n\n")
        f.write("Source: %s  (N=%d)\n" % (ZERO_CSV, N))
        f.write("Counts: + = %d, - = %d  (%.1f%% / %.1f%%)\n\n"
                % (n_plus, n_minus, 100.0 * n_plus / N, 100.0 * n_minus / N))
        f.write("Test 1 (binomial, two-sided, H0: p=0.5):\n")
        f.write("  z = %+.4f, p = %.4f\n\n" % (z_bin, p_bin))
        f.write("Test 2 (Wald-Wolfowitz runs):\n")
        f.write("  observed = %d, expected = %.2f +/- %.2f\n"
                % (runs, mu_R, sigma_R))
        f.write("  z = %+.4f, p = %.4f\n\n" % (z_r, p_r))
        f.write("Test 3 (autocorrelation, 95%% band = +/- %.4f):\n" % ci_95)
        f.write("  lag    r_lag\n")
        for h, r in enumerate(autocorrs, start=1):
            mark = "  *" if abs(r) > ci_95 else ""
            f.write("  %3d  %+8.4f%s\n" % (h, r, mark))
        f.write("\n")
        f.write("Test 4 (spacing correlations, 95%% band = +/- %.4f):\n" % ci_corr)
        f.write("  Pearson(sign,   spacing) = %+.4f\n" % pearson_sign_spacing)
        f.write("  Pearson(change, spacing) = %+.4f\n" % pearson_change_spacing)
        f.write("  mean spacing | sign = + : %.4f\n" % spacing_at_plus.mean())
        f.write("  mean spacing | sign = - : %.4f\n" % spacing_at_minus.mean())
        f.write("  mean spacing | change   : %.4f\n" % spacing_at_change.mean())
        f.write("  mean spacing | same     : %.4f\n" % spacing_at_nochange.mean())
    print("Stats written: %s" % txt_path)

    # ---- Plots --------------------------------------------------------------
    plot_path = 'sign_plots.png'
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # 1. Sign vs gamma (strip)
    ax = axes[0, 0]
    plus_idx = np.where(signs == +1)[0]
    minus_idx = np.where(signs == -1)[0]
    ax.scatter(gammas[plus_idx], np.ones(len(plus_idx)),
               marker='+', color='steelblue', s=50,
               label='+  (n=%d)' % n_plus)
    ax.scatter(gammas[minus_idx], -np.ones(len(minus_idx)),
               marker='_', color='darkorange', s=50,
               label='-  (n=%d)' % n_minus)
    ax.set_yticks([-1, +1])
    ax.set_yticklabels(['-', '+'])
    ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel('gamma_k')
    ax.set_title('Sign vs gamma (visual periodicity check)')
    ax.legend(loc='center right')
    ax.grid(True, alpha=0.3)

    # 2. Autocorrelation
    ax = axes[0, 1]
    ax.bar(range(1, N_LAGS + 1), autocorrs, color='steelblue',
           alpha=0.85, edgecolor='black', linewidth=0.4)
    ax.axhline(ci_95, color='red', linestyle='--', alpha=0.6,
               label='95%% white-noise band (+/- %.4f)' % ci_95)
    ax.axhline(-ci_95, color='red', linestyle='--', alpha=0.6)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xlabel('lag h')
    ax.set_ylabel('r_h')
    ax.set_title('Autocorrelation function (lags 1..%d)' % N_LAGS)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 3. Spacing histogram by sign
    ax = axes[1, 0]
    bins = np.linspace(spacings.min(), spacings.max(), 22)
    ax.hist(spacing_at_plus, bins=bins, alpha=0.6, color='steelblue',
            label='sign[k] = + (n=%d)' % len(spacing_at_plus),
            edgecolor='black', linewidth=0.4)
    ax.hist(spacing_at_minus, bins=bins, alpha=0.6, color='darkorange',
            label='sign[k] = - (n=%d)' % len(spacing_at_minus),
            edgecolor='black', linewidth=0.4)
    ax.set_xlabel('spacing  gamma_{k+1} - gamma_k')
    ax.set_ylabel('count')
    ax.set_title('Spacing by sign\n(Pearson sign-spacing = %+.3f)'
                 % pearson_sign_spacing)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Spacing histogram by sign-change
    ax = axes[1, 1]
    ax.hist(spacing_at_change, bins=bins, alpha=0.6, color='crimson',
            label='sign change (n=%d)' % len(spacing_at_change),
            edgecolor='black', linewidth=0.4)
    ax.hist(spacing_at_nochange, bins=bins, alpha=0.6, color='darkgreen',
            label='same sign (n=%d)' % len(spacing_at_nochange),
            edgecolor='black', linewidth=0.4)
    ax.set_xlabel('spacing  gamma_{k+1} - gamma_k')
    ax.set_ylabel('count')
    ax.set_title('Spacing by sign-change\n(Pearson change-spacing = %+.3f)'
                 % pearson_change_spacing)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=120)
    plt.close()
    print("Plot written: %s" % plot_path)
    print()

    # -----------------------------------------------------------------
    # Verdict (exploratory; report all four tests; no overall pass/fail)
    # -----------------------------------------------------------------
    print("=" * 76)
    print("Task 4 (sign distribution): exploratory, no pass/fail.")
    findings = []
    if p_bin < 0.05:
        findings.append("marginal split deviates from 50/50 (binomial p=%.4f)" % p_bin)
    if p_r < 0.05:
        findings.append("runs test rejects randomness (z=%+.3f, p=%.4f)" % (z_r, p_r))
    if n_outside > 2 * expected_outside:
        findings.append("autocorrelation: %d/%d lags outside white-noise band"
                        % (n_outside, N_LAGS))
    if abs(pearson_sign_spacing) > ci_corr:
        findings.append("sign-spacing correlation %+.3f exceeds CI" % pearson_sign_spacing)
    if abs(pearson_change_spacing) > ci_corr:
        findings.append("change-spacing correlation %+.3f exceeds CI"
                        % pearson_change_spacing)
    if findings:
        print("Findings (%d):" % len(findings))
        for f in findings:
            print("  - %s" % f)
    else:
        print("No findings: the sign sequence is consistent with i.i.d. fair-coin")
        print("draws across all four tests. Paper 196 records the rule-out.")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
