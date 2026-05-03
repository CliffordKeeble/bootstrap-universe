# off_line_test.py - Lock-breaking confirmation off the critical line
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a follow-up Task 2
#
# Per Mr Code Brief Addendum 3 (3 May 2026, CinC), Task 2.
#
# Purpose
# -------
# Promote Paper 125 sec.5's FRAMEWORK claim ("the golden phase lock breaks
# off the critical line") to OBSERVED, by computing the phase residual
# arg(P(chi_2)/P(chi_3)) - {nearest locked value} at three Re(s) values for
# the first 10 zeta zeros:
#   - s = 0.4 + i*gamma_k  (off-line, lower)
#   - s = 0.5 + i*gamma_k  (on-line baseline; should match phase_lock_test.py)
#   - s = 0.6 + i*gamma_k  (off-line, upper)
#
# Pass condition (from brief): for at least 9 of the 10 zeros, the residuals
# at Re(s) = 0.4 AND Re(s) = 0.6 both exceed BREAK_THRESHOLD = 0.01.
# Expected: 10/10 with off-line residuals of order 0.1-1.5 rad (much larger
# than the dps=50 on-line floor of ~10^-48). Average residual for a
# uniformly random phase from the nearest locked value is pi/4 ~ 0.79.
#
# Output: off_line_results.csv with per-zero residuals at the three Re(s);
# off_line_residuals.png with two panels (residual vs Re(s) showing the
# U-shape lock-breaking, and residual vs |Re(s)-0.5| testing left/right
# symmetry of the off-line behaviour).

import csv
import sys

import mpmath as mp
from mpmath import mp as _mp_ctx, mpc, mpf, fabs

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from character_projections import P
from phase_lock_test import locked_targets, classify_phase

N_ZEROS = 10
RE_VALUES = [mpf("0.4"), mpf("0.5"), mpf("0.6")]
BREAK_THRESHOLD = 0.01  # residual above this == lock decisively broken


def compute_residual_at(s):
    """Compute (phase, magnitude, direction, residual) at complex s."""
    p2 = P(s, 2)
    p3 = P(s, 3)
    ratio = p2 / p3
    phase = mp.arg(ratio)
    target_minus, target_plus = locked_targets()
    direction, residual = classify_phase(phase, target_minus, target_plus)
    return phase, fabs(ratio), direction, residual


def main():
    _mp_ctx.dps = 50

    print("off_line_test - lock-breaking confirmation off the critical line")
    print("=" * 76)
    print("dps = %d, N_ZEROS = %d, Re(s) in %s, break threshold = %s"
          % (_mp_ctx.dps, N_ZEROS,
             [float(re) for re in RE_VALUES], BREAK_THRESHOLD))
    print()

    print("Per-zero residuals:")
    print()
    print("|  k | gamma_k          | res(Re=0.4)  | res(Re=0.5)  | res(Re=0.6)  | broken? |")
    print("|----|------------------|--------------|--------------|--------------|---------|")

    results = []
    n_broken = 0

    for k in range(1, N_ZEROS + 1):
        rho = _mp_ctx.zetazero(k)
        gamma_k = rho.imag
        per_re = {}
        for re_val in RE_VALUES:
            s = mpc(re_val, gamma_k)
            phase, mag, direction, residual = compute_residual_at(s)
            per_re[float(re_val)] = {
                'phase': phase, 'magnitude': mag,
                'direction': direction, 'residual': residual,
            }
        results.append({'k': k, 'gamma': gamma_k, 'per_re': per_re})

        r04 = per_re[0.4]['residual']
        r05 = per_re[0.5]['residual']
        r06 = per_re[0.6]['residual']
        broken = (r04 > BREAK_THRESHOLD and r06 > BREAK_THRESHOLD)
        if broken:
            n_broken += 1

        print("| %2d | %16.12f | %12s | %12s | %12s | %7s |"
              % (k, float(gamma_k),
                 mp.nstr(r04, 4), mp.nstr(r05, 4), mp.nstr(r06, 4),
                 'YES' if broken else 'NO'))

    print()

    # ---- Summary statistics --------------------------------------------
    on_line_residuals = [r['per_re'][0.5]['residual'] for r in results]
    off_line_residuals_04 = [r['per_re'][0.4]['residual'] for r in results]
    off_line_residuals_06 = [r['per_re'][0.6]['residual'] for r in results]
    all_off_line = off_line_residuals_04 + off_line_residuals_06

    on_line_max = max(on_line_residuals)
    off_line_min = min(all_off_line)
    off_line_max = max(all_off_line)
    off_line_mean = sum(all_off_line) / mpf(len(all_off_line))

    # Symmetry check: are residuals at Re=0.4 close to those at Re=0.6 per zero?
    sym_diffs = [fabs(r04 - r06) for r04, r06 in
                 zip(off_line_residuals_04, off_line_residuals_06)]
    sym_max_diff = max(sym_diffs)
    sym_max_rel = max(fabs(r04 - r06) / ((r04 + r06) / 2)
                       for r04, r06 in
                       zip(off_line_residuals_04, off_line_residuals_06))

    print("Summary:")
    print("  Zeros where lock breaks at Re=0.4 AND Re=0.6: %d / %d"
          % (n_broken, N_ZEROS))
    print("  Max residual on-line  (Re=0.5):  %s   [should be at dps=50 floor]"
          % mp.nstr(on_line_max, 4))
    print("  Min residual off-line (Re in {0.4, 0.6}): %s   [should >> %s]"
          % (mp.nstr(off_line_min, 4), BREAK_THRESHOLD))
    print("  Max residual off-line (Re in {0.4, 0.6}): %s"
          % mp.nstr(off_line_max, 4))
    print("  Mean residual off-line:                  %s   [pi/4 ~ 0.785 expected]"
          % mp.nstr(off_line_mean, 4))
    print("  Off-line / on-line ratio (orders of magnitude): %s"
          % mp.nstr(off_line_min / on_line_max, 4))
    print()
    print("  Symmetry check (per zero, |res(Re=0.4) - res(Re=0.6)|):")
    print("    Max absolute diff:  %s" % mp.nstr(sym_max_diff, 4))
    print("    Max relative diff:  %s   [perfect symmetry would be 0]"
          % mp.nstr(sym_max_rel, 4))
    print()

    # ---- CSV output -----------------------------------------------------
    csv_path = 'off_line_results.csv'
    with open(csv_path, 'w', newline='') as f:
        f.write("# Off-line lock-breaking test (Layer 3a follow-up Task 2)\n")
        f.write("# dps=%d  N_ZEROS=%d  Re(s) in {0.4, 0.5, 0.6}\n"
                % (_mp_ctx.dps, N_ZEROS))
        f.write("# Pass: residual > %s at Re(s) in {0.4, 0.6} for >= %d zeros\n"
                % (BREAK_THRESHOLD, N_ZEROS - 1))
        w = csv.writer(f)
        w.writerow(['k', 'gamma', 'residual_at_0.4', 'residual_at_0.5',
                    'residual_at_0.6', 'magnitude_at_0.4',
                    'magnitude_at_0.5', 'magnitude_at_0.6', 'broken'])
        for r in results:
            r04 = r['per_re'][0.4]['residual']
            r05 = r['per_re'][0.5]['residual']
            r06 = r['per_re'][0.6]['residual']
            m04 = r['per_re'][0.4]['magnitude']
            m05 = r['per_re'][0.5]['magnitude']
            m06 = r['per_re'][0.6]['magnitude']
            broken = (r04 > BREAK_THRESHOLD and r06 > BREAK_THRESHOLD)
            w.writerow([
                r['k'], mp.nstr(r['gamma'], 30),
                mp.nstr(r04, 6), mp.nstr(r05, 6), mp.nstr(r06, 6),
                mp.nstr(m04, 6), mp.nstr(m05, 6), mp.nstr(m06, 6),
                'true' if broken else 'false',
            ])
    print("CSV written: %s  (%d rows + 3 header lines)" % (csv_path, N_ZEROS))

    # ---- Plot -----------------------------------------------------------
    # Two panels:
    #   Left  - residual vs Re(s), 10 U-shape lines, log-y. Visualises the
    #           lock-breaking continuum: low at 0.5, high at 0.4 and 0.6.
    #   Right - residual vs |Re(s) - 0.5|, scatter coloured by zero. Tests
    #           left/right symmetry of off-line residuals (expected: pairs
    #           at x=0.1 cluster vertically per zero).
    plot_path = 'off_line_residuals.png'
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    re_floats = [float(re) for re in RE_VALUES]
    cmap = plt.cm.viridis

    ax = axes[0]
    for i, r in enumerate(results):
        residuals = [max(float(r['per_re'][re]['residual']), 1e-50)
                     for re in re_floats]
        color = cmap(i / max(N_ZEROS - 1, 1))
        ax.semilogy(re_floats, residuals, '-o',
                    color=color, label='rho_%d' % r['k'], alpha=0.85)
    ax.axhline(BREAK_THRESHOLD, color='red', linestyle='--', alpha=0.6,
               label='break threshold = 0.01')
    ax.set_xlabel('Re(s)')
    ax.set_ylabel('phase residual (radians, log)')
    ax.set_title('Phase lock vs Re(s) at first 10 zeta zeros\n'
                 '(U-shape: locked on-line, broken off-line)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='center right', ncol=2)

    ax = axes[1]
    distances = [abs(float(re) - 0.5) for re in RE_VALUES]
    for i, r in enumerate(results):
        residuals = [max(float(r['per_re'][re]['residual']), 1e-50)
                     for re in re_floats]
        color = cmap(i / max(N_ZEROS - 1, 1))
        ax.semilogy(distances, residuals, 'o',
                    color=color, label='rho_%d' % r['k'], alpha=0.85,
                    markersize=8)
    ax.axhline(BREAK_THRESHOLD, color='red', linestyle='--', alpha=0.6)
    ax.set_xlabel('|Re(s) - 0.5|')
    ax.set_ylabel('phase residual (radians, log)')
    ax.set_title('Residual vs distance from critical line\n'
                 '(symmetry test: pairs at x=0.1 stack per zero)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=120)
    plt.close()
    print("Plot written: %s" % plot_path)
    print()

    # ---- Verdict --------------------------------------------------------
    print("=" * 76)
    overall_pass = n_broken >= (N_ZEROS - 1)
    print("Task 2 (off-line lock-breaking): %s"
          % ("PASS - lock breaks decisively off-line; framework claim promoted to OBSERVED"
             if overall_pass
             else "FAIL - lock-breaking weaker than expected, investigate"))


if __name__ == "__main__":
    main()
