# density_transition_robustness.py - robustness for density_transition_scan
# Bootstrap Universe Programme - density transition brief follow-up
#
# Standalone:
#   python density_transition_robustness.py
#
# Three robustness questions about the chi_5 envelope dip near tau=137:
#
#   Q1: Where does the tau=137 dip rank among LOCAL MINIMA of the predicted
#       envelope in [130, 150]?  If it is the deepest, that is content.  If
#       it is the median local minimum, the alignment is coincidental.
#
#   Q2: Window-width sensitivity. Does the dip survive across a range of
#       rolling-RMS window widths?
#
#   Q3: Truncation sensitivity. Does the dip survive as we vary the number
#       of L(chi_5) zeros from 10..25?
#
# Outputs:
#   density_transition_robustness.csv  - tabulated robustness sweep
#   density_transition_robustness.md   - inline summary

import csv
from pathlib import Path

import numpy as np

from lchi5_zeros import get_zeros
from density_lchi5 import LCHI4_ZEROS_25
from density_transition_scan import B_norm, rolling_rms

TAU_PROBES = {
    'tau_alpha (137)':   137.0,
    'tau_digit (140.5)': 140.5,
    'tau_e5 (142)':      142.0,
}


def find_local_minima(tau_grid, env, window_pts=5):
    """Return list of (tau, env) for local minima (env[i] is min in
    [i-w..i+w] for half-window w). Endpoints excluded."""
    w = window_pts // 2
    minima = []
    for i in range(w, len(env) - w):
        if env[i] == env[i - w:i + w + 1].min() and env[i] < env[i - 1]:
            minima.append((float(tau_grid[i]), float(env[i])))
    return minima


def find_local_maxima(tau_grid, env, window_pts=5):
    w = window_pts // 2
    maxima = []
    for i in range(w, len(env) - w):
        if env[i] == env[i - w:i + w + 1].max() and env[i] > env[i - 1]:
            maxima.append((float(tau_grid[i]), float(env[i])))
    return maxima


def compute_envelope_grid(tau_grid, gammas, window_tau, dtau):
    B = B_norm(tau_grid, gammas)
    window_pts = int(round(2.0 * window_tau / dtau))
    if window_pts % 2 == 0:
        window_pts += 1
    return rolling_rms(B, window_pts)


def main():
    here = Path(__file__).parent
    outdir = here

    # --- Get zeros ---
    print('Refining 25 chi_5 zeros at dps=30...')
    g5_mp = get_zeros(n=25, dps=30, validate=True, verbose=False)
    g5_full = np.array([float(z) for z in g5_mp], dtype=np.float64)
    g4_full = np.array(LCHI4_ZEROS_25, dtype=np.float64)

    dtau = 0.05
    tau_min, tau_max = 130.0, 150.0
    n_scan = int(round((tau_max - tau_min) / dtau)) + 1
    tau_scan = np.linspace(tau_min, tau_max, n_scan)

    print('=' * 72)
    print('ROBUSTNESS Q1: rank of programme-tau dips among local minima of envelope')
    print('=' * 72)
    print('  scan window: tau in [130, 150], dtau = 0.05, window_tau = 1.0\n')

    rows = []

    env_chi5 = compute_envelope_grid(tau_scan, g5_full, window_tau=1.0, dtau=dtau)
    env_chi4 = compute_envelope_grid(tau_scan, g4_full, window_tau=1.0, dtau=dtau)

    for chi_name, env in [('chi_5', env_chi5), ('chi_4', env_chi4)]:
        minima = find_local_minima(tau_scan, env, window_pts=21)  # 1.0 tau-units half-w
        maxima = find_local_maxima(tau_scan, env, window_pts=21)
        minima_sorted = sorted(minima, key=lambda x: x[1])         # deepest first
        print(f'\n  {chi_name}: {len(minima)} local minima, {len(maxima)} local maxima '
              f'in scan window')
        if minima_sorted:
            print(f'    deepest 5 local minima:')
            for i, (t, v) in enumerate(minima_sorted[:5], 1):
                print(f'      #{i}  tau = {t:6.2f}  env = {v:.4f}')

        for label, t_target in TAU_PROBES.items():
            # find nearest local minimum within +/- 1 tau-unit
            nearby = [(t, v) for t, v in minima if abs(t - t_target) <= 1.0]
            if not nearby:
                rank = None
                nearest = None
            else:
                nearest = min(nearby, key=lambda x: abs(x[0] - t_target))
                # rank of this minimum among all minima (1 = deepest)
                rank = next(i for i, m in enumerate(minima_sorted, 1)
                            if m == nearest)
            print(f'    nearest local min to {label}:', end=' ')
            if nearest is None:
                print('  none within +/- 1 tau-unit')
            else:
                print(f'tau = {nearest[0]:6.2f}, env = {nearest[1]:.4f}, '
                      f'rank {rank}/{len(minima)}')
            rows.append({
                'character': chi_name, 'probe': label,
                'tau_target': t_target,
                'tau_min_nearest': nearest[0] if nearest else None,
                'env_min': nearest[1] if nearest else None,
                'rank': rank,
                'n_minima': len(minima),
            })

    # --- Q2: window-width sensitivity for chi_5 dip near tau=137 ---
    print('\n' + '=' * 72)
    print('ROBUSTNESS Q2: window-width sensitivity (chi_5 envelope at programme tau)')
    print('=' * 72)
    window_widths = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    print(f'  {"w_tau":>6} {"env(137)":>10} {"env(140.5)":>11} {"env(142)":>10} '
          f'{"baseline_mean":>14} {"baseline_std":>13}')
    # baseline grid for window-sensitivity
    n_base = int(round((200 - 20) / dtau)) + 1
    tau_base = np.linspace(20.0, 200.0, n_base)
    win_rows = []
    for w in window_widths:
        env_scan = compute_envelope_grid(tau_scan, g5_full, w, dtau)
        env_base = compute_envelope_grid(tau_base, g5_full, w, dtau)
        m, s = float(env_base.mean()), float(env_base.std(ddof=1))
        vals = {}
        for label, t in TAU_PROBES.items():
            idx = int(np.argmin(np.abs(tau_scan - t)))
            vals[label] = float(env_scan[idx])
        print(f'  {w:>6.2f} {vals["tau_alpha (137)"]:>10.4f} '
              f'{vals["tau_digit (140.5)"]:>11.4f} {vals["tau_e5 (142)"]:>10.4f} '
              f'{m:>14.4f} {s:>13.4f}')
        win_rows.append({'window_tau': w,
                         'env_137': vals['tau_alpha (137)'],
                         'env_140p5': vals['tau_digit (140.5)'],
                         'env_142': vals['tau_e5 (142)'],
                         'baseline_mean': m,
                         'baseline_std': s,
                         'z_137': (vals['tau_alpha (137)'] - m) / s,
                         'z_140p5': (vals['tau_digit (140.5)'] - m) / s,
                         'z_142': (vals['tau_e5 (142)'] - m) / s})
    print(f'\n  z-scores at tau=137:')
    for r in win_rows:
        print(f"    w_tau = {r['window_tau']:.2f}: z = {r['z_137']:+.2f}")

    # --- Q3: truncation sensitivity ---
    print('\n' + '=' * 72)
    print('ROBUSTNESS Q3: zero-count truncation (chi_5 envelope at tau=137)')
    print('=' * 72)
    Ns = [10, 15, 20, 25]
    print(f'  {"N_zeros":>8} {"env(137)":>10} {"baseline_mean":>14} '
          f'{"baseline_std":>13} {"z_137":>9}')
    trunc_rows = []
    for N in Ns:
        g_trunc = g5_full[:N]
        env_scan = compute_envelope_grid(tau_scan, g_trunc, 1.0, dtau)
        env_base = compute_envelope_grid(tau_base, g_trunc, 1.0, dtau)
        m, s = float(env_base.mean()), float(env_base.std(ddof=1))
        idx = int(np.argmin(np.abs(tau_scan - 137.0)))
        v = float(env_scan[idx])
        z = (v - m) / s
        print(f'  {N:>8d} {v:>10.4f} {m:>14.4f} {s:>13.4f} {z:>+9.2f}')
        trunc_rows.append({'N_zeros': N, 'env_137': v,
                           'baseline_mean': m, 'baseline_std': s, 'z_137': z})

    # --- Save CSV ---
    csv_out = outdir / 'density_transition_robustness.csv'
    with open(csv_out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['# Q1: rank of nearest local minimum (within +/- 1 tau-unit) '
                    'among all local minima'])
        w.writerow(['Q', 'character', 'probe', 'tau_target', 'tau_min_nearest',
                    'env_min', 'rank', 'n_minima'])
        for r in rows:
            w.writerow(['Q1', r['character'], r['probe'], r['tau_target'],
                        r['tau_min_nearest'], r['env_min'],
                        r['rank'], r['n_minima']])
        w.writerow([])
        w.writerow(['# Q2: window-width sensitivity (chi_5)'])
        w.writerow(['Q', 'window_tau', 'env_137', 'env_140p5', 'env_142',
                    'baseline_mean', 'baseline_std', 'z_137',
                    'z_140p5', 'z_142'])
        for r in win_rows:
            w.writerow(['Q2', r['window_tau'], r['env_137'],
                        r['env_140p5'], r['env_142'],
                        r['baseline_mean'], r['baseline_std'],
                        r['z_137'], r['z_140p5'], r['z_142']])
        w.writerow([])
        w.writerow(['# Q3: zero-count truncation (chi_5)'])
        w.writerow(['Q', 'N_zeros', 'env_137', 'baseline_mean',
                    'baseline_std', 'z_137'])
        for r in trunc_rows:
            w.writerow(['Q3', r['N_zeros'], r['env_137'],
                        r['baseline_mean'], r['baseline_std'], r['z_137']])
    print(f'\n  wrote {csv_out.name}')

    print('\n' + '=' * 72)
    print('ROBUSTNESS COMPLETE')
    print('=' * 72)


if __name__ == '__main__':
    main()
