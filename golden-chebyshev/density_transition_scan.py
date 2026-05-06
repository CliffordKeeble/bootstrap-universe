# density_transition_scan.py - explicit-formula scan over tau in [130, 150]
# Bootstrap Universe Programme - density transition brief (5 May 2026)
#
# Standalone:
#   python density_transition_scan.py [--tau-min 130] [--tau-max 150]
#
# Computes the normalised Chebyshev bias
#     B_chi(x) / sqrt(x) = psi_pred(x, chi) / sqrt(x)
#                       = -2 Sum_n [(1/2) cos(g_n tau) + g_n sin(g_n tau)] / (1/4 + g_n^2)
# at tau = ln x, for chi_5 and chi_4 using the first 25 zeros of L(s, chi).
#
# Status: this is the EXPLICIT-FORMULA PREDICTION over tau in [130, 150]. No
# direct prime enumeration is possible at these scales (x in [e^130, e^150]).
# The curve is by construction a smooth quasi-periodic function in tau; the
# only "fingerprints" detectable here are envelope/beat features arising from
# interference of the lowest L(s, chi) zeros. Such features are statements
# about the L-function zero spectrum, not about a prime-distribution
# transition (no observation possible at these tau).
#
# Programme-relevant tau:
#   tau_alpha = 137    (PNT density 1/ln N = 1/137 = alpha)
#   tau_digit = 140.5  (10^61 ~ e^140.5, first digital boundary)
#   tau_e5    = 142    (e^5 past alpha-density, closure scale)
#
# Diagnostics:
#   - Rolling-RMS envelope of B_norm over scan and baseline windows
#   - Baseline mean / std for envelope amplitude
#   - z-scores at programme-relevant tau (descriptive, NOT a strict p-value:
#     the baseline envelope is quasi-periodic, not Gaussian)
#   - Top-K envelope extrema in scan window
#   - Cross-character (chi_4) for comparison

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from lchi5_zeros import get_zeros
from density_lchi5 import LCHI4_ZEROS_25


# ---------- Programme-relevant tau ----------

TAU_ALPHA = 137.0       # PNT density 1/ln N = 1/137
TAU_DIGIT = 140.5       # 10^61 ~ e^140.5 first digital boundary
TAU_E5    = 142.0       # e^5 past alpha-density

PROGRAMME_TAU = [
    (TAU_ALPHA, 'tau_alpha (e^137, alpha-density)'),
    (TAU_DIGIT, 'tau_digit (10^61 ~ e^140.5)'),
    (TAU_E5,    'tau_e5    (e^142, closure scale)'),
]


# ---------- Spectral prediction ----------

def B_norm(tau, gammas):
    """
    Normalised Chebyshev bias amplitude:
        B_chi(x) / sqrt(x) at x = e^tau, leading explicit-formula term
        B_norm(tau) = -2 Sum_n [(1/2) cos(g_n tau) + g_n sin(g_n tau)] / (1/4 + g_n^2)
    Real-character sign convention follows psi_pred / sqrt(x) as in
    explicit_formula_lchi5.predict_psi_chi5.

    tau    : array-like, shape (M,) or scalar.
    gammas : array-like, shape (N,), positive imaginary parts of zeros.
    Returns: array shape (M,).
    """
    tau_arr = np.atleast_1d(tau).astype(np.float64)
    g = np.asarray(gammas, dtype=np.float64)
    angle = np.outer(tau_arr, g)                       # (M, N)
    contrib = (0.5 * np.cos(angle) + g * np.sin(angle)) / (0.25 + g * g)
    return -2.0 * contrib.sum(axis=1)


def rolling_rms(y, window_pts):
    """
    RMS in a centred window of `window_pts` samples (must be odd).
    Reflective edge padding so boundary samples aren't biased low.
    """
    if window_pts <= 1:
        return np.abs(y).copy()
    if window_pts % 2 == 0:
        raise ValueError('window_pts must be odd')
    pad = window_pts // 2
    y_pad = np.pad(y, pad, mode='reflect')
    kernel = np.ones(window_pts, dtype=np.float64) / window_pts
    rms_sq = np.convolve(y_pad ** 2, kernel, mode='valid')
    return np.sqrt(rms_sq)


# ---------- Probe diagnostics ----------

def envelope_at(tau_grid, env, tau_target):
    """Envelope value at the grid point nearest tau_target."""
    if not (tau_grid[0] <= tau_target <= tau_grid[-1]):
        return None
    idx = int(np.argmin(np.abs(tau_grid - tau_target)))
    return idx, float(tau_grid[idx]), float(env[idx])


def local_stats(tau_grid, env, tau_target, dtau_local):
    """Local-window stats around tau_target."""
    mask = np.abs(tau_grid - tau_target) <= dtau_local
    if not mask.any():
        return None
    band = env[mask]
    return {
        'n': int(mask.sum()),
        'mean': float(band.mean()),
        'std': float(band.std(ddof=1)) if mask.sum() > 1 else float('nan'),
        'min': float(band.min()),
        'max': float(band.max()),
    }


def top_k_extrema(tau_grid, env, k=5):
    """Top-k envelope peaks and bottom-k troughs (by value)."""
    idx_sort = np.argsort(env)
    troughs = [(float(tau_grid[i]), float(env[i])) for i in idx_sort[:k]]
    peaks = [(float(tau_grid[i]), float(env[i])) for i in idx_sort[-k:][::-1]]
    return peaks, troughs


# ---------- Plotting ----------

def plot_scan(tau_scan, B5_scan, B4_scan, env5_scan, env4_scan,
              m5, s5, m4, s4, args, outdir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    # Panel 1: B_norm raw oscillations
    ax = axes[0]
    ax.axhline(0, color='gray', lw=0.5)
    ax.plot(tau_scan, B5_scan, color='steelblue', lw=0.8,
            label=r'$\chi_5$: $B_{\chi_5}(x)/\sqrt{x}$')
    ax.plot(tau_scan, B4_scan, color='crimson', lw=0.8, alpha=0.7,
            label=r'$\chi_4$: $B_{\chi_4}(x)/\sqrt{x}$')
    for tt, _ in PROGRAMME_TAU:
        if args.tau_min <= tt <= args.tau_max:
            ax.axvline(tt, color='k', ls=':', alpha=0.4)
    ax.set_ylabel(r'$B_\chi(x)/\sqrt{x}$')
    ax.set_title(r'Chebyshev bias (normalised), $\tau = \ln x \in '
                 fr'[{args.tau_min}, {args.tau_max}]$, first {args.n_zeros} zeros')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: rolling-RMS envelope vs baseline mean/2-sigma band
    ax = axes[1]
    ax.plot(tau_scan, env5_scan, color='steelblue', lw=1.2,
            label=fr'$\chi_5$ env (RMS $\pm{args.window_tau}$)')
    ax.plot(tau_scan, env4_scan, color='crimson', lw=1.2, alpha=0.85,
            label=fr'$\chi_4$ env')
    ax.axhline(m5, color='steelblue', ls='--', lw=0.8, alpha=0.7)
    ax.fill_between(tau_scan, m5 - 2 * s5, m5 + 2 * s5, color='steelblue',
                    alpha=0.10, label=r'$\chi_5$ baseline mean $\pm 2\sigma$')
    ax.axhline(m4, color='crimson', ls='--', lw=0.8, alpha=0.7)
    ax.fill_between(tau_scan, m4 - 2 * s4, m4 + 2 * s4, color='crimson',
                    alpha=0.08, label=r'$\chi_4$ baseline mean $\pm 2\sigma$')
    for tt, lbl in PROGRAMME_TAU:
        if args.tau_min <= tt <= args.tau_max:
            ax.axvline(tt, color='k', ls=':', alpha=0.5)
            ax.text(tt, ax.get_ylim()[1] * 0.97,
                    lbl.split(' ')[0],
                    rotation=90, va='top', ha='right', fontsize=8, alpha=0.7)
    ax.set_xlabel(r'$\tau = \ln x$')
    ax.set_ylabel('envelope (RMS)')
    ax.set_title(fr'Envelope amplitude, baseline = $\tau \in '
                 fr'[{args.baseline_min}, {args.baseline_max}]$')
    ax.legend(loc='best', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out = outdir / 'density_transition_scan.png'
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_baseline_context(tau_base, env5_base, env4_base, m5, m4, args, outdir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(tau_base, env5_base, color='steelblue', lw=0.7, alpha=0.85,
            label=r'$\chi_5$ envelope')
    ax.plot(tau_base, env4_base, color='crimson', lw=0.7, alpha=0.7,
            label=r'$\chi_4$ envelope')
    ax.axhline(m5, color='steelblue', ls='--', lw=0.8, alpha=0.7)
    ax.axhline(m4, color='crimson', ls='--', lw=0.8, alpha=0.7)
    ax.axvspan(args.tau_min, args.tau_max, color='gray', alpha=0.18,
               label='scan window [130, 150]')
    for tt, _ in PROGRAMME_TAU:
        ax.axvline(tt, color='k', ls=':', alpha=0.4)
    ax.set_xlabel(r'$\tau = \ln x$')
    ax.set_ylabel('envelope (RMS)')
    ax.set_title(fr'Envelope across baseline window $\tau \in '
                 fr'[{args.baseline_min}, {args.baseline_max}]$')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = outdir / 'density_transition_baseline.png'
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description='Explicit-formula scan over tau in [130, 150] for '
                    'chi_5 and chi_4, with envelope diagnostics.')
    parser.add_argument('--tau-min', type=float, default=130.0)
    parser.add_argument('--tau-max', type=float, default=150.0)
    parser.add_argument('--dtau', type=float, default=0.05,
                        help='Sample spacing in tau (default 0.05)')
    parser.add_argument('--baseline-min', type=float, default=20.0)
    parser.add_argument('--baseline-max', type=float, default=200.0)
    parser.add_argument('--window-tau', type=float, default=1.0,
                        help='Rolling-RMS half-width in tau-units (default 1.0)')
    parser.add_argument('--probe-tau', type=float, default=2.0,
                        help='Probe local-stats half-width (default 2.0)')
    parser.add_argument('--n-zeros', type=int, default=25)
    parser.add_argument('--dps', type=int, default=30)
    parser.add_argument('--outdir', type=str, default=None)
    args = parser.parse_args()

    here = Path(__file__).parent
    outdir = Path(args.outdir) if args.outdir else here
    outdir.mkdir(parents=True, exist_ok=True)

    print('=' * 72)
    print('DENSITY TRANSITION SCAN  (Mr Code brief, 5 May 2026)')
    print('  Explicit-formula prediction of B_chi(x)/sqrt(x) over tau = ln x')
    print(f'  Scan window:     tau in [{args.tau_min}, {args.tau_max}], '
          f'dtau = {args.dtau}')
    print(f'  Baseline window: tau in [{args.baseline_min}, {args.baseline_max}]')
    print(f'  N zeros = {args.n_zeros} per character')
    print('=' * 72)

    # --- Zeros ---
    print('\n--- Zero data ---')
    print(f'  chi_5: refining {args.n_zeros} zeros at dps = {args.dps}...')
    t0 = time.time()
    g5_mp = get_zeros(n=args.n_zeros, dps=args.dps, validate=True, verbose=False)
    g5 = np.array([float(z) for z in g5_mp], dtype=np.float64)
    print(f'    done in {time.time() - t0:.1f}s; '
          f'gamma_1 = {g5[0]:.4f}, gamma_N = {g5[-1]:.4f}')

    g4 = np.array(LCHI4_ZEROS_25, dtype=np.float64)
    if args.n_zeros < len(g4):
        g4 = g4[:args.n_zeros]
    print(f'  chi_4: using {len(g4)} LMFDB zeros at 4-decimal precision')
    print(f'    gamma_1 = {g4[0]:.4f}, gamma_N = {g4[-1]:.4f}')

    # --- Tau grids ---
    n_scan = int(round((args.tau_max - args.tau_min) / args.dtau)) + 1
    tau_scan = np.linspace(args.tau_min, args.tau_max, n_scan)

    n_base = int(round((args.baseline_max - args.baseline_min) / args.dtau)) + 1
    tau_base = np.linspace(args.baseline_min, args.baseline_max, n_base)

    print(f'\n--- Grids ---')
    print(f'  scan:     {n_scan} samples over tau in '
          f'[{args.tau_min}, {args.tau_max}]')
    print(f'  baseline: {n_base} samples over tau in '
          f'[{args.baseline_min}, {args.baseline_max}]')

    # --- Compute B_norm ---
    print('\n--- Computing B_norm = B_chi(x)/sqrt(x) ---')
    t0 = time.time()
    B5_scan = B_norm(tau_scan, g5)
    B5_base = B_norm(tau_base, g5)
    B4_scan = B_norm(tau_scan, g4)
    B4_base = B_norm(tau_base, g4)
    print(f'  done in {time.time() - t0:.2f}s')

    # --- Envelope (rolling RMS) ---
    window_pts = int(round(2.0 * args.window_tau / args.dtau))
    if window_pts % 2 == 0:
        window_pts += 1
    print(f'\n--- Rolling-RMS envelope ---')
    print(f'  window full-width = {2 * args.window_tau} tau-units = '
          f'{window_pts} samples')

    env5_scan = rolling_rms(B5_scan, window_pts)
    env5_base = rolling_rms(B5_base, window_pts)
    env4_scan = rolling_rms(B4_scan, window_pts)
    env4_base = rolling_rms(B4_base, window_pts)

    # --- Baseline statistics ---
    m5 = float(env5_base.mean())
    s5 = float(env5_base.std(ddof=1))
    m4 = float(env4_base.mean())
    s4 = float(env4_base.std(ddof=1))

    print(f'\n--- Baseline envelope statistics (tau in '
          f'[{args.baseline_min}, {args.baseline_max}]) ---')
    print(f'  chi_5: mean = {m5:.4f}, std = {s5:.4f}, '
          f'range [{env5_base.min():.4f}, {env5_base.max():.4f}]')
    print(f'  chi_4: mean = {m4:.4f}, std = {s4:.4f}, '
          f'range [{env4_base.min():.4f}, {env4_base.max():.4f}]')
    print(f'  note: baseline excludes scan window? {"yes" if args.baseline_max < args.tau_min or args.baseline_min > args.tau_max else "no (overlaps)"}')

    # --- Probe at programme-relevant tau ---
    print(f'\n--- Probe at programme-relevant tau ---')
    print(f'  z = (env(tau) - baseline_mean) / baseline_std')
    print(f'  z is DESCRIPTIVE: baseline envelope is quasi-periodic, not Gaussian.')

    probe_records = []   # for CSV
    for chi_name, env_scan, m, s in [
        ('chi_5', env5_scan, m5, s5),
        ('chi_4', env4_scan, m4, s4),
    ]:
        print(f'\n  {chi_name}:')
        for t_target, label in PROGRAMME_TAU:
            r = envelope_at(tau_scan, env_scan, t_target)
            if r is None:
                continue
            idx, t_actual, val = r
            z = (val - m) / s if s > 0 else float('nan')
            ls = local_stats(tau_scan, env_scan, t_target, args.probe_tau)
            print(f'    {label:<32s}: env = {val:.4f}, z = {z:+.2f}  '
                  f'(local-{args.probe_tau} mean = {ls["mean"]:.4f}, '
                  f'range [{ls["min"]:.4f}, {ls["max"]:.4f}])')
            probe_records.append({
                'character': chi_name,
                'tau_target': t_target,
                'tau_actual': t_actual,
                'env': val,
                'baseline_mean': m,
                'baseline_std': s,
                'z': z,
                'local_mean': ls['mean'],
                'local_std': ls['std'],
                'local_min': ls['min'],
                'local_max': ls['max'],
            })

    # --- Top-K extrema in scan window ---
    print(f'\n--- Envelope extrema in scan window [{args.tau_min}, {args.tau_max}] ---')
    for chi_name, env_scan in [('chi_5', env5_scan), ('chi_4', env4_scan)]:
        peaks, troughs = top_k_extrema(tau_scan, env_scan, k=5)
        print(f'  {chi_name}: top-5 peaks    {[(round(t, 2), round(v, 3)) for t, v in peaks]}')
        print(f'  {chi_name}: top-5 troughs  {[(round(t, 2), round(v, 3)) for t, v in troughs]}')

    # --- Compare envelope at programme-tau across the FULL baseline ---
    # How extreme is the scan-window envelope at programme tau, vs the entire
    # baseline tau distribution?  Quantile-based, distribution-free.
    print(f'\n--- Quantile of programme-tau envelope vs baseline distribution ---')
    print(f'  rank in baseline = fraction of baseline samples with env <= probe value')
    for chi_name, env_scan, env_base in [
        ('chi_5', env5_scan, env5_base),
        ('chi_4', env4_scan, env4_base),
    ]:
        print(f'\n  {chi_name}:')
        for t_target, label in PROGRAMME_TAU:
            r = envelope_at(tau_scan, env_scan, t_target)
            if r is None:
                continue
            _, _, val = r
            rank = float(np.mean(env_base <= val))
            print(f'    {label:<32s}: env = {val:.4f}, '
                  f'baseline-rank = {rank:.3f}  '
                  f'({"extreme" if rank < 0.05 or rank > 0.95 else "typical"})')

    # --- CSVs ---
    csv_scan = outdir / 'density_transition_scan.csv'
    with open(csv_scan, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tau', 'B_norm_chi5', 'env_rms_chi5',
                    'B_norm_chi4', 'env_rms_chi4'])
        for t, b5, e5, b4, e4 in zip(tau_scan, B5_scan, env5_scan,
                                     B4_scan, env4_scan):
            w.writerow([f'{t:.4f}', f'{b5:.6f}', f'{e5:.6f}',
                        f'{b4:.6f}', f'{e4:.6f}'])
    print(f'\n  wrote {csv_scan.name}')

    csv_base = outdir / 'density_transition_baseline.csv'
    with open(csv_base, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tau', 'B_norm_chi5', 'env_rms_chi5',
                    'B_norm_chi4', 'env_rms_chi4'])
        for t, b5, e5, b4, e4 in zip(tau_base, B5_base, env5_base,
                                     B4_base, env4_base):
            w.writerow([f'{t:.4f}', f'{b5:.6f}', f'{e5:.6f}',
                        f'{b4:.6f}', f'{e4:.6f}'])
    print(f'  wrote {csv_base.name}')

    csv_probe = outdir / 'density_transition_probes.csv'
    with open(csv_probe, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['character', 'tau_target', 'tau_actual', 'env',
                    'baseline_mean', 'baseline_std', 'z',
                    'local_mean', 'local_std', 'local_min', 'local_max'])
        for r in probe_records:
            w.writerow([r['character'], f'{r["tau_target"]:.4f}',
                        f'{r["tau_actual"]:.4f}', f'{r["env"]:.6f}',
                        f'{r["baseline_mean"]:.6f}', f'{r["baseline_std"]:.6f}',
                        f'{r["z"]:+.4f}', f'{r["local_mean"]:.6f}',
                        f'{r["local_std"]:.6f}', f'{r["local_min"]:.6f}',
                        f'{r["local_max"]:.6f}'])
    print(f'  wrote {csv_probe.name}')

    # --- Plots ---
    p1 = plot_scan(tau_scan, B5_scan, B4_scan, env5_scan, env4_scan,
                   m5, s5, m4, s4, args, outdir)
    if p1:
        print(f'  wrote {p1.name}')
    p2 = plot_baseline_context(tau_base, env5_base, env4_base, m5, m4,
                               args, outdir)
    if p2:
        print(f'  wrote {p2.name}')

    print('\n' + '=' * 72)
    print('SCAN COMPLETE')
    print('=' * 72)


if __name__ == '__main__':
    main()
