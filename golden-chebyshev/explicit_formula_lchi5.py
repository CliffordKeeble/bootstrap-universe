# explicit_formula_lchi5.py - Layer 2: spectral reconstruction of E(x) via L(chi_5) zeros
# Bootstrap Universe Programme - Paper 196 candidate
#
# Status flag: DERIVED.
#
# Standalone:
#   python explicit_formula_lchi5.py [--n-zeros 25] [--dps 30] [--sweep]
#
# Compares the explicit-formula prediction
#     E(x) ~ (sqrt(x) - psi(x, chi_5)) / log(x)
#     psi(x, chi_5) ~ -2 Re Sum_n  x^(1/2 + i gamma_n) / (1/2 + i gamma_n)
# against the Layer 1 direct prime count (chebyshev_bias_data.csv).
#
# Derivation (one line):
#   psi(x, chi) = theta(x, chi) + psi(sqrt(x), chi^2) + ...    [chi^2 = chi_0 for chi_5]
#   So theta(x, chi_5) ~ psi(x, chi_5) - sqrt(x).
#   pi(x, chi) = theta(x, chi)/log(x) + integral_2^x theta/(t log^2 t) dt   [partial summation]
#   E(x) = -pi(x, chi_5) -> E ~ (sqrt(x) - psi_pred)/log(x) at leading order.
# The dropped integral is O(sqrt(x)/log^2(x)); it shows up as a small
# +1/log(x)-scale systematic offset in normalised units. Reported, not absorbed.

import argparse
import csv
import time
from pathlib import Path

import numpy as np
import mpmath
from mpmath import mp, mpf

from lchi5_zeros import get_zeros


# ---------- Layer 1 CSV reader ----------

def load_layer1(csv_path):
    """Returns (xs, E_actual, E_norm_actual) as numpy arrays."""
    xs, e_actual, e_norm_actual = [], [], []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        next(reader)  # header
        for row in reader:
            xs.append(int(row[0]))
            e_actual.append(int(row[3]))
            e_norm_actual.append(float(row[4]))
    return (np.array(xs, dtype=np.int64),
            np.array(e_actual, dtype=np.int64),
            np.array(e_norm_actual, dtype=np.float64))


# ---------- Spectral prediction (vectorised float64) ----------
#
# After zeros are refined to many digits via lchi5_zeros.get_zeros, we cast
# them to float64 for the prediction loop. Each gamma is held to ~16 digits;
# at log(x) <= 18.5 (x <= 10^8) the phase gamma*log(x) is computed with
# ~14 effective digits, far better than the explicit-formula truncation
# error from using only N=25 zeros. Float64 is the right precision here.

def predict_psi_chi5(xs, gammas):
    """
    psi(x, chi_5) ~ -2 * sqrt(x) * Sum_n [(1/2)*cos(g_n log x) + g_n sin(g_n log x)] / (1/4 + g_n^2)
    xs, gammas: float arrays. Returns float64 array same shape as xs.
    """
    xs = xs.astype(np.float64)
    log_xs = np.log(xs)
    sqrt_xs = np.sqrt(xs)
    g = np.asarray(gammas, dtype=np.float64)
    angle = np.outer(log_xs, g)              # (Nx, Ngammas)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    denom = 0.25 + g * g                     # (Ngammas,)
    contrib = (0.5 * cos_a + g * sin_a) / denom
    return -2.0 * sqrt_xs * contrib.sum(axis=1)


def predict_E(xs, gammas):
    """E_pred(x) = (sqrt(x) - psi_pred(x, chi_5)) / log(x). Float64 array."""
    xs_f = xs.astype(np.float64)
    return (np.sqrt(xs_f) - predict_psi_chi5(xs, gammas)) / np.log(xs_f)


def normalise(xs, e):
    """E / (sqrt(x)/log(x)) = E * log(x) / sqrt(x). Granville-Martin units."""
    x = xs.astype(np.float64)
    return e * np.log(x) / np.sqrt(x)


# ---------- Diagnostics ----------

def summary_stats(xs, e_norm_actual, e_norm_pred, x_threshold=1e6):
    """Stats on samples with x > threshold (where leading-order should dominate)."""
    mask = xs > x_threshold
    if not mask.any():
        return None
    a = e_norm_actual[mask]
    p = e_norm_pred[mask]
    resid = p - a
    return {
        'n_samples': int(mask.sum()),
        'corr': float(np.corrcoef(a, p)[0, 1]) if len(a) > 1 else float('nan'),
        'mean_actual': float(np.mean(a)),
        'mean_pred': float(np.mean(p)),
        'mean_residual': float(np.mean(resid)),
        'rms_residual': float(np.sqrt(np.mean(resid**2))),
        'max_abs_residual': float(np.max(np.abs(resid))),
        'rms_actual': float(np.sqrt(np.mean(a**2))),
    }


def ratio_table(xs, e_actual, e_pred, e_norm_actual, e_norm_pred, decades=None):
    """One row per decade; nearest-x lookup. Helps surface scale errors."""
    if decades is None:
        decades = [10**k for k in range(2, 9)]
    rows = []
    for x_target in decades:
        idx = int(np.argmin(np.abs(xs - x_target)))
        x = int(xs[idx])
        a = int(e_actual[idx])
        p = float(e_pred[idx])
        ratio = p / a if a != 0 else float('inf')
        rows.append((x, a, p, ratio, e_norm_actual[idx], e_norm_pred[idx]))
    return rows


def offset_diag(xs, e_norm_actual, e_norm_pred):
    """
    Diagnose the systematic offset implied by the dropped partial-summation
    integral. That integral predicts an extra +c/log(x) bias on E_norm with
    c ~ 2 (leading order of integral_2^x dt/(sqrt(t) log^2 t)).
    Fit (e_norm_actual - e_norm_pred) ~ c / log(x); report c.
    """
    x = xs.astype(np.float64)
    inv_log = 1.0 / np.log(x)
    diff = e_norm_actual - e_norm_pred       # actual - pred = missing piece
    mask = x > 1e3                           # drop the very-low-x noise band
    A = inv_log[mask].reshape(-1, 1)
    b = diff[mask]
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = b - A.flatten() * c[0]
    return float(c[0]), float(np.sqrt(np.mean(resid**2)))


def n_sweep(xs, e_norm_actual, gammas_full, n_values=(5, 10, 15, 20, 25)):
    out = {}
    for N in n_values:
        if N > len(gammas_full):
            continue
        e_pred = predict_E(xs, gammas_full[:N])
        e_norm_pred = normalise(xs, e_pred)
        stats = summary_stats(xs, e_norm_actual, e_norm_pred)
        out[N] = (stats, e_norm_pred)
    return out


# ---------- Plotting ----------

def plot_overlay(xs, e_norm_actual, e_norm_pred, out_path, n_zeros):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print(f'  (matplotlib not available; skipping {out_path.name})')
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axhline(1, color='red', linewidth=0.5, linestyle=':',
               label='Chebyshev bias offset = 1')
    ax.semilogx(xs, e_norm_actual, color='steelblue', linewidth=0.9,
                label='direct count (Layer 1)')
    ax.semilogx(xs, e_norm_pred, color='crimson', linewidth=0.9,
                label=f'explicit formula (N={n_zeros} zeros)')
    ax.set_xlabel('x')
    ax.set_ylabel(r'$E(x) / (\sqrt{x}/\ln x)$')
    ax.set_title(r'Layer 2: explicit-formula reconstruction of $E(x)$ for $q=5$')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def plot_sweep(xs, e_norm_actual, sweep_results, out_path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.semilogx(xs, e_norm_actual, color='black', linewidth=1.4,
                label='direct count (Layer 1)')
    cmap = plt.get_cmap('viridis')
    Ns = sorted(sweep_results.keys())
    for i, N in enumerate(Ns):
        _, e_norm_pred = sweep_results[N]
        ax.semilogx(xs, e_norm_pred,
                    color=cmap(i / max(1, len(Ns) - 1)),
                    linewidth=0.7, alpha=0.9, label=f'N = {N}')
    ax.set_xlabel('x')
    ax.set_ylabel(r'$E(x) / (\sqrt{x}/\ln x)$')
    ax.set_title('Layer 2 sweep: spectral reconstruction sharpens with N')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description='Layer 2: explicit-formula reconstruction of E(x) via L(chi_5) zeros')
    parser.add_argument('--n-zeros', type=int, default=25,
                        help='Number of L(chi_5) zeros to use (default 25)')
    parser.add_argument('--dps', type=int, default=30,
                        help='mpmath precision for zero refinement (default 30)')
    parser.add_argument('--layer1-csv', type=str, default='chebyshev_bias_data.csv',
                        help='Path to Layer 1 CSV (relative to this script)')
    parser.add_argument('--sweep', action='store_true',
                        help='Also run the N-sweep diagnostic (N=5,10,15,20,25)')
    parser.add_argument('--outdir', type=str, default=None,
                        help='Output directory (default: alongside this script)')
    args = parser.parse_args()

    here = Path(__file__).parent
    outdir = Path(args.outdir) if args.outdir else here
    outdir.mkdir(parents=True, exist_ok=True)
    csv_in = here / args.layer1_csv

    print('=' * 70)
    print('LAYER 2: explicit-formula reconstruction of E(x) for q=5')
    print(f'  N zeros = {args.n_zeros}, dps = {args.dps}')
    print('=' * 70)

    print(f'\nLoading Layer 1 data from {csv_in.name}...')
    xs, e_actual, e_norm_actual = load_layer1(csv_in)
    print(f'  {len(xs)} x-samples, range [{xs[0]:,}, {xs[-1]:,}]')

    print(f'\nRefining {args.n_zeros} zeros of L(s, chi_5)...')
    t0 = time.time()
    gammas_mp = get_zeros(n=args.n_zeros, dps=args.dps, validate=True, verbose=True)
    gammas = np.array([float(g) for g in gammas_mp], dtype=np.float64)
    print(f'  done in {time.time() - t0:.1f}s; gamma_1={gammas[0]:.4f}, '
          f'gamma_N={gammas[-1]:.4f}')

    print(f'\nComputing E_pred at {len(xs)} sample points...')
    t0 = time.time()
    e_pred = predict_E(xs, gammas)
    e_norm_pred = normalise(xs, e_pred)
    print(f'  done in {time.time() - t0:.3f}s')

    # CSV output
    csv_out = outdir / 'explicit_formula_data.csv'
    with open(csv_out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['x', 'E_actual', 'E_predicted', 'E_norm_actual',
                    'E_norm_predicted', 'residual_norm'])
        for i in range(len(xs)):
            w.writerow([int(xs[i]), int(e_actual[i]), f'{e_pred[i]:.4f}',
                        f'{e_norm_actual[i]:.6f}', f'{e_norm_pred[i]:.6f}',
                        f'{e_norm_pred[i] - e_norm_actual[i]:.6f}'])
    print(f'\n  wrote {csv_out.name}')

    plot_overlay(xs, e_norm_actual, e_norm_pred,
                 outdir / 'explicit_formula.png', args.n_zeros)
    print(f'  wrote explicit_formula.png')

    # Summary stats (high-x window)
    stats = summary_stats(xs, e_norm_actual, e_norm_pred)
    print('\n--- Summary stats (x > 10^6) ---')
    if stats:
        print(f'  samples:             {stats["n_samples"]}')
        print(f'  Pearson correlation: {stats["corr"]:+.4f}')
        print(f'  mean E_norm actual:  {stats["mean_actual"]:+.4f}')
        print(f'  mean E_norm pred:    {stats["mean_pred"]:+.4f}')
        print(f'  mean residual:       {stats["mean_residual"]:+.4f}  (pred - actual)')
        print(f'  RMS residual:        {stats["rms_residual"]:.4f}')
        print(f'  max |residual|:      {stats["max_abs_residual"]:.4f}')
        print(f'  RMS actual:          {stats["rms_actual"]:.4f}')

    # Per-decade ratio table
    print('\n--- Ratio diagnostic (one row per decade) ---')
    print(f'  {"x":>12} {"E_act":>8} {"E_pred":>10} {"ratio":>8} '
          f'{"En_act":>9} {"En_pred":>9}')
    for x, a, p, ratio, na, np_ in ratio_table(xs, e_actual, e_pred,
                                               e_norm_actual, e_norm_pred):
        print(f'  {x:>12,} {a:>+8d} {p:>+10.2f} {ratio:>8.3f} '
              f'{na:>+9.4f} {np_:>+9.4f}')

    # Systematic offset fit: (actual - pred) ~ c / log(x)
    c_fit, fit_rms = offset_diag(xs, e_norm_actual, e_norm_pred)
    print('\n--- Systematic offset fit ---')
    print(f'  Model: (E_norm_actual - E_norm_pred) = c / log(x)')
    print(f'  c_fit = {c_fit:+.4f}   (theory predicts +2 from partial-summation integral)')
    print(f'  RMS residual after fit: {fit_rms:.4f}')

    # Sweep
    if args.sweep:
        print('\n--- N-sweep diagnostic ---')
        sweep = n_sweep(xs, e_norm_actual, gammas)
        print(f'  {"N":>3} {"corr":>8} {"mean_resid":>11} {"rms_resid":>10} '
              f'{"max|res|":>10}')
        for N in sorted(sweep.keys()):
            s = sweep[N][0]
            print(f'  {N:>3d} {s["corr"]:>+8.4f} {s["mean_residual"]:>+11.4f} '
                  f'{s["rms_residual"]:>10.4f} {s["max_abs_residual"]:>10.4f}')
        plot_sweep(xs, e_norm_actual, sweep, outdir / 'explicit_formula_sweep.png')
        print(f'  wrote explicit_formula_sweep.png')

    print('\n' + '=' * 70)
    print('LAYER 2 COMPLETE')
    print('=' * 70)


if __name__ == '__main__':
    main()
