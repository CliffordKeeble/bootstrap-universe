# density_lchi5.py - Chebyshev-bias logarithmic density for q=5, three paths.
# Bootstrap Universe Programme - Paper 196 candidate.
#
# Standalone:
#   python density_lchi5.py [--n-zeros 25] [--n-mc 1000000] [--seed 20260426]
#
# Computes the logarithmic-measure density
#     delta(5) := lim_{X->inf} (1/ln X) * mu_log{ x <= X : E(x) > 0 }
# i.e. the long-run logarithmic-measure fraction of x at which the q=5
# Chebyshev race shows the stubborn lead E(x) > 0.
#
# Three independent paths:
#
#   Path 1 (OBSERVED) - empirical from Layer 1 direct count.
#     Samples are log-spaced; the natural finite-X analogue of the RS
#     density is the logarithmic-measure fraction over [x_min, X], which
#     reduces (for log-uniform samples) to the simple sample fraction.
#
#   Path 2 (DERIVED, under GRH + LI) - Rubinstein-Sarnak construction.
#     Y = c + 2 * Sum_n cos(theta_n) / |rho_n|
#     theta_n iid uniform [0, 2pi], c = 1 (squares-vs-non-squares offset),
#     |rho_n| = sqrt(1/4 + gamma_n^2). delta = P(Y > 0) by Monte Carlo.
#
#   Path 3 (OBSERVED if found) - LMFDB / literature lookup.
#
# Reports all three. Convergence (or lack of it) is the result.

import argparse
import csv
import time
from pathlib import Path

import numpy as np

from lchi5_zeros import get_zeros


# ---------- Path 1: empirical density from Layer 1 CSV ----------

def _load_layer1(csv_path):
    xs, e_actual = [], []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            xs.append(int(row[0]))
            e_actual.append(int(row[3]))
    return np.array(xs, dtype=np.int64), np.array(e_actual, dtype=np.int64)


def path1_empirical(layer1_csv_path):
    """
    Empirical logarithmic-measure density of {E(x) > 0} from the Layer 1
    log-spaced sample. Returns a dict with two normalisations:

      delta_local      = (1/ln(X/x_min)) * mu_log{x in [x_min, X]: E(x)>0}
                       = local logarithmic-measure fraction over the sample
                         interval; the right finite-X analogue of the RS
                         density (best comparator to Paths 2/3).

      delta_rs_strict  = (1/ln X) * mu_log{x in [x_min, X]: E(x)>0}
                       = the brief's literal (1/ln X) prefactor; equal to
                         delta_local * ln(X/x_min)/ln(X). Reported for
                         transparency; smaller than delta_local because
                         the [2, x_min] portion is unmeasured.

    Log-measure integral is approximated by trapezoidal weights in ln(x).
    Binomial standard error is reported as a *lower bound* on uncertainty
    (it ignores E(x) autocorrelation across samples).
    """
    xs, e_actual = _load_layer1(layer1_csv_path)
    x_min, x_max = int(xs[0]), int(xs[-1])
    log_x = np.log(xs.astype(np.float64))

    # Trapezoidal weights in ln(x): w_0 = (ln x_1 - ln x_0)/2; ...
    w = np.empty_like(log_x)
    w[1:-1] = (log_x[2:] - log_x[:-2]) / 2.0
    w[0] = (log_x[1] - log_x[0]) / 2.0
    w[-1] = (log_x[-1] - log_x[-2]) / 2.0

    indicator = (e_actual > 0).astype(np.float64)

    log_measure_pos = float(np.sum(w * indicator))     # mu_log{E>0}
    log_measure_total = float(np.sum(w))               # ~ ln(X/x_min)

    delta_local = log_measure_pos / log_measure_total
    delta_rs_strict = log_measure_pos / float(np.log(x_max))

    n = len(xs)
    n_pos = int(np.sum(e_actual > 0))
    sample_fraction = n_pos / n
    binom_se = float(np.sqrt(sample_fraction * (1 - sample_fraction) / n))

    return {
        'x_min': x_min,
        'x_max': x_max,
        'n_samples': n,
        'n_positive': n_pos,
        'sample_fraction': sample_fraction,
        'log_measure_pos': log_measure_pos,
        'log_measure_total': log_measure_total,
        'delta_local': delta_local,
        'delta_rs_strict': delta_rs_strict,
        'binom_se': binom_se,
        'binom_ci95_local_lo': max(0.0, delta_local - 1.96 * binom_se),
        'binom_ci95_local_hi': min(1.0, delta_local + 1.96 * binom_se),
    }


# ---------- Path 2: Rubinstein-Sarnak Monte Carlo ----------
#
# Under GRH + LI, the limiting distribution of the normalised bias
# E(x) * log(x) / sqrt(x) under the logarithmic measure on x is
#
#     Y = c + 2 * Sum_n  cos(Theta_n) / |rho_n|
#
# with Theta_n iid uniform on [0, 2pi], |rho_n| = sqrt(1/4 + gamma_n^2),
# and c = 1 for the q=5 quadratic (squares-vs-non-squares) race.
#
# Derivation sketch: from Layer 2 we have, after amplitude-phase combine,
#   E_norm(x) = c + 2 Sum_n cos(gamma_n log x - phi_n) / |rho_n|
# where phi_n = arctan(2 gamma_n). Under LI, (gamma_n log x - phi_n) mod 2pi
# becomes equidistributed on the infinite-dim torus as x runs through the
# log-uniform measure, hence the iid uniform Theta_n. (Rubinstein-Sarnak 1994,
# Granville-Martin 2006 §3.) The c = 1 offset is the Chebyshev bias term
# already exhibited and validated in Layer 2.

def path2_rs_montecarlo(gammas, c=1.0, n_mc=1_000_000, seed=20260426,
                        chunk=None):
    """
    Monte Carlo estimate of delta = P(Y > 0) under the RS limiting
    distribution. Returns a dict of point estimate, binomial CI, and
    Y-distribution moments. Chunks the MC to bound memory.
    """
    rng = np.random.default_rng(seed)
    g = np.asarray(gammas, dtype=np.float64)
    rho_abs = np.sqrt(0.25 + g * g)
    coeffs = 2.0 / rho_abs                              # shape (N_zeros,)

    if chunk is None:
        # cap rough memory at ~100 MB float64: chunk * N * 8 bytes
        chunk = max(1, min(n_mc, 100_000_000 // (8 * len(g) + 1)))
    n_pos = 0
    Y_sum = 0.0
    Y_sq_sum = 0.0
    Y_min = float('inf')
    Y_max = float('-inf')
    done = 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        thetas = rng.uniform(0.0, 2.0 * np.pi, size=(m, len(g)))
        Y = c + np.cos(thetas) @ coeffs                  # shape (m,)
        n_pos += int(np.sum(Y > 0))
        Y_sum += float(Y.sum())
        Y_sq_sum += float((Y * Y).sum())
        if Y.size:
            Y_min = min(Y_min, float(Y.min()))
            Y_max = max(Y_max, float(Y.max()))
        done += m

    delta = n_pos / n_mc
    se = float(np.sqrt(max(delta * (1 - delta), 0.0) / n_mc))
    Y_mean = Y_sum / n_mc
    Y_var = Y_sq_sum / n_mc - Y_mean * Y_mean

    return {
        'n_zeros': len(g),
        'c': c,
        'n_mc': n_mc,
        'n_positive': n_pos,
        'delta': delta,
        'se': se,
        'ci95_lo': max(0.0, delta - 1.96 * se),
        'ci95_hi': min(1.0, delta + 1.96 * se),
        'Y_mean': Y_mean,
        'Y_std': float(np.sqrt(max(Y_var, 0.0))),
        'Y_min': Y_min,
        'Y_max': Y_max,
        'gamma_min': float(g[0]),
        'gamma_max': float(g[-1]),
    }


def path2_n_sweep(gammas_full, c=1.0, n_mc=1_000_000, seed=20260426,
                  Ns=(10, 15, 20, 25)):
    """Truncation sensitivity: same MC seed, vary N. Same RNG state per N is
    not enforced (each call reseeds), so sweep results are independent samples,
    which is what we want for honest CIs."""
    out = {}
    for N in Ns:
        if N > len(gammas_full):
            continue
        out[N] = path2_rs_montecarlo(gammas_full[:N], c=c, n_mc=n_mc,
                                     seed=seed + N)
    return out


# ---------- Path 3: LMFDB / literature lookup ----------
# (TODO: implement in next commit)

def path3_lookup():
    raise NotImplementedError('path3_lookup: pending implementation')


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description='Chebyshev-bias logarithmic density for q=5, three paths')
    parser.add_argument('--n-zeros', type=int, default=25,
                        help='Number of L(chi_5) zeros for Path 2 (default 25)')
    parser.add_argument('--n-mc', type=int, default=1_000_000,
                        help='Monte Carlo samples for Path 2 (default 1e6)')
    parser.add_argument('--seed', type=int, default=20260426,
                        help='RNG seed for Path 2 (default 20260426)')
    parser.add_argument('--dps', type=int, default=30,
                        help='mpmath precision for zero refinement (default 30)')
    parser.add_argument('--layer1-csv', type=str, default='chebyshev_bias_data.csv',
                        help='Layer 1 CSV filename (relative to this script)')
    parser.add_argument('--outdir', type=str, default=None,
                        help='Output directory (default: alongside this script)')
    args = parser.parse_args()

    here = Path(__file__).parent
    outdir = Path(args.outdir) if args.outdir else here
    outdir.mkdir(parents=True, exist_ok=True)
    csv_in = here / args.layer1_csv

    print('=' * 70)
    print('CHEBYSHEV-BIAS LOGARITHMIC DENSITY for q = 5  (three paths)')
    print('=' * 70)

    print(f'\n--- Path 1 (OBSERVED): empirical from Layer 1 sample ---')
    print(f'  source: {csv_in.name}')
    p1 = path1_empirical(csv_in)
    print(f'  range:               [{p1["x_min"]:,}, {p1["x_max"]:,}]')
    print(f'  samples:             {p1["n_samples"]} '
          f'({p1["n_positive"]} with E(x) > 0)')
    print(f'  log-measure window:  ln(X/x_min) = {p1["log_measure_total"]:.4f},  '
          f'ln X = {np.log(p1["x_max"]):.4f}')
    print(f'  sample fraction:     {p1["sample_fraction"]:.4f}')
    print(f'  delta_local:         {p1["delta_local"]:.4f}  '
          f'[95% binom CI {p1["binom_ci95_local_lo"]:.4f}, '
          f'{p1["binom_ci95_local_hi"]:.4f}]')
    print(f'  delta_rs_strict:     {p1["delta_rs_strict"]:.4f}  '
          f'(equals delta_local * ln(X/x_min)/ln X)')
    print(f'  note: binomial SE assumes independence; E(x) is autocorrelated,')
    print(f'        so true uncertainty is larger than {p1["binom_se"]:.4f}.')

    print(f'\n--- Path 2 (DERIVED, GRH+LI): Rubinstein-Sarnak Monte Carlo ---')
    print(f'  refining {args.n_zeros} zeros at dps={args.dps}...')
    t0 = time.time()
    gammas_mp = get_zeros(n=args.n_zeros, dps=args.dps, validate=True,
                          verbose=False)
    gammas = np.array([float(g) for g in gammas_mp], dtype=np.float64)
    print(f'    done in {time.time() - t0:.1f}s; '
          f'gamma_1={gammas[0]:.4f}, gamma_N={gammas[-1]:.4f}')

    print(f'  running MC (n_mc={args.n_mc:,}, seed={args.seed})...')
    t0 = time.time()
    p2 = path2_rs_montecarlo(gammas, c=1.0, n_mc=args.n_mc, seed=args.seed)
    print(f'    done in {time.time() - t0:.1f}s')
    print(f'  Y distribution: mean={p2["Y_mean"]:+.4f}, std={p2["Y_std"]:.4f}, '
          f'range [{p2["Y_min"]:+.3f}, {p2["Y_max"]:+.3f}]')
    print(f'  delta_RS = {p2["delta"]:.5f}  '
          f'[95% binom CI {p2["ci95_lo"]:.5f}, {p2["ci95_hi"]:.5f}]  '
          f'(SE = {p2["se"]:.5f})')

    print(f'\n  N-sweep (truncation sensitivity, same n_mc each):')
    sweep = path2_n_sweep(gammas, c=1.0, n_mc=args.n_mc, seed=args.seed)
    print(f'  {"N":>3} {"delta":>8} {"se":>8} {"Y_std":>8} {"Y_min":>9}')
    for N in sorted(sweep.keys()):
        s = sweep[N]
        print(f'  {N:>3d} {s["delta"]:>8.5f} {s["se"]:>8.5f} '
              f'{s["Y_std"]:>8.4f} {s["Y_min"]:>+9.3f}')

    # TODO: Path 3, write CSV, convergence summary.
    print('\n(Path 3 pending in next commit.)')


if __name__ == '__main__':
    main()
