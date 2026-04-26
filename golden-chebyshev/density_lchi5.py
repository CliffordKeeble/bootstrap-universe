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
# (TODO: implement in next commit)

def path2_rs_montecarlo(gammas, c=1.0, n_mc=1_000_000, seed=20260426):
    raise NotImplementedError('path2_rs_montecarlo: pending implementation')


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

    # TODO: Path 2, Path 3, write CSV, convergence summary.
    print('\n(Paths 2 and 3 pending in subsequent commits.)')


if __name__ == '__main__':
    main()
