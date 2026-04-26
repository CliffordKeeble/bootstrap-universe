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
# (TODO: implement in next commit)

def path1_empirical(layer1_csv_path):
    raise NotImplementedError('path1_empirical: pending implementation')


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

    # TODO: Path 1, Path 2, Path 3, write CSV, print convergence summary.
    raise NotImplementedError('density_lchi5: scaffold only; paths pending')


if __name__ == '__main__':
    main()
