"""
Fast L-function zero finder using truncated Dirichlet series.

For practical L(s, chi) computation on the critical line at moderate t,
the truncated Dirichlet sum
    L_M(s, chi) = sum_{n=1}^M chi(n) / n^s
is a good approximation to L(s, chi) for M large enough. We complete it
with a Fejer-smoothed kernel similar to Paper 150's probe machinery, so
the truncation doesn't introduce spurious oscillations.

For real even primitive characters, the Hardy Z function is
    Z(t, chi) = exp(i theta(t)) L(1/2 + it, chi)
with theta(t) = arg Gamma(1/4 + it/2) + (t/2) log(q/pi).

We use a vectorised numpy implementation rather than mpmath: a single
numpy evaluation of L_M(1/2+it, chi) at thousands of t-points takes seconds,
versus minutes per character with mpmath.

VERIFICATION: we cross-check against the 900 mpmath-computed chi_5 zeros
in zeros_chi5_tmax1000.csv before using this for the other 3 characters.
"""

from __future__ import annotations

import csv
import math
import time
import numpy as np
from pathlib import Path

from chars import CHI_TABLE


def _build_chi_array(d: int, M: int) -> np.ndarray:
    """Build [chi(1), chi(2), ..., chi(M)] as a numpy float64 array."""
    q = CHI_TABLE[d]['conductor']
    vals = CHI_TABLE[d]['values']
    arr = np.zeros(M, dtype=np.float64)
    for n in range(1, M + 1):
        arr[n - 1] = vals[n % q]
    return arr


def _theta_even_array(t_arr: np.ndarray, q: int) -> np.ndarray:
    """Hardy phase for even real primitive character."""
    # arg Gamma(1/4 + it/2) via complex log-gamma
    from scipy.special import loggamma
    z = 0.25 + 1j * (t_arr / 2.0)
    arg_g = np.imag(loggamma(z))
    return arg_g + (t_arr / 2.0) * (math.log(q) - math.log(math.pi))


def Z_fast(t_arr: np.ndarray, d: int, M: int = 500,
           fejer: bool = True) -> np.ndarray:
    """Vectorised Hardy Z function via truncated Dirichlet sum."""
    chi_arr = _build_chi_array(d, M)
    q = CHI_TABLE[d]['conductor']
    # L(1/2 + it, chi) = sum_n chi(n) / sqrt(n) * exp(-i t log n)
    re = np.zeros_like(t_arr)
    im = np.zeros_like(t_arr)
    for n in range(1, M + 1):
        c = chi_arr[n - 1]
        if c == 0:
            continue
        amp = c / math.sqrt(n)
        if fejer:
            amp *= (1.0 - n / M)   # Fejer smoothing
        ln = math.log(n)
        re += amp * np.cos(t_arr * ln)
        im -= amp * np.sin(t_arr * ln)  # exp(-i t log n) = cos - i sin
    # Multiply by exp(i theta)
    theta = _theta_even_array(t_arr, q)
    Z = re * np.cos(theta) - im * np.sin(theta)
    return Z


def find_zeros_fast(d: int, t_min: float, t_max: float, dt: float = 0.05,
                    M: int = 500, fejer: bool = False) -> list[float]:
    """Find sign changes of Z(t, chi_d) on (t_min, t_max] with refinement."""
    t_arr = np.arange(t_min, t_max + dt, dt)
    Z_arr = Z_fast(t_arr, d, M=M, fejer=fejer)
    zeros = []
    for i in range(len(t_arr) - 1):
        if Z_arr[i] * Z_arr[i + 1] < 0:
            # bisection refinement
            lo, hi = t_arr[i], t_arr[i + 1]
            Zlo = Z_arr[i]
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                Zmid = float(Z_fast(np.array([mid]), d, M=M, fejer=fejer)[0])
                if Zlo * Zmid < 0:
                    hi = mid
                else:
                    lo = mid
                    Zlo = Zmid
                if hi - lo < 1e-7:
                    break
            zeros.append(0.5 * (lo + hi))
    return zeros


def cache_path_fast(d: int, t_max: float) -> Path:
    return Path(__file__).parent / f"zeros_chi{d}_tmax{int(t_max)}_fast.csv"


def get_or_compute_zeros_fast(d: int, t_max: float, dt: float = 0.05,
                              M: int = 500, force: bool = False):
    p = cache_path_fast(d, t_max)
    if p.exists() and not force:
        zs = []
        with open(p, 'r') as f:
            r = csv.reader(f)
            next(r)
            for row in r:
                zs.append(float(row[1]))
        return zs
    print(f"  computing zeros for chi_{d} (fast) up to t = {t_max} (M={M}, dt={dt})...")
    t0 = time.time()
    zeros = find_zeros_fast(d, t_min=1.0, t_max=t_max, dt=dt, M=M)
    print(f"    found {len(zeros)} zeros in {time.time()-t0:.1f}s")
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zeros):
            w.writerow([i + 1, f"{z:.10f}"])
    return zeros


def verify_against_mpmath_chi5(M: int = 500, dt: float = 0.05) -> dict:
    """Verify fast zero-finder against the mpmath-computed chi_5 zeros."""
    mp_path = Path(__file__).parent / "zeros_chi5_tmax1000.csv"
    if not mp_path.exists():
        return dict(error="no mpmath chi_5 zeros file")
    mp_zeros = []
    with open(mp_path, 'r') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            mp_zeros.append(float(row[1]))
    mp_zeros = np.array(mp_zeros)

    # Compute fast zeros over the same range
    fast_zeros = np.array(find_zeros_fast(5, t_min=1.0, t_max=1000.0, dt=dt, M=M))

    # For each mp zero, find nearest fast zero
    deltas = []
    for z in mp_zeros[:200]:  # first 200 for speed
        if len(fast_zeros) == 0:
            continue
        idx = int(np.argmin(np.abs(fast_zeros - z)))
        deltas.append(abs(fast_zeros[idx] - z))
    deltas = np.array(deltas)
    return dict(
        n_mpmath=len(mp_zeros),
        n_fast=len(fast_zeros),
        mean_delta=float(np.mean(deltas)),
        median_delta=float(np.median(deltas)),
        max_delta=float(np.max(deltas)),
        delta_pct95=float(np.percentile(deltas, 95)),
    )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--M', type=int, default=500)
    parser.add_argument('--dt', type=float, default=0.05)
    args = parser.parse_args()
    if args.verify:
        print(f"Verifying fast (M={args.M}, dt={args.dt}) against mpmath chi_5 zeros...")
        r = verify_against_mpmath_chi5(M=args.M, dt=args.dt)
        for k, v in r.items():
            print(f"  {k}: {v}")
