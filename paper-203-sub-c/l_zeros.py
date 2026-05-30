"""
Compute zeros of L(s, chi_d) for d in {3, 5, 7, 13} using mpmath.

L(s, chi) is built via mpmath.dirichlet (which uses Hurwitz zeta internally
for analytic continuation). Zeros are located by sign change of the Hardy
function Z(t, chi), then refined by bisection.

For each character chi_d, all four are EVEN (chi(-1) = +1), so the
appropriate Hardy phase is:
    theta(t) = arg Gamma(1/4 + it/2) - (t/2) log(q/pi)
Z(t) = exp(i theta(t)) L(1/2 + it, chi)  is real (for real characters).

Sign changes of Z are zeros of L.

Caching: zeros saved to CSV files in this directory.
"""

from __future__ import annotations

import csv
import math
import os
import time
from pathlib import Path
import mpmath as mp

from chars import chi, CHI_TABLE

mp.mp.dps = 20

HALF = mp.mpf('0.5')
LOG_PI = mp.log(mp.pi)


def _char_coeffs(d):
    """Build [chi(0), chi(1), ..., chi(q-1)] for mpmath.dirichlet."""
    q = CHI_TABLE[d]['conductor']
    return [mp.mpf(chi(d, n)) for n in range(q)]


def L_chi(d, s, coeffs=None):
    """Compute L(s, chi_d) via mpmath.dirichlet (Hurwitz expansion)."""
    if coeffs is None:
        coeffs = _char_coeffs(d)
    # mpmath.dirichlet: sum_{n=1}^inf a_n / n^s where a_n = coeffs[n mod q]
    return mp.dirichlet(s, coeffs)


def theta_even(t, q):
    """Hardy phase for an EVEN primitive real character of conductor q.

    Lambda(s,chi) = (q/pi)^(s/2) Gamma(s/2) L(s,chi) satisfies
    Lambda(s,chi) = eps(chi) Lambda(1-s,chi). On s = 1/2 + it:
        Lambda = (q/pi)^(1/4) · |Gamma| · exp(i [arg Gamma(1/4+it/2)
                                                 + (t/2) log(q/pi)]) · L

    For Lambda real (up to sign eps), the bracketed phase cancels phase(L).
    Hardy Z(t) = exp(i·theta(t)) · L(1/2+it, chi) is real-valued.
    """
    g = mp.loggamma(mp.mpc(0.25, t / 2))
    arg_g = g.imag
    return arg_g + (t / 2) * (mp.log(q) - LOG_PI)


def Z_chi(t, d, coeffs=None, q=None):
    """Hardy Z function for L(s, chi_d). Real-valued for even real character."""
    if coeffs is None:
        coeffs = _char_coeffs(d)
    if q is None:
        q = CHI_TABLE[d]['conductor']
    s = mp.mpc(HALF, t)
    L = mp.dirichlet(s, coeffs)
    th = theta_even(t, q)
    # Z = e^{i theta} L should be real (up to numerical noise)
    Z = mp.exp(mp.mpc(0, th)) * L
    return Z.real, abs(Z.imag)  # return real and imag-magnitude (sanity)


def find_zeros(d, t_max, dt=0.5, refine=True):
    """Find zeros of L(1/2 + it, chi_d) for t in (0, t_max]."""
    coeffs = _char_coeffs(d)
    q = CHI_TABLE[d]['conductor']
    zeros = []
    t = mp.mpf('1.0')
    prev_t = t
    prev_Z, _ = Z_chi(t, d, coeffs, q)
    n_steps = 0
    while t <= t_max:
        t += dt
        Z, _ = Z_chi(t, d, coeffs, q)
        n_steps += 1
        if prev_Z * Z < 0:
            # sign change → zero in (prev_t, t)
            if refine:
                # Bisection refinement
                lo, hi = prev_t, t
                Zlo = prev_Z
                for _ in range(40):
                    mid = (lo + hi) / 2
                    Zmid, _ = Z_chi(mid, d, coeffs, q)
                    if Zlo * Zmid < 0:
                        hi = mid
                    else:
                        lo = mid
                        Zlo = Zmid
                    if hi - lo < mp.mpf('1e-8'):
                        break
                zeros.append(float((lo + hi) / 2))
            else:
                zeros.append(float((prev_t + t) / 2))
        prev_t = t
        prev_Z = Z
    return zeros


def cache_path(d, t_max):
    return Path(__file__).parent / f"zeros_chi{d}_tmax{int(t_max)}.csv"


def get_or_compute_zeros(d, t_max, dt=0.5, force=False):
    p = cache_path(d, t_max)
    if p.exists() and not force:
        zeros = []
        with open(p, 'r') as f:
            r = csv.reader(f)
            next(r)
            for row in r:
                zeros.append(float(row[1]))
        return zeros
    print(f"  computing zeros for chi_{d} up to t = {t_max} ...")
    t0 = time.time()
    zeros = find_zeros(d, t_max, dt=dt)
    print(f"    found {len(zeros)} zeros in {time.time()-t0:.1f}s")
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zeros):
            w.writerow([i + 1, f"{z:.10f}"])
    return zeros


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--d', type=int, default=5)
    parser.add_argument('--tmax', type=float, default=300.0)
    parser.add_argument('--dt', type=float, default=0.5)
    args = parser.parse_args()
    zeros = get_or_compute_zeros(args.d, args.tmax, args.dt, force=True)
    print(f"Found {len(zeros)} zeros of L(s, chi_{args.d}) in (0, {args.tmax}]")
    print(f"First few: {zeros[:5]}")
    print(f"Last few:  {zeros[-5:]}")
