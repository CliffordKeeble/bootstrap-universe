"""
L-function zero finder for complex characters.

For chi a complex Dirichlet character mod q, L(s, chi) is complex-valued
on the critical line. Zeros of L are local minima of |L|^2.

We use the truncated Dirichlet sum (M=1000) matching the validated
fast method from Sub C-1. For complex chi, the Hardy Z function is
not naturally real, so we just find minima of |L|.

Validation: cross-check L(s, chi_2 mod 7) first zero against mpmath
via Hurwitz zeta evaluation.
"""

from __future__ import annotations

import csv
import math
import time
from pathlib import Path
import numpy as np

from chi2_mod7 import CHI_2_MOD_7, chi2


def _build_chi2_array(M: int) -> np.ndarray:
    """Build [chi(1), chi(2), ..., chi(M)] for chi_2 mod 7."""
    arr = np.zeros(M, dtype=np.complex128)
    for n in range(1, M + 1):
        arr[n - 1] = CHI_2_MOD_7[n % 7]
    return arr


def L_complex_truncated(t_arr: np.ndarray, M: int = 1000) -> np.ndarray:
    """L(1/2 + it, chi_2 mod 7) via truncated Dirichlet sum, vectorised."""
    chi_arr = _build_chi2_array(M)
    L = np.zeros_like(t_arr, dtype=np.complex128)
    for n in range(1, M + 1):
        c = chi_arr[n - 1]
        if c == 0:
            continue
        amp = c / math.sqrt(n)
        ln = math.log(n)
        L += amp * np.exp(-1j * t_arr * ln)
    return L


def find_zeros_complex(t_min: float, t_max: float, dt: float = 0.05,
                       M: int = 1000) -> list[float]:
    """Find local minima of |L|^2; refine with bisection-on-Re*Im sign."""
    t_arr = np.arange(t_min, t_max + dt, dt)
    L = L_complex_truncated(t_arr, M=M)
    abs2 = (L.real ** 2 + L.imag ** 2)

    # Find local minima (three-point) where |L|^2 is also small
    # We want sign changes of both Re(L) and Im(L); a "true zero" has
    # both close to 0.
    zeros = []
    # Use a percentile-of-|L|^2 threshold to filter to small values
    thr = np.percentile(abs2, 5.0)
    for i in range(1, len(t_arr) - 1):
        if abs2[i] < abs2[i - 1] and abs2[i] < abs2[i + 1] and abs2[i] < thr:
            # Refine: find the local minimum more precisely via parabolic interp
            x0, x1, x2 = t_arr[i - 1], t_arr[i], t_arr[i + 1]
            y0, y1, y2 = abs2[i - 1], abs2[i], abs2[i + 1]
            # Parabolic vertex
            denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
            if abs(denom) > 1e-15:
                A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
                B = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / denom
                if A > 1e-15:
                    t_min_local = -B / (2 * A)
                    if x0 < t_min_local < x2:
                        zeros.append(float(t_min_local))
                        continue
            zeros.append(float(x1))
    # Dedupe: minima within 0.1 of each other → keep one
    deduped = []
    for z in zeros:
        if not deduped or z - deduped[-1] >= 0.1:
            deduped.append(z)
    return deduped


def cache_path(t_max: float) -> Path:
    return Path(__file__).parent / f"zeros_chi2_mod7_tmax{int(t_max)}.csv"


def get_or_compute(t_max: float = 500.0, dt: float = 0.05, M: int = 1000,
                   force: bool = False) -> list[float]:
    p = cache_path(t_max)
    if p.exists() and not force:
        zs = []
        with open(p, 'r') as f:
            r = csv.reader(f)
            next(r)
            for row in r:
                zs.append(float(row[1]))
        return zs
    print(f"  computing L(s, chi_2 mod 7) zeros up to t={t_max} (M={M}, dt={dt})...")
    t0 = time.time()
    zeros = find_zeros_complex(1.0, t_max, dt=dt, M=M)
    print(f"    found {len(zeros)} zero candidates in {time.time()-t0:.1f}s")
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zeros):
            w.writerow([i + 1, f"{z:.10f}"])
    return zeros


def validate_against_mpmath(zeros: list[float], n_check: int = 5) -> dict:
    """Cross-check first few zeros against mpmath Hurwitz-based L evaluation."""
    try:
        import mpmath as mp
    except ImportError:
        return dict(error="mpmath not available")
    mp.mp.dps = 20
    chi_vals = [mp.mpc(CHI_2_MOD_7[a % 7].real, CHI_2_MOD_7[a % 7].imag) for a in range(7)]
    q = 7

    def L_mp(t):
        """L(1/2 + it, chi_2) via mpmath Hurwitz."""
        s = mp.mpc('0.5', t)
        # L(s, chi) = q^(-s) * sum_{a=1}^{q-1} chi(a) * zeta(s, a/q)
        result = mp.mpc(0)
        for a in range(1, q):
            c = chi_vals[a]
            if c != 0:
                result += c * mp.zeta(s, mp.mpf(a) / q)
        return result / mp.mpf(q) ** s

    out = []
    for i, t in enumerate(zeros[:n_check]):
        L_val = L_mp(t)
        abs_L = abs(L_val)
        out.append((i + 1, float(t), float(abs_L)))
    return dict(checks=out)


if __name__ == '__main__':
    print("=" * 60)
    print("L(s, chi_2 mod 7) cubic character zero finder")
    print("=" * 60)
    zs = get_or_compute(t_max=500.0, dt=0.05, M=1000, force=True)
    print(f"\n{len(zs)} zeros found")
    print(f"First 10: {[f'{z:.3f}' for z in zs[:10]]}")
    print(f"Last 3:   {[f'{z:.3f}' for z in zs[-3:]]}")

    print("\nValidation: |L| at first 5 candidate zeros via mpmath...")
    v = validate_against_mpmath(zs, n_check=5)
    if 'checks' in v:
        for n, t, abs_L in v['checks']:
            print(f"  zero {n}: t = {t:.4f}, |L(1/2 + it, chi_2)| = {abs_L:.6f}")
