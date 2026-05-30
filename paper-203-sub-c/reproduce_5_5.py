"""
Reproduction of Paper 150 v2.0 (5,5) baseline at N=1000 Riemann zeros.

Per brief and pre-registration: if z is not approximately -8.14 with |z| >= 5,
the pipeline has a methodological problem and we STOP.

This uses the EXACT Paper 150 v2.0 setup:
- N_TERMS = 5000 (Dirichlet sum truncation)
- alpha = phi (golden angle)
- norm: |Re^2 - 5*Im^2| (d = 5)
- target: Riemann zeta zeros (not L(chi_5))
- W = 1.0
- N_null = 1000, seed = 42
- N_zeros = 1000
"""

import csv
import time
import numpy as np
import mpmath as mp

from probe import compute_probe, find_minima, match_to_targets, mc_null

mp.mp.dps = 20

N_TERMS = 5000
N_ZEROS = 1000
WINDOW = 1.0
N_NULL = 1000
SEED = 42
DT = 0.008
T_PAD = 5.0


def get_riemann_zeros(n: int, cache_path: str = "riemann_zeros_1000.csv"):
    import os
    if os.path.exists(cache_path):
        zs = []
        with open(cache_path, 'r') as f:
            r = csv.reader(f)
            next(r, None)
            for row in r:
                zs.append(float(row[1]))
        if len(zs) >= n:
            return np.array(zs[:n])
    print(f"  fetching {n} Riemann zeros...")
    t0 = time.time()
    zs = []
    for k in range(1, n + 1):
        zs.append(float(mp.zetazero(k).imag))
        if k % 100 == 0:
            print(f"    {k}/{n} (t_{k} = {zs[-1]:.2f}, {time.time()-t0:.0f}s)")
    with open(cache_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zs):
            w.writerow([i + 1, f"{z:.10f}"])
    return np.array(zs)


def main():
    print("=" * 70)
    print("PIPELINE-SOUNDNESS REPRODUCTION: (5, 5)_Riemann at N=1000 zeros")
    print(f"Target from Paper 150 v2.0: z = -8.14 (run at SHA c2b9d0d)")
    print("=" * 70)

    print(f"\nFetching first {N_ZEROS} Riemann zeros...")
    zeros = get_riemann_zeros(N_ZEROS)
    print(f"  Got {len(zeros)} zeros, range [{zeros[0]:.2f}, {zeros[-1]:.2f}]")

    t_min = max(1.0, zeros[0] - T_PAD)
    t_max = zeros[-1] + T_PAD
    print(f"\nProbe grid: t in [{t_min:.1f}, {t_max:.1f}], dt={DT}")
    t_array = np.arange(t_min, t_max, DT)
    print(f"  Grid points: {len(t_array)}")

    print(f"\nComputing probe |Re^2 - 5*Im^2|, N_terms={N_TERMS}...")
    t0 = time.time()
    re, im, probe = compute_probe(t_array, d=5, N=N_TERMS)
    print(f"  Done in {time.time()-t0:.1f}s")

    minima_t, minima_v = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    density = len(minima_t) / (t_max - t_min)
    print(f"\n  Probe minima: {len(minima_t)}, density {density:.3f}/unit")

    matches = match_to_targets(minima_t, zeros, window=WINDOW)
    signal_mean = float(np.mean([m['delta'] for m in matches]))
    print(f"  Signal mean delta: {signal_mean:.6f}")

    print(f"\nMC null ({N_NULL} trials, seed={SEED})...")
    t0 = time.time()
    null = mc_null(zeros, n_minima=len(minima_t),
                   t_min=t_min, t_max=t_max, window=WINDOW,
                   n_trials=N_NULL, seed=SEED)
    print(f"  Done in {time.time()-t0:.1f}s")
    print(f"  Null mean delta: {null['mean']:.6f}, std: {null['std']:.6f}")

    z = (signal_mean - null['mean']) / null['std']
    eff = (null['mean'] - signal_mean) / null['mean']

    print(f"\n  z-overall: {z:.4f}")
    print(f"  effect size: {eff:.4f}")
    print(f"  Paper 150 v2.0 target: z = -8.14")
    print(f"  Ratio (z / -8.14):     {z / -8.14:.3f}")

    if abs(z) >= 5.0:
        print(f"\n  PIPELINE-SOUND: |z| >= 5; reproduction successful.")
    else:
        print(f"\n  *** PIPELINE FAIL *** |z| < 5; STOP. ***")

    return z, signal_mean, null['mean'], null['std']


if __name__ == '__main__':
    main()
