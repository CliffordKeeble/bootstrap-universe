"""
Sub C-2: effect size epsilon(N) at the (5, 5) cell across N.

eps(N) = (null_mean_delta - signal_mean_delta) / null_mean_delta

The (5, 5) cell is the Paper 150 v2.0 baseline. We use the d=5 probe
against the first N Riemann zeros (this is the Paper 150 v2.0 (5,5)_Riemann
cell — the headline z = -8.14 at N=1000, z = -22.91 at N=10^4 result).

Sub C-2's question: is eps(N) flat across N, or does it decline?
"""

from __future__ import annotations

import csv
import time
import numpy as np
import mpmath as mp

from probe import compute_probe, find_minima, match_to_targets, mc_null

mp.mp.dps = 20

N_TERMS_PROBE = 5000
WINDOW = 1.0
N_NULL = 1000
SEED = 42
DT = 0.008
T_PAD = 5.0

import sys
N_VALUES = [100, 1000, 10000] if '--with10k' in sys.argv else [100, 1000]


def get_riemann_zeros(n: int, cache_path: str = "riemann_zeros_10000.csv"):
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
        existing = zs
    else:
        existing = []
    # need to fetch
    print(f"  fetching Riemann zeros {len(existing)+1}..{n}...")
    t0 = time.time()
    zs = list(existing)
    for k in range(len(existing) + 1, n + 1):
        zs.append(float(mp.zetazero(k).imag))
        if k % 500 == 0:
            print(f"    {k}/{n} (t_{k} = {zs[-1]:.2f}, {time.time()-t0:.0f}s)")
    with open(cache_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zs):
            w.writerow([i + 1, f"{z:.10f}"])
    return np.array(zs)


def run_one_N(N_zeros: int, all_zeros: np.ndarray) -> dict:
    zeros = all_zeros[:N_zeros]
    t_min = max(1.0, zeros[0] - T_PAD)
    t_max = zeros[-1] + T_PAD
    t_array = np.arange(t_min, t_max, DT)
    re, im, probe = compute_probe(t_array, d=5, N=N_TERMS_PROBE)
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    matches = match_to_targets(minima_t, zeros, window=WINDOW)
    signal_mean = float(np.mean([m['delta'] for m in matches]))
    null = mc_null(zeros, n_minima=len(minima_t), t_min=t_min, t_max=t_max,
                   window=WINDOW, n_trials=N_NULL, seed=SEED)
    z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
    eff = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0
    return dict(
        N_zeros=N_zeros, n_minima=len(minima_t),
        t_min=t_min, t_max=t_max,
        signal_mean=signal_mean,
        null_mean=null['mean'], null_std=null['std'],
        z=float(z), effect_size=float(eff),
    )


def main():
    print("=" * 70)
    print("Sub C-2: effect size epsilon(N) at (5, 5)_Riemann")
    print("=" * 70)

    print(f"\nFetching up to {max(N_VALUES)} Riemann zeros...")
    all_zeros = get_riemann_zeros(max(N_VALUES))
    print(f"  Got {len(all_zeros)} zeros, last = {all_zeros[-1]:.2f}")

    results = []
    for N in N_VALUES:
        print(f"\n--- N = {N} ---")
        t0 = time.time()
        r = run_one_N(N, all_zeros)
        print(f"  Done in {time.time()-t0:.1f}s")
        print(f"  signal_mean = {r['signal_mean']:.6f}, null_mean = {r['null_mean']:.6f}")
        print(f"  z = {r['z']:.4f}, effect_size = {r['effect_size']:.4f}")
        results.append(r)

    # Write CSV
    with open('effect_size_vs_N.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    print("\nWrote effect_size_vs_N.csv")

    # Verdict
    print("\n=== Sub C-2 verdict ===")
    eps_values = [r['effect_size'] for r in results]
    eps_max = max(eps_values)
    eps_min = min(eps_values)
    print(f"  eps range: [{eps_min:.4f}, {eps_max:.4f}]")
    if eps_max > 0 and eps_min / eps_max >= 0.8:
        print("  FLAT: effect size stable within +/- 20% across N — no horizon problem")
    elif eps_values == sorted(eps_values, reverse=True):
        print("  DECLINING: effect size monotonically decreasing with N — horizon problem real")
    else:
        print("  MIXED/STABLE-THEN-DECLINING: interpretation deferred")

    return results


if __name__ == '__main__':
    main()
