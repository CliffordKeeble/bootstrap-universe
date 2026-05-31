"""
Sub F follow-up — RvM-random draw distribution across 10 seeds.

Reuses Sub F machinery; only the RvM seed varies.

Per PRE_REGISTRATION.md, committed before any new seed is run.
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path
import numpy as np

# Reuse Sub F + Sub C-1 machinery
SUB_F_DIR = Path(__file__).parent.parent / "paper-150-sub-f"
SUB_C_DIR = Path(__file__).parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))
sys.path.insert(0, str(SUB_F_DIR))
from probe import compute_probe, find_minima, match_to_targets, mc_null
from sub_f import sample_rvm_random


# Pre-registered constants (frozen — identical to Sub F)
N_TARGET = 2792
N_TERMS = 5000
WINDOW = 1.0
N_NULL = 1000
NULL_SEED = 42
DT = 0.008
T_PAD = 5.0
D_PROBE = 5

# Riemann control z from Sub F
Z_RIEMANN = -12.5629
EFF_RIEMANN = 0.2335

# 10 RvM-random seeds (frozen at pre-registration)
RVM_SEEDS = [20260530, 20260531, 20260601, 20260602, 20260603,
             20260604, 20260605, 20260606, 20260607, 20260608]


def load_riemann_first_n(n: int = N_TARGET) -> np.ndarray:
    p = SUB_C_DIR / "riemann_zeros_10000.csv"
    zs = []
    with open(p, 'r') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            zs.append(float(row[1]))
    return np.array(zs[:n])


def main():
    print("=" * 70)
    print("Sub F follow-up: RvM-random distribution across 10 seeds")
    print(f"Riemann control (Sub F): z = {Z_RIEMANN:.4f}, eff = {EFF_RIEMANN:.4f}")
    print("=" * 70)

    # T_max from Riemann cache
    riemann = load_riemann_first_n(N_TARGET)
    T_max = float(riemann[-1])
    print(f"\nRiemann cache: {len(riemann)} zeros, T_max = {T_max:.2f}")

    # Compute probe Re/Im once (target-independent)
    # Use a t-range that covers all RvM draws (which extend up to T_max + T_PAD)
    t_min = 1.0
    t_max_grid = T_max + 10.0
    t_array = np.arange(t_min, t_max_grid, DT)
    print(f"\nProbe grid: t in [{t_min}, {t_max_grid:.1f}], {len(t_array)} pts")
    t0 = time.time()
    re, im, _ = compute_probe(t_array, d=D_PROBE, N=N_TERMS)
    print(f"Probe Re/Im in {time.time()-t0:.1f}s (one-shot, reused across seeds)")
    probe = np.sqrt(np.abs(re * re - float(D_PROBE) * im * im))
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    print(f"Minima: {len(minima_t)} (density {len(minima_t)/(t_max_grid-t_min):.3f}/u)")

    # For each RvM seed
    rows = []
    print(f"\nRunning {len(RVM_SEEDS)} RvM-random seeds...")
    for seed in RVM_SEEDS:
        t0 = time.time()
        rvm = sample_rvm_random(N_TARGET, T_max + T_PAD, seed=seed)
        # restrict to probe range
        rvm_in = rvm[(rvm >= t_min) & (rvm <= t_max_grid)]
        matches = match_to_targets(minima_t, rvm_in, window=WINDOW)
        sig = float(np.mean([m['delta'] for m in matches]))
        null = mc_null(rvm_in, n_minima=len(minima_t),
                       t_min=t_min, t_max=t_max_grid, window=WINDOW,
                       n_trials=N_NULL, seed=NULL_SEED)
        z = (sig - null['mean']) / null['std'] if null['std'] > 0 else 0.0
        eff = (null['mean'] - sig) / null['mean'] if null['mean'] > 0 else 0.0
        n_matched = sum(1 for m in matches if m['matched'])
        row = dict(seed=seed, n_targets=len(rvm_in),
                   n_matched=n_matched,
                   signal_mean=sig, null_mean=null['mean'], null_std=null['std'],
                   z=z, effect=eff)
        rows.append(row)
        print(f"  seed={seed}: z = {z:+.4f}, eff = {eff:.4f}, "
              f"n_matched = {n_matched} ({time.time()-t0:.1f}s)")

    # Distribution stats
    zs = np.array([r['z'] for r in rows])
    effs = np.array([r['effect'] for r in rows])
    abs_zs = np.abs(zs)

    print(f"\n=== z distribution (10 seeds) ===")
    print(f"  mean = {zs.mean():+.4f}")
    print(f"  std  = {zs.std():+.4f}")
    print(f"  min  = {zs.min():+.4f}")
    print(f"  max  = {zs.max():+.4f}")
    print(f"  median = {np.median(zs):+.4f}")

    print(f"\n=== effect-size distribution (10 seeds) ===")
    print(f"  mean = {effs.mean():.4f}")
    print(f"  std  = {effs.std():.4f}")
    print(f"  min  = {effs.min():.4f}")
    print(f"  max  = {effs.max():.4f}")
    print(f"  median = {np.median(effs):.4f}")

    # K-sigma comparison
    K = (abs(Z_RIEMANN) - abs_zs.mean()) / abs_zs.std() if abs_zs.std() > 0 else float('inf')
    pct = float(np.mean(abs_zs >= abs(Z_RIEMANN))) * 100
    print(f"\n=== Comparison to Riemann control ===")
    print(f"  Riemann |z| = {abs(Z_RIEMANN):.4f}")
    print(f"  RvM   |z| mean = {abs_zs.mean():.4f}, std = {abs_zs.std():.4f}")
    print(f"  K-sigma separation: K = ({abs(Z_RIEMANN):.2f} - {abs_zs.mean():.2f}) / {abs_zs.std():.3f} = {K:.3f}")
    print(f"  Riemann |z| >= max(RvM |z|)? {abs(Z_RIEMANN) >= abs_zs.max()}")
    print(f"  Fraction of RvM draws with |z| >= Riemann: {pct:.1f}%")

    if K < 1:
        verdict = "CONFIRM SUSTAIN (K < 1)"
    elif K < 3:
        verdict = "PARTIAL DISTINCTION (1 <= K < 3)"
    else:
        verdict = "REJECT SUSTAIN (K >= 3)"
    print(f"\nVerdict: {verdict}")

    # Write CSV
    with open('sub_f_followup_rvm_distribution.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['seed', 'n_targets', 'n_matched', 'signal_mean',
                    'null_mean', 'null_std', 'z', 'effect_size'])
        for r in rows:
            w.writerow([r['seed'], r['n_targets'], r['n_matched'],
                        r['signal_mean'], r['null_mean'], r['null_std'],
                        r['z'], r['effect']])
    print("Wrote sub_f_followup_rvm_distribution.csv")
    return dict(rows=rows, K=K, verdict=verdict,
                abs_z_mean=abs_zs.mean(), abs_z_std=abs_zs.std())


if __name__ == '__main__':
    main()
