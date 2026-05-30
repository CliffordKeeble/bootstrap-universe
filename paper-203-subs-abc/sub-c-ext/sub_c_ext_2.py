"""
Sub C-ext-2: cubic L-function detection test.

Probe: |Re^2 - 5*Im^2| (Paper 150 canonical, d_probe = 5).
Target: zeros of L(s, chi_2 mod 7) (cubic complex character).

Pre-registered thresholds:
  |z| < 3: NOT DETECTED (mechanism is real-quadratic-specific)
  |z| >= 5: DETECTED (mechanism is generic to L-functions)
  3 <= |z| < 5: AMBIGUOUS
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path
import numpy as np

# Add Sub C-1 directory to path
SUB_C_DIR = Path(__file__).parent.parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))

from probe import compute_probe, find_minima, match_to_targets, mc_null
from l_zeros_complex import get_or_compute as get_cubic_zeros


def load_zeros():
    return np.array(get_cubic_zeros(t_max=500.0, dt=0.05, M=1000, force=False))


def run():
    print("=" * 70)
    print("Sub C-ext-2: cubic L(chi_2 mod 7) detection test")
    print("=" * 70)

    zeros = load_zeros()
    print(f"\nLoaded {len(zeros)} cubic L-zeros, range [{zeros[0]:.2f}, {zeros[-1]:.2f}]")

    # Probe parameters (matching Sub C-1)
    N_TERMS = 5000
    DT = 0.008
    T_PAD = 5.0
    WINDOW = 1.0
    N_NULL = 1000
    SEED = 42
    D_PROBE = 5

    t_min = max(1.0, float(zeros[0]) - T_PAD)
    t_max = float(zeros[-1]) + T_PAD
    t_array = np.arange(t_min, t_max, DT)
    print(f"Probe grid: t in [{t_min:.1f}, {t_max:.1f}], {len(t_array)} pts")

    t0 = time.time()
    re, im, probe = compute_probe(t_array, d=D_PROBE, N=N_TERMS)
    print(f"Probe |Re^2 - 5*Im^2| computed in {time.time()-t0:.1f}s")
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    print(f"Minima: {len(minima_t)} (density {len(minima_t)/(t_max-t_min):.3f}/unit)")

    matches = match_to_targets(minima_t, zeros, window=WINDOW)
    signal_mean = float(np.mean([m['delta'] for m in matches]))
    print(f"Signal mean delta: {signal_mean:.6f}")

    t0 = time.time()
    null = mc_null(zeros, n_minima=len(minima_t),
                   t_min=t_min, t_max=t_max, window=WINDOW,
                   n_trials=N_NULL, seed=SEED)
    print(f"Null computed in {time.time()-t0:.1f}s")
    print(f"Null mean: {null['mean']:.6f}, std: {null['std']:.6f}")

    z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
    eff = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0

    print(f"\n  z-score: {z:.4f}")
    print(f"  effect size: {eff:.4f}")

    z_abs = abs(z)
    if z_abs < 3:
        verdict = "NOT DETECTED (mechanism is real-quadratic-specific)"
    elif z_abs >= 5:
        verdict = "DETECTED (mechanism is generic to L-functions)"
    else:
        verdict = "AMBIGUOUS"
    print(f"\nVerdict: {verdict}")

    # Comparison reference: Sub C-1 (5, 7) cell with real L(chi_28) zeros
    # gave z = -5.50, effect = 0.282. The cubic L(chi_2 mod 7) test uses
    # the same probe but a different target.
    print(f"\nReference comparison (Sub C-1 (5,7) cell, real chi_28):")
    print(f"  z = -5.50, effect = 0.282 — DETECTED")
    print(f"\nThis test (Sub C-ext-2, cubic chi_2 mod 7):")
    print(f"  z = {z:.4f}, effect = {eff:.4f} — {verdict.split()[0]}")

    with open('sub_c_ext_2_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['target', 'n_targets', 'n_minima', 'signal_mean',
                    'null_mean', 'null_std', 'z', 'effect_size', 'verdict'])
        w.writerow(['L(chi_2 mod 7) cubic', len(zeros), len(minima_t),
                    signal_mean, null['mean'], null['std'], z, eff, verdict])
    print("\nWrote sub_c_ext_2_results.csv")
    return dict(z=z, effect=eff, verdict=verdict)


if __name__ == '__main__':
    run()
