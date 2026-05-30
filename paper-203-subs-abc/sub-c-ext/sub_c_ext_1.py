"""
Sub C-ext-1: positive-definite probe control.

Probe: probe_pos(t) = sqrt(Re^2(Z_phi(t)) + 5 * Im^2(Z_phi(t)))
Target: zeros of L(s, chi_5) cached from Sub C-1.
Compare z-score against Sub C-1's (5,5) cell (z = -8.90, eps = 0.298).

Pre-registered thresholds:
  |z| < 3: CONFIRMS-INDEFINITENESS-ESSENTIAL
  |z| >= 5: REFUTES-INDEFINITENESS-CLAIM
  3 <= |z| < 5: AMBIGUOUS
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path
import numpy as np

# Add Sub C-1 directory to path for imports
SUB_C_DIR = Path(__file__).parent.parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))

from probe import compute_probe as compute_probe_indef  # noqa
from probe import find_minima, match_to_targets, mc_null  # noqa

PHI = (1.0 + math.sqrt(5.0)) / 2.0


def compute_probe_pos_def(t_array: np.ndarray, d: int = 5, N: int = 5000,
                           alpha: float = PHI):
    """Positive-definite probe: sqrt(Re^2 + d*Im^2). No null cone."""
    t = np.asarray(t_array, dtype=np.float64)
    re = np.zeros_like(t)
    im = np.zeros_like(t)
    for n in range(1, N + 1):
        theta_n = 2.0 * math.pi * ((n * alpha) % 1.0)
        sigma_n = 1.0 - n / N
        amp = sigma_n / math.sqrt(n)
        phase = t * math.log(n) + theta_n
        re += amp * np.cos(phase)
        im += amp * np.sin(phase)
    probe = np.sqrt(re * re + float(d) * im * im)
    return re, im, probe


def load_chi5_zeros(extended: bool = False) -> np.ndarray:
    if extended:
        p = SUB_C_DIR / "zeros_chi5_tmax3000_fast.csv"
    else:
        p = SUB_C_DIR / "zeros_chi5_tmax1000_fast.csv"
    zs = []
    with open(p, 'r') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            zs.append(float(row[1]))
    return np.array(zs)


def run_cell(probe_fn, label: str, target_zeros: np.ndarray,
             N_terms: int = 5000, dt: float = 0.008, t_pad: float = 5.0,
             window: float = 1.0, n_null: int = 1000, seed: int = 42):
    t_min = max(1.0, float(target_zeros[0]) - t_pad)
    t_max = float(target_zeros[-1]) + t_pad
    t_array = np.arange(t_min, t_max, dt)
    print(f"  [{label}] grid: t in [{t_min:.1f}, {t_max:.1f}], {len(t_array)} pts")
    t0 = time.time()
    re, im, probe = probe_fn(t_array)
    print(f"  [{label}] probe computed in {time.time()-t0:.1f}s")
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    density = len(minima_t) / (t_max - t_min)
    print(f"  [{label}] minima: {len(minima_t)} (density {density:.3f}/unit)")
    matches = match_to_targets(minima_t, target_zeros, window=window)
    signal_mean = float(np.mean([m['delta'] for m in matches]))
    t0 = time.time()
    null = mc_null(target_zeros, n_minima=len(minima_t),
                   t_min=t_min, t_max=t_max, window=window,
                   n_trials=n_null, seed=seed)
    print(f"  [{label}] null computed in {time.time()-t0:.1f}s")
    z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
    eff = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0
    return dict(
        label=label, n_targets=len(target_zeros), n_minima=len(minima_t),
        signal_mean=signal_mean, null_mean=null['mean'], null_std=null['std'],
        z=float(z), effect=float(eff),
    )


def main():
    print("=" * 70)
    print("Sub C-ext-1: positive-definite probe control")
    print("=" * 70)

    import sys
    extended = '--extended' in sys.argv
    zeros = load_chi5_zeros(extended=extended)
    print(f"Loaded {len(zeros)} chi_5 zeros from Sub C-1 cache "
          f"({'extended t=3000' if extended else 'baseline t=1000'})")
    print(f"Range: [{zeros[0]:.2f}, {zeros[-1]:.2f}]")

    # Indefinite baseline (matches Sub C-1 (5,5) cell)
    print("\n--- Indefinite baseline (matches Sub C-1 (5,5) cell) ---")
    r_indef = run_cell(lambda t: compute_probe_indef(t, d=5, N=5000),
                       "indefinite |Re^2 - 5*Im^2|", zeros)
    print(f"  z = {r_indef['z']:.4f}, effect = {r_indef['effect']:.4f}")

    # Positive-definite test
    print("\n--- Positive-definite test |Re^2 + 5*Im^2| ---")
    r_pos = run_cell(lambda t: compute_probe_pos_def(t, d=5, N=5000),
                     "positive-definite", zeros)
    print(f"  z = {r_pos['z']:.4f}, effect = {r_pos['effect']:.4f}")

    print("\n=== Sub C-ext-1 verdict ===")
    z_abs = abs(r_pos['z'])
    if z_abs < 3:
        verdict = "CONFIRMS-INDEFINITENESS-ESSENTIAL"
    elif z_abs >= 5:
        verdict = "REFUTES-INDEFINITENESS-CLAIM"
    else:
        verdict = "AMBIGUOUS"
    print(f"  Positive-definite z = {r_pos['z']:.4f} -> {verdict}")
    print(f"  Indefinite baseline z = {r_indef['z']:.4f} (effect {r_indef['effect']:.4f})")

    # Save results
    import csv as _csv
    with open('sub_c_ext_1_results.csv', 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(['probe', 'n_targets', 'n_minima', 'signal_mean',
                    'null_mean', 'null_std', 'z', 'effect_size'])
        for r in (r_indef, r_pos):
            w.writerow([r['label'], r['n_targets'], r['n_minima'],
                        r['signal_mean'], r['null_mean'], r['null_std'],
                        r['z'], r['effect']])
    print("Wrote sub_c_ext_1_results.csv")
    return r_indef, r_pos


if __name__ == '__main__':
    main()
