"""
Run the 4x4 (d_probe, d_target) matrix.

Compute each probe once (per d_probe), then match against each of 4 target
L-zero sets. Per-cell: signal mean delta, null MC, z-score, effect size,
PASS/FAIL/AMBIGUOUS per pre-registered thresholds.
"""

from __future__ import annotations

import csv
import math
import time
import numpy as np

from probe import compute_probe, find_minima, match_to_targets, mc_null
from l_fast import cache_path_fast as cache_path
from precompute_l_zeros import T_MAX_PER_D

DISCRIMINANTS = [3, 5, 7, 13]

N_TERMS_PROBE = 5000
DT = 0.008
WINDOW = 1.0
N_NULL = 1000
SEED = 42
T_PAD = 5.0


def load_l_zeros(d):
    p = cache_path(d, T_MAX_PER_D[d])
    if not p.exists():
        raise FileNotFoundError(f"L-zero cache missing for d={d}: {p}")
    zs = []
    with open(p, 'r') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            zs.append(float(row[1]))
    return np.array(zs)


def threshold_verdict(z):
    a = abs(z)
    if a >= 5.0:
        return 'PASS'
    if a >= 3.0:
        return 'AMBIGUOUS'
    return 'FAIL'


def run_matrix(probe_t_max_extra=20.0):
    # Load all L-zero sets
    print("Loading L-function zeros...")
    target_zeros = {}
    for d in DISCRIMINANTS:
        z = load_l_zeros(d)
        target_zeros[d] = z
        print(f"  d={d}: {len(z)} zeros, range [{z[0]:.2f}, {z[-1]:.2f}]")

    # Determine the maximum t-range we need (covers all target zero ranges)
    t_max_global = max(z[-1] for z in target_zeros.values()) + probe_t_max_extra
    t_min_global = max(0.5, min(z[0] for z in target_zeros.values()) - 2.0)
    print(f"\nGlobal probe t-range: [{t_min_global:.1f}, {t_max_global:.1f}]")

    # Compute probe for each d_probe ONCE
    print(f"\nComputing probes (N_terms={N_TERMS_PROBE}, dt={DT})...")
    t_array = np.arange(t_min_global, t_max_global, DT)
    print(f"  Grid points: {len(t_array)}")
    probes = {}
    for d in DISCRIMINANTS:
        t0 = time.time()
        re, im, probe = compute_probe(t_array, d=d, N=N_TERMS_PROBE)
        probes[d] = probe
        print(f"  d_probe={d}: probe computed in {time.time()-t0:.1f}s")

    # For each (d_probe, d_target), match and compute z
    print(f"\n=== 4x4 MATRIX ===")
    print(f"  rows = d_probe (norm |Re^2 - d*Im^2|)")
    print(f"  cols = d_target (zeros of L(chi_{{d}}))\n")
    results = []
    header = f"  {'d_probe \\ d_target':>20} | " + " | ".join(f"d'={d:2d}" for d in DISCRIMINANTS)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for d_probe in DISCRIMINANTS:
        row_strs = []
        for d_tgt in DISCRIMINANTS:
            tgts = target_zeros[d_tgt]
            # Restrict to the t-range covered by the probe; pad slightly
            t_lo = max(t_min_global, tgts[0] - T_PAD)
            t_hi = min(t_max_global, tgts[-1] + T_PAD)
            mask = (t_array >= t_lo) & (t_array <= t_hi)
            t_local = t_array[mask]
            probe_local = probes[d_probe][mask]
            minima_t, _ = find_minima(t_local, probe_local,
                                       percentile=5.0, dedup_window=0.3)
            matches = match_to_targets(minima_t, tgts, window=WINDOW)
            signal_mean = float(np.mean([m['delta'] for m in matches]))
            null = mc_null(tgts, n_minima=len(minima_t),
                           t_min=t_lo, t_max=t_hi, window=WINDOW,
                           n_trials=N_NULL, seed=SEED)
            z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
            eff = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0
            verdict = threshold_verdict(z)
            results.append({
                'd_probe': d_probe, 'd_target': d_tgt,
                'n_targets': len(tgts), 'n_minima': len(minima_t),
                't_lo': t_lo, 't_hi': t_hi,
                'signal_mean': signal_mean,
                'null_mean': null['mean'],
                'null_std': null['std'],
                'z': z, 'effect_size': eff,
                'verdict': verdict,
            })
            tag = '*' if d_probe == d_tgt else ' '
            row_strs.append(f"{tag}z={z:+6.2f}{tag}")
        print(f"  {f'd={d_probe:2d}':>20} | " + " | ".join(row_strs))

    # Write CSV
    out_csv = 'matrix_N1000.csv'
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"\nWrote {out_csv}")

    # Diagonal summary
    print(f"\n=== Diagonal verdict ===")
    diag_passes = 0
    diag_fails = 0
    for r in results:
        if r['d_probe'] == r['d_target']:
            print(f"  (d={r['d_probe']}, d={r['d_target']}): "
                  f"z = {r['z']:+.4f}, effect = {r['effect_size']:+.4f}, "
                  f"n_zeros = {r['n_targets']}, n_minima = {r['n_minima']}, "
                  f"verdict = {r['verdict']}")
            if r['verdict'] == 'PASS':
                diag_passes += 1
            elif r['verdict'] == 'FAIL':
                diag_fails += 1

    print(f"\n=== Off-diagonal summary ===")
    off_passes = 0
    off_fails = 0
    for r in results:
        if r['d_probe'] != r['d_target']:
            if r['verdict'] == 'PASS':
                off_passes += 1
            elif r['verdict'] == 'FAIL':
                off_fails += 1
    print(f"  off-diagonal PASS count: {off_passes} / 12")
    print(f"  off-diagonal FAIL count: {off_fails} / 12")
    print(f"\n  diagonal PASS count: {diag_passes} / 4")
    print(f"  diagonal FAIL count: {diag_fails} / 4")

    # Overall verdict
    if diag_passes == 4 and off_passes == 0:
        verdict = 'DIAGONAL (clean): all diagonal PASS, all off-diagonal FAIL'
    elif diag_passes == 4 and off_passes == 12:
        verdict = 'MATRIX: all cells PASS'
    elif diag_fails >= 3:
        verdict = 'DIAGONAL-FAIL: 3+ diagonal cells FAIL'
    elif diag_passes >= 2 and off_passes >= 1:
        verdict = 'MIXED'
    else:
        verdict = 'NEEDS REVIEW'
    print(f"\n  Overall: {verdict}")
    return results, verdict


if __name__ == '__main__':
    run_matrix()
