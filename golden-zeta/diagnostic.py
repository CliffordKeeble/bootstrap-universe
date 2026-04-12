# diagnostic.py
# Per-zero forensic diagnostic table
# Paper 150 investigation step 5
# Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 12 2026

"""
Standalone: python diagnostic.py
Produces per_zero.csv with forensic detail for all 100 zeros.
"""

import csv
import numpy as np
from zeta_phi import (compute_zeta_phi, find_minima, get_riemann_zeros,
                      match_zeros)

N_TERMS = 5000
T_MIN, T_MAX, DT = 1.0, 237.0, 0.008
WINDOW = 1.0
N_NULL_TRIALS = 200
OUTPUT_FILE = 'golden-zeta/per_zero.csv'

print("=" * 70)
print("PER-ZERO DIAGNOSTIC")
print("=" * 70)

# ── Compute real minima ───────────────────────────────────────────────────────

t_array = np.arange(T_MIN, T_MAX, DT)
re, im, golden_norm = compute_zeta_phi(t_array, N=N_TERMS)
minima_t, minima_v = find_minima(t_array, golden_norm, percentile=5, dedup_window=0.3)

# ── Get Riemann zeros ─────────────────────────────────────────────────────────

riemann_zeros = get_riemann_zeros(100)
in_range = riemann_zeros[riemann_zeros <= T_MAX]
n_zeros = len(in_range)

# ── Match ─────────────────────────────────────────────────────────────────────

matches = match_zeros(minima_t, in_range, window=WINDOW)

# ── Riemann zero spacings ────────────────────────────────────────────────────

spacings = np.diff(in_range)
# Pad to same length (last zero gets spacing from previous)
spacings = np.append(spacings, spacings[-1] if len(spacings) > 0 else np.nan)

# ── Null: mean nearest-fake-minimum delta per zero ────────────────────────────

print(f"Computing null reference ({N_NULL_TRIALS} trials)...")
rng = np.random.default_rng(42)
null_deltas = np.zeros(n_zeros)

for trial in range(N_NULL_TRIALS):
    fake = np.sort(rng.uniform(T_MIN, T_MAX, size=len(minima_t)))
    for i, z in enumerate(in_range):
        dists = np.abs(fake - z)
        null_deltas[i] += np.min(dists)

null_deltas /= N_NULL_TRIALS

# ── Build CSV ─────────────────────────────────────────────────────────────────

print(f"Writing {OUTPUT_FILE}...")

rows = []
for i, m in enumerate(matches):
    row = {
        'n': m['n'],
        'riemann_zero_t': f"{m['zero']:.6f}",
        'matched_minimum_t': f"{m['min_t']:.6f}" if not np.isnan(m['min_t']) else '',
        'delta': f"{m['delta']:.6f}" if not np.isnan(m['delta']) else '',
        'local_spacing': f"{spacings[i]:.6f}",
        'delta_over_spacing': f"{m['delta']/spacings[i]:.4f}" if not np.isnan(m['delta']) and spacings[i] > 0 else '',
        'null_mean_delta': f"{null_deltas[i]:.6f}",
        'delta_vs_null': f"{m['delta']/null_deltas[i]:.4f}" if not np.isnan(m['delta']) and null_deltas[i] > 0 else '',
        'matched_W1': 'Y' if m['matched'] else 'N',
    }
    rows.append(row)

fieldnames = ['n', 'riemann_zero_t', 'matched_minimum_t', 'delta',
              'local_spacing', 'delta_over_spacing', 'null_mean_delta',
              'delta_vs_null', 'matched_W1']

with open(OUTPUT_FILE, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"\n  {len(rows)} rows written to {OUTPUT_FILE}")

# ── Summary stats ─────────────────────────────────────────────────────────────

deltas = np.array([float(r['delta']) for r in rows if r['delta']])
null_d = np.array([float(r['null_mean_delta']) for r in rows])
ratios = deltas / null_d[:len(deltas)]

print(f"\n  delta/null ratio: mean={ratios.mean():.3f}, median={np.median(ratios):.3f}")
print(f"  (ratio < 1 means ζ_φ is closer than random)")

# Per bin
for lo, hi in [(0, 80), (80, 160), (160, 237)]:
    mask = (in_range[:len(deltas)] >= lo) & (in_range[:len(deltas)] < hi)
    if mask.any():
        r = ratios[mask]
        print(f"  [{lo:3.0f},{hi:3.0f}): mean ratio={r.mean():.3f}, "
              f"median={np.median(r):.3f}, "
              f"<1: {np.sum(r<1)}/{len(r)}")

print("\n" + "=" * 70)
print("DIAGNOSTIC COMPLETE")
print("=" * 70)
