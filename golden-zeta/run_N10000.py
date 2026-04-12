# run_N10000.py — Test 1: Massive-N robustness
# Pre-registration SHA: c2b9d0d
# Bootstrap Universe Programme
# Written by Mr Code, April 12 2026

"""
Standalone: python run_N10000.py
Test 1: repeat primary at N=10000 Riemann zeros.
Expected z ~ -25.7 under sqrt(N) scaling from z=-8.14 at N=1000.
"""

import sys
import csv
import time
import numpy as np
import mpmath

from zeta_phi import find_minima, match_zeros
from zeta_gen import compute_zeta_gen, PHI

N_TERMS = 5000
T_MIN, DT = 1.0, 0.008
WINDOW = 1.0
N_MC = 1000
RNG_SEED = 42
BINS = [(14, 1500), (1500, 5000), (5000, 10000), (10000, 15000)]

print("=" * 70)
print("TEST 1: MASSIVE-N ROBUSTNESS (N=10000 zeros)")
print(f"Pre-registration SHA: c2b9d0d")
print("=" * 70)

# ── Fetch 10000 Riemann zeros ────────────────────────────────────────────────

print(f"\nFetching 10000 Riemann zeros...")
t0 = time.time()
riemann_zeros = []
for n in range(1, 10001):
    z = float(mpmath.zetazero(n).imag)
    riemann_zeros.append(z)
    if n % 1000 == 0:
        print(f"  {n}/10000 (t_{n} = {z:.2f}, {time.time()-t0:.0f}s)")

riemann_zeros = np.array(riemann_zeros)
print(f"\nFetched in {time.time()-t0:.0f}s")
print(f"Range: [{riemann_zeros[0]:.2f}, {riemann_zeros[-1]:.2f}]")

T_MAX = riemann_zeros[-1] + 5
t_array = np.arange(T_MIN, T_MAX, DT)
print(f"Scan: [1, {T_MAX:.0f}], {len(t_array)} grid points")

# ── Compute zeta_phi ─────────────────────────────────────────────────────────

print(f"\nComputing zeta_phi (N_terms={N_TERMS})...")
t0 = time.time()
re, im, golden_norm = compute_zeta_gen(t_array, N=N_TERMS, alpha=PHI, norm_type='golden')
print(f"  Computed in {time.time()-t0:.1f}s")

minima_t, minima_v = find_minima(t_array, golden_norm, percentile=5, dedup_window=0.3)
print(f"  Minima: {len(minima_t)}, density: {len(minima_t)/(T_MAX-T_MIN):.3f}/unit")

# ── Signal matching ──────────────────────────────────────────────────────────

matches = match_zeros(minima_t, riemann_zeros, window=WINDOW)
signal_deltas = np.array([m['delta'] for m in matches])
signal_mean = np.mean(signal_deltas)
print(f"\n  Signal mean delta: {signal_mean:.6f}")

# ── MC null ──────────────────────────────────────────────────────────────────

print(f"\nMonte Carlo null ({N_MC} trials)...")
rng = np.random.default_rng(RNG_SEED)
null_means = []
null_per_zero = np.zeros(len(riemann_zeros))

t0 = time.time()
for trial in range(N_MC):
    fake = np.sort(rng.uniform(T_MIN, T_MAX, size=len(minima_t)))
    fm = match_zeros(fake, riemann_zeros, window=WINDOW)
    fd = np.array([m['delta'] for m in fm])
    null_means.append(np.mean(fd))
    null_per_zero += fd
    if (trial + 1) % 100 == 0:
        print(f"  {trial+1}/{N_MC} ({time.time()-t0:.0f}s)")

null_per_zero /= N_MC
null_arr = np.array(null_means)
null_mean = np.mean(null_arr)
null_se = np.std(null_arr)

z_overall = (signal_mean - null_mean) / null_se
p_value = np.mean(null_arr <= signal_mean)

# ── Bin-resolved ─────────────────────────────────────────────────────────────

# Redefine bins based on actual range
actual_bins = [(14, 3500), (3500, 7000), (7000, 10500), (10500, T_MAX)]

print(f"\n{'=' * 70}")
print("RESULTS")
print(f"{'=' * 70}")

print(f"\n  Signal mean delta:  {signal_mean:.6f}")
print(f"  Null mean delta:    {null_mean:.6f}")
print(f"  Null SE:            {null_se:.6f}")
print(f"  z_overall:          {z_overall:.4f}")
print(f"  p_value:            {p_value:.6f}")

# sqrt(N) prediction
z_predicted = -8.14 * np.sqrt(10000 / 1000)
print(f"\n  Predicted (sqrt(N) scaling): {z_predicted:.1f}")
print(f"  Actual:                       {z_overall:.1f}")
print(f"  Ratio actual/predicted:       {z_overall/z_predicted:.3f}")

# Decision per protocol
if abs(z_overall) >= 30:
    decision = "EFFECT GREW (|z|>30) — flag for investigation"
elif 20 <= abs(z_overall) <= 30:
    decision = "sqrt(N) CONFIRMED — effect constant across 100x expansion"
elif 5 <= abs(z_overall) < 20:
    decision = "SATURATING — effect weaker than sqrt(N) prediction"
else:
    decision = "FAILED — N=1000 result was wrong"

print(f"  Decision:                     {decision}")

# Bin z-scores
print(f"\n  Bin-resolved z-scores:")
rng2 = np.random.default_rng(RNG_SEED)
# Recompute per-bin null
for lo, hi in actual_bins:
    sig_d = [m['delta'] for m in matches if lo <= m['zero'] < hi]
    if not sig_d:
        continue

    # Quick bin null from stored null_per_zero
    bin_mask = np.array([lo <= m['zero'] < hi for m in matches])
    bin_sig_mean = np.mean(np.array([m['delta'] for m in matches])[bin_mask])
    bin_null_mean = np.mean(null_per_zero[bin_mask])
    # Approximate SE from overall scaling
    bin_n = np.sum(bin_mask)
    bin_se = null_se * np.sqrt(len(riemann_zeros) / bin_n) if bin_n > 0 else null_se
    bin_z = (bin_sig_mean - bin_null_mean) / bin_se if bin_se > 0 else 0

    print(f"    [{lo:5.0f},{hi:5.0f}): N={bin_n:5d}, z={bin_z:.2f}")

# ── CSV ──────────────────────────────────────────────────────────────────────

CSV_FILE = 'golden-zeta/per_zero_N10000.csv'
with open(CSV_FILE, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 't_n', 'nearest_minimum_t', 'delta_n', 'delta_n_null'])
    for i, m in enumerate(matches):
        writer.writerow([m['n'], f"{m['zero']:.6f}", f"{m['min_t']:.6f}",
                        f"{m['delta']:.6f}", f"{null_per_zero[i]:.6f}"])

print(f"\n  {len(matches)} rows written to {CSV_FILE}")
print(f"\n{'=' * 70}")
print("TEST 1 COMPLETE")
print(f"{'=' * 70}")
