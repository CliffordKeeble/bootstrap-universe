# null_test.py
# Null hypothesis test: does ζ_φ matching beat uniform random minima?
# Paper 150 investigation step 3
# Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 12 2026

"""
Standalone: python null_test.py
Generates 1000 sets of fake minima at the same density as ζ_φ's real minima,
matches each set to Riemann zeros, and compares mean-Δ distribution to actual.
Reports z-score / empirical p-value per range bin.
"""

import numpy as np
from zeta_phi import (compute_zeta_phi, find_minima, get_riemann_zeros,
                      match_zeros, bin_matches)

# ── Parameters ────────────────────────────────────────────────────────────────

N_TERMS = 5000
T_MIN, T_MAX, DT = 1.0, 237.0, 0.008
WINDOW = 1.0
N_TRIALS = 1000
BINS = [(0, 80), (80, 160), (160, 237)]

# ── Reproduce real minima ─────────────────────────────────────────────────────

print("=" * 70)
print("NULL HYPOTHESIS TEST")
print("=" * 70)

print("\nComputing zeta_phi and finding real minima...")
t_array = np.arange(T_MIN, T_MAX, DT)
re, im, golden_norm = compute_zeta_phi(t_array, N=N_TERMS)
minima_t, minima_v = find_minima(t_array, golden_norm, percentile=5, dedup_window=0.3)
real_density = len(minima_t) / (T_MAX - T_MIN)

print(f"  Real minima: {len(minima_t)}, density: {real_density:.3f}/unit")

# Also compute density per bin
for lo, hi in BINS:
    count = np.sum((minima_t >= lo) & (minima_t < hi))
    density = count / (hi - lo)
    print(f"  [{lo:3.0f},{hi:3.0f}): {count} minima, density {density:.3f}/unit")

# ── Real matching ─────────────────────────────────────────────────────────────

print("\nFetching Riemann zeros...")
riemann_zeros = get_riemann_zeros(100)
in_range = riemann_zeros[riemann_zeros <= T_MAX]

real_matches = match_zeros(minima_t, in_range, window=WINDOW)

# Per-bin real mean delta
real_binned = {}
for lo, hi in BINS:
    deltas = [m['delta'] for m in real_matches
              if lo <= m['zero'] < hi and not np.isnan(m['delta'])]
    real_binned[(lo, hi)] = {
        'mean_delta': np.mean(deltas) if deltas else np.nan,
        'median_delta': np.median(deltas) if deltas else np.nan,
        'match_rate': np.mean([m['delta'] <= WINDOW for m in real_matches
                               if lo <= m['zero'] < hi]) if deltas else np.nan,
        'count': len(deltas)
    }

# ── Window sensitivity (Step 2) ──────────────────────────────────────────────

print("\n" + "=" * 70)
print("WINDOW SENSITIVITY ANALYSIS")
print("=" * 70)

windows = [1.0, 0.8, 0.5, 0.3, 0.2]

print(f"\n{'W':>5s}  {'Range':>12s}  {'N':>4s}  {'Match':>5s}  "
      f"{'Rate':>6s}  {'Mean_D':>7s}  {'Med_D':>7s}")
print("-" * 60)

for W in windows:
    matches_w = match_zeros(minima_t, in_range, window=W)
    for lo, hi in BINS:
        in_bin = [m for m in matches_w if lo <= m['zero'] < hi]
        deltas = [m['delta'] for m in in_bin if not np.isnan(m['delta'])]
        matched = sum(1 for m in in_bin if m['matched'])
        n = len(in_bin)
        if n == 0:
            continue
        print(f"{W:5.1f}  [{lo:3.0f},{hi:3.0f})  {n:4d}  {matched:5d}  "
              f"{matched/n:6.1%}  {np.mean(deltas):7.4f}  {np.median(deltas):7.4f}")
    print()

# ── Null hypothesis Monte Carlo ──────────────────────────────────────────────

print("=" * 70)
print(f"MONTE CARLO NULL TEST ({N_TRIALS} trials)")
print("=" * 70)

rng = np.random.default_rng(42)

# Use per-bin density for more accurate null
null_mean_deltas = {b: [] for b in BINS}
null_match_rates = {b: [] for b in BINS}

for trial in range(N_TRIALS):
    # Generate fake minima at same density as real
    n_fake = len(minima_t)
    fake_minima = rng.uniform(T_MIN, T_MAX, size=n_fake)
    fake_minima.sort()

    fake_matches = match_zeros(fake_minima, in_range, window=WINDOW)

    for lo, hi in BINS:
        deltas = [m['delta'] for m in fake_matches
                  if lo <= m['zero'] < hi and not np.isnan(m['delta'])]
        matched = sum(1 for m in fake_matches
                      if lo <= m['zero'] < hi and m['matched'])
        n = sum(1 for m in fake_matches if lo <= m['zero'] < hi)

        if deltas:
            null_mean_deltas[(lo, hi)].append(np.mean(deltas))
            null_match_rates[(lo, hi)].append(matched / n if n else 0)

    if (trial + 1) % 200 == 0:
        print(f"  Trial {trial+1}/{N_TRIALS}")

# ── Results ───────────────────────────────────────────────────────────────────

print(f"\n{'Range':>12s}  {'Real_mD':>8s}  {'Null_mD':>8s}  {'Null_sd':>8s}  "
      f"{'z-score':>8s}  {'p-value':>8s}  {'Verdict':>10s}")
print("-" * 75)

for lo, hi in BINS:
    real_md = real_binned[(lo, hi)]['mean_delta']
    null_arr = np.array(null_mean_deltas[(lo, hi)])

    if len(null_arr) == 0:
        continue

    null_mean = np.mean(null_arr)
    null_std = np.std(null_arr)

    if null_std > 0:
        z_score = (real_md - null_mean) / null_std
    else:
        z_score = 0

    # Empirical p-value: fraction of null trials with mean_delta <= real
    p_value = np.mean(null_arr <= real_md)

    if z_score < -3:
        verdict = "STRONG"
    elif z_score < -2:
        verdict = "GOOD"
    elif z_score < -1:
        verdict = "WEAK"
    else:
        verdict = "CHANCE"

    print(f"[{lo:3.0f},{hi:3.0f})  {real_md:8.4f}  {null_mean:8.4f}  "
          f"{null_std:8.4f}  {z_score:8.2f}  {p_value:8.4f}  {verdict:>10s}")

# Also report match rates
print(f"\n{'Range':>12s}  {'Real_rate':>10s}  {'Null_rate':>10s}  "
      f"{'Null_sd':>8s}  {'z-score':>8s}")
print("-" * 55)

for lo, hi in BINS:
    real_rate = real_binned[(lo, hi)]['match_rate']
    null_arr = np.array(null_match_rates[(lo, hi)])

    if len(null_arr) == 0:
        continue

    null_mean = np.mean(null_arr)
    null_std = np.std(null_arr)

    if null_std > 0:
        z = (real_rate - null_mean) / null_std
    else:
        z = 0

    print(f"[{lo:3.0f},{hi:3.0f})  {real_rate:10.1%}  {null_mean:10.1%}  "
          f"{null_std:8.4f}  {z:8.2f}")

print("\nInterpretation:")
print("  z < -3: ζ_φ matching is many sigma better than chance (STRONG signal)")
print("  z ~ -1 to -2: modest improvement over chance (WEAK signal)")
print("  z ~ 0: indistinguishable from chance (no signal)")
print("  Negative z = ζ_φ is CLOSER than random; positive z = WORSE than random")

print("\n" + "=" * 70)
print("NULL TEST COMPLETE")
print("=" * 70)
