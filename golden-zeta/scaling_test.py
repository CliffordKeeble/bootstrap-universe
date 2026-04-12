# scaling_test.py
# Scaling in N: does matching improve at N=20000?
# Paper 150 investigation step 4
# Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 12 2026

"""
Standalone: python scaling_test.py
Recomputes ζ_φ at N=5000, 10000, 20000 and compares match quality.
"""

import numpy as np
from zeta_phi import (compute_zeta_phi, find_minima, get_riemann_zeros,
                      match_zeros, bin_matches)

BINS = [(0, 80), (80, 160), (160, 237)]
T_MIN, T_MAX, DT = 1.0, 237.0, 0.008
WINDOW = 1.0

print("=" * 70)
print("SCALING TEST: N = 5000, 10000, 20000")
print("=" * 70)

t_array = np.arange(T_MIN, T_MAX, DT)

print("\nFetching Riemann zeros...")
riemann_zeros = get_riemann_zeros(100)
in_range = riemann_zeros[riemann_zeros <= T_MAX]

for N in [5000, 10000, 20000]:
    print(f"\n{'─' * 50}")
    print(f"N = {N}")
    print(f"{'─' * 50}")

    print("  Computing zeta_phi...")
    re, im, golden_norm = compute_zeta_phi(t_array, N=N)

    minima_t, minima_v = find_minima(t_array, golden_norm,
                                      percentile=5, dedup_window=0.3)
    print(f"  Minima found: {len(minima_t)}, "
          f"density: {len(minima_t)/(T_MAX-T_MIN):.3f}/unit")

    matches = match_zeros(minima_t, in_range, window=WINDOW)
    all_deltas = [m['delta'] for m in matches if not np.isnan(m['delta'])]

    print(f"  Overall: mean_delta={np.mean(all_deltas):.4f}, "
          f"median={np.median(all_deltas):.4f}")

    binned = bin_matches(matches, bins=BINS)

    print(f"  {'Range':>12s}  {'N':>4s}  {'Rate':>6s}  "
          f"{'Mean_D':>7s}  {'Med_D':>7s}  {'Max_D':>7s}")
    for (lo, hi), stats in sorted(binned.items()):
        print(f"  [{lo:3.0f},{hi:3.0f})  {stats['count']:4d}  "
              f"{stats['match_rate']:6.1%}  "
              f"{stats['mean_delta']:7.4f}  {stats['median_delta']:7.4f}  "
              f"{stats['max_delta']:7.4f}")

print(f"\n{'─' * 50}")
print("INTERPRETATION")
print(f"{'─' * 50}")
print("""
If mean_delta at t>80 shrinks substantially with N:
  => finite-N artefact, Paper 150's conjecture supported.

If mean_delta at t>80 stays ~0.3-0.5 regardless of N:
  => not a finite-N artefact, the gap is structural.
""")

print("=" * 70)
print("SCALING TEST COMPLETE")
print("=" * 70)
