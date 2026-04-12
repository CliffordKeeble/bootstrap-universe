# run_norm_decoupling.py — Test 3: Norm decoupling
# Pre-registration SHA: c2b9d0d
# Bootstrap Universe Programme
# Written by Mr Code, April 12 2026

"""
Standalone: python run_norm_decoupling.py
Test 3a: {n*phi} phase + circular norm (Re^2 + Im^2)
Test 3b: {n*sqrt(2)} phase + golden norm (|Re^2 - 5*Im^2|)
Tests whether the effect requires both golden components together.
"""

import csv
import math
import time
import numpy as np
import mpmath

from zeta_phi import find_minima, match_zeros
from zeta_gen import compute_zeta_gen, run_test, print_result, PHI

N_TERMS = 5000
T_MIN, DT = 1.0, 0.008
N_MC = 1000
RNG_SEED = 42

print("=" * 70)
print("TEST 3: NORM DECOUPLING")
print(f"Pre-registration SHA: c2b9d0d")
print("=" * 70)

# ── Fetch 1000 Riemann zeros ────────────────────────────────────────────────

print(f"\nFetching 1000 Riemann zeros...")
t0 = time.time()
riemann_zeros = []
for n in range(1, 1001):
    z = float(mpmath.zetazero(n).imag)
    riemann_zeros.append(z)
    if n % 200 == 0:
        print(f"  {n}/1000")
riemann_zeros = np.array(riemann_zeros)
T_MAX = riemann_zeros[-1] + 5
t_array = np.arange(T_MIN, T_MAX, DT)
print(f"  Fetched in {time.time()-t0:.0f}s")

# ── Test 3a: golden phase + circular norm ────────────────────────────────────

print(f"\n{'─' * 50}")
print("Test 3a: phi phase + CIRCULAR norm")
print(f"{'─' * 50}")

t0 = time.time()
r3a = run_test(riemann_zeros, t_array, N_terms=N_TERMS, alpha=PHI,
               norm_type='circular', n_mc=N_MC, rng_seed=RNG_SEED,
               label="3a: phi+circular")
print(f"  Completed in {time.time()-t0:.1f}s")
print_result(r3a)

csv_file = 'golden-zeta/per_zero_circular_norm.csv'
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 't_n', 'nearest_minimum_t', 'delta_n', 'delta_n_null'])
    for i, m in enumerate(r3a['matches']):
        writer.writerow([m['n'], f"{m['zero']:.6f}", f"{m['min_t']:.6f}",
                        f"{m['delta']:.6f}", f"{r3a['null_per_zero'][i]:.6f}"])
print(f"  Written to {csv_file}")

# ── Test 3b: sqrt(2) phase + golden norm ─────────────────────────────────────

print(f"\n{'─' * 50}")
print("Test 3b: sqrt(2) phase + GOLDEN norm")
print(f"{'─' * 50}")

t0 = time.time()
r3b = run_test(riemann_zeros, t_array, N_terms=N_TERMS, alpha=math.sqrt(2),
               norm_type='golden', n_mc=N_MC, rng_seed=RNG_SEED,
               label="3b: sqrt2+golden")
print(f"  Completed in {time.time()-t0:.1f}s")
print_result(r3b)

csv_file = 'golden-zeta/per_zero_sqrt2_goldenNorm.csv'
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 't_n', 'nearest_minimum_t', 'delta_n', 'delta_n_null'])
    for i, m in enumerate(r3b['matches']):
        writer.writerow([m['n'], f"{m['zero']:.6f}", f"{m['min_t']:.6f}",
                        f"{m['delta']:.6f}", f"{r3b['null_per_zero'][i]:.6f}"])
print(f"  Written to {csv_file}")

# ── Also run the baseline (phi + golden) for direct comparison ───────────────

print(f"\n{'─' * 50}")
print("Baseline: phi phase + GOLDEN norm")
print(f"{'─' * 50}")

t0 = time.time()
r_base = run_test(riemann_zeros, t_array, N_terms=N_TERMS, alpha=PHI,
                  norm_type='golden', n_mc=N_MC, rng_seed=RNG_SEED,
                  label="baseline: phi+golden")
print(f"  Completed in {time.time()-t0:.1f}s")
print_result(r_base)

# ── Summary ──────────────────────────────────────────────────────────────────

print(f"\n{'=' * 70}")
print("TEST 3 SUMMARY")
print(f"{'=' * 70}")

print(f"\n  {'Test':>25s}  {'Phase':>8s}  {'Norm':>10s}  {'z':>8s}")
print("  " + "-" * 55)
print(f"  {'baseline':>25s}  {'phi':>8s}  {'golden':>10s}  {r_base['z_overall']:8.4f}")
print(f"  {'3a':>25s}  {'phi':>8s}  {'circular':>10s}  {r3a['z_overall']:8.4f}")
print(f"  {'3b':>25s}  {'sqrt2':>8s}  {'golden':>10s}  {r3b['z_overall']:8.4f}")

# Decision per protocol
z_base = abs(r_base['z_overall'])
z_3a = abs(r3a['z_overall'])
z_3b = abs(r3b['z_overall'])

print()
if z_3a < 3 and z_3b < 3:
    print("  RESULT: Both 3a and 3b weak (|z| < 3).")
    print("  Effect requires BOTH golden phase AND golden norm together.")
    print("  Consistent with Dedekind bridge hypothesis.")
elif z_3a > 5 and z_3b < 3:
    print("  RESULT: 3a strong, 3b weak.")
    print("  Golden PHASE alone is sufficient; norm is cosmetic.")
elif z_3a < 3 and z_3b > 5:
    print("  RESULT: 3a weak, 3b strong.")
    print("  Golden NORM alone is sufficient; any equidistributed phase works.")
    print("  Major reframing needed.")
elif z_3a > 5 and z_3b > 5:
    print("  RESULT: Both 3a and 3b strong.")
    print("  Both components carry independent signal. Unexpected.")
else:
    print(f"  RESULT: Mixed (3a z={r3a['z_overall']:.2f}, 3b z={r3b['z_overall']:.2f})")

print(f"\n{'=' * 70}")
print("TEST 3 COMPLETE")
print(f"{'=' * 70}")
