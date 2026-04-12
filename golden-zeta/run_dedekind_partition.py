# run_dedekind_partition.py — Test 4: Dedekind factorisation partition
# Pre-registration SHA: c2b9d0d
# Bootstrap Universe Programme
# Written by Mr Code, April 12 2026

"""
Standalone: python run_dedekind_partition.py
Test 4: classify zeta_phi minima against both:
  (a) Riemann zeta zeros
  (b) L(chi_5) zeros
Tests whether the Dedekind bridge is the mechanism.
"""

import csv
import math
import time
import numpy as np
import mpmath
from mpmath import mp, mpf, exp, pi, log, sqrt, fabs, cos, sin

from zeta_phi import find_minima, match_zeros
from zeta_gen import compute_zeta_gen, PHI

N_TERMS = 5000
T_MIN, DT = 1.0, 0.008
N_MC = 1000
RNG_SEED = 42

print("=" * 70)
print("TEST 4: DEDEKIND FACTORISATION PARTITION")
print(f"Pre-registration SHA: c2b9d0d")
print("=" * 70)

# ── Compute L(chi_5) zeros ───────────────────────────────────────────────────
#
# chi_5 is the Kronecker symbol (n/5):
#   chi_5(1)=1, chi_5(2)=-1, chi_5(3)=-1, chi_5(4)=1, chi_5(0)=0
# Period 5.
#
# L(s, chi_5) = sum_{n=1}^{inf} chi_5(n) / n^s
#
# Find zeros on the critical line s = 1/2 + it by computing
# |L(1/2+it, chi_5)| and finding its minima.

def chi5(n):
    """Kronecker symbol (n/5) = Legendre symbol for p=5."""
    r = n % 5
    if r == 0: return 0
    if r == 1 or r == 4: return 1
    return -1  # r == 2 or r == 3

def compute_L_chi5_vectorised(t_array, N_sum=10000):
    """
    Compute |L(1/2 + it, chi_5)| for array of t values.
    Vectorised: outer loop over n, inner numpy over all t values.
    """
    t = np.asarray(t_array, dtype=np.float64)
    re = np.zeros_like(t)
    im = np.zeros_like(t)

    for n in range(1, N_sum + 1):
        c = chi5(n)
        if c == 0:
            continue
        sigma = 1.0 - n / N_sum  # Fejer
        amp = c * sigma / math.sqrt(n)
        phase = t * math.log(n)
        re += amp * np.cos(phase)
        im += amp * np.sin(phase)

    return np.sqrt(re**2 + im**2)


print("\nComputing L(chi_5) on critical line to find zeros...")
print("(Vectorised Fejer-smoothed Dirichlet series, N=10000 terms)")

# Scan a range comparable to our Riemann zeros ([14, 1420])
t_scan = np.arange(1.0, 1425.0, 0.01)

t0 = time.time()
L_values = compute_L_chi5_vectorised(t_scan, N_sum=10000)
print(f"  L(chi_5) computed: {len(t_scan)} points in {time.time()-t0:.0f}s")

# Find zeros as minima
L_minima_t, L_minima_v = find_minima(t_scan, L_values, percentile=2, dedup_window=0.1)

# Filter to very small values (actual zeros should be near 0)
zero_threshold = np.percentile(L_minima_v, 50)  # keep deepest half
L_zeros = L_minima_t[L_minima_v < zero_threshold]

print(f"\n  L(chi_5) candidate zeros found: {len(L_zeros)}")
print(f"  Range: [{L_zeros[0]:.2f}, {L_zeros[-1]:.2f}]" if len(L_zeros) > 0 else "  No zeros found!")

# Sanity check: first few L(chi_5) zeros should be near known values
# LMFDB: first zero of L(s, chi_5) with chi_5 = Kronecker(5,.) is around t ~ 6.64
print(f"  First 10 candidate zeros: {L_zeros[:10].round(2)}")

t_Lcompute = time.time() - t0
print(f"  Computed in {t_Lcompute:.0f}s")

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
print(f"  Fetched in {time.time()-t0:.0f}s")

# ── Compute zeta_phi ────────────────────────────────────────────────────────

T_MAX = riemann_zeros[-1] + 5
t_array = np.arange(T_MIN, T_MAX, DT)

print(f"\nComputing zeta_phi...")
re, im, golden_norm = compute_zeta_gen(t_array, N=N_TERMS, alpha=PHI, norm_type='golden')
minima_t, minima_v = find_minima(t_array, golden_norm, percentile=5, dedup_window=0.3)
print(f"  Minima: {len(minima_t)}")

# ── Test 4a: zeta_phi minima vs Riemann zeros ───────────────────────────────

print(f"\n{'─' * 50}")
print("Test 4a: zeta_phi minima vs RIEMANN zeros")
print(f"{'─' * 50}")

matches_zeta = match_zeros(minima_t, riemann_zeros, window=1.0)
sig_deltas_zeta = np.array([m['delta'] for m in matches_zeta])
sig_mean_zeta = np.mean(sig_deltas_zeta)

rng = np.random.default_rng(RNG_SEED)
null_means_zeta = []
for trial in range(N_MC):
    fake = np.sort(rng.uniform(T_MIN, T_MAX, size=len(minima_t)))
    fm = match_zeros(fake, riemann_zeros, window=1.0)
    null_means_zeta.append(np.mean([m['delta'] for m in fm]))

null_arr_zeta = np.array(null_means_zeta)
z_zeta = (sig_mean_zeta - np.mean(null_arr_zeta)) / np.std(null_arr_zeta)

print(f"  Signal mean delta: {sig_mean_zeta:.6f}")
print(f"  Null mean delta:   {np.mean(null_arr_zeta):.6f}")
print(f"  z_overall (zeta):  {z_zeta:.4f}")

# ── Test 4b: zeta_phi minima vs L(chi_5) zeros ──────────────────────────────

print(f"\n{'─' * 50}")
print("Test 4b: zeta_phi minima vs L(chi_5) zeros")
print(f"{'─' * 50}")

# Use L_zeros in comparable range
L_zeros_in_range = L_zeros[L_zeros <= T_MAX]
print(f"  L(chi_5) zeros in range: {len(L_zeros_in_range)}")

if len(L_zeros_in_range) > 0:
    matches_L = match_zeros(minima_t, L_zeros_in_range, window=1.0)
    sig_deltas_L = np.array([m['delta'] for m in matches_L])
    sig_mean_L = np.mean(sig_deltas_L)

    rng2 = np.random.default_rng(RNG_SEED)
    null_means_L = []
    null_per_zero_L = np.zeros(len(L_zeros_in_range))
    for trial in range(N_MC):
        fake = np.sort(rng2.uniform(T_MIN, T_MAX, size=len(minima_t)))
        fm = match_zeros(fake, L_zeros_in_range, window=1.0)
        fd = np.array([m['delta'] for m in fm])
        null_means_L.append(np.mean(fd))
        null_per_zero_L += fd

    null_per_zero_L /= N_MC
    null_arr_L = np.array(null_means_L)
    z_L = (sig_mean_L - np.mean(null_arr_L)) / np.std(null_arr_L)

    print(f"  Signal mean delta: {sig_mean_L:.6f}")
    print(f"  Null mean delta:   {np.mean(null_arr_L):.6f}")
    print(f"  z_overall (L):     {z_L:.4f}")

    # Write CSV
    csv_file = 'golden-zeta/per_zero_Lchi5.csv'
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['n', 't_n', 'nearest_minimum_t', 'delta_n', 'delta_n_null'])
        for i, m in enumerate(matches_L):
            writer.writerow([m['n'], f"{m['zero']:.6f}", f"{m['min_t']:.6f}",
                            f"{m['delta']:.6f}", f"{null_per_zero_L[i]:.6f}"])
    print(f"  Written to {csv_file}")
else:
    z_L = 0.0
    print("  No L(chi_5) zeros found in range — cannot test.")

# ── Summary ──────────────────────────────────────────────────────────────────

print(f"\n{'=' * 70}")
print("TEST 4 SUMMARY")
print(f"{'=' * 70}")

print(f"\n  {'Target':>20s}  {'N_zeros':>8s}  {'z_overall':>10s}")
print("  " + "-" * 42)
print(f"  {'Riemann zeros':>20s}  {len(riemann_zeros):8d}  {z_zeta:10.4f}")
print(f"  {'L(chi_5) zeros':>20s}  {len(L_zeros_in_range):8d}  {z_L:10.4f}")

print()
z_z = abs(z_zeta)
z_l = abs(z_L)

if z_z > 5 and z_l > 5:
    print("  RESULT: Both factors tightened.")
    print("  zeta_phi detects BOTH zeta and L(chi_5) zeros.")
    print("  Dedekind bridge IS the mechanism. DERIVED + OBSERVED.")
elif z_z > 5 and z_l < 3:
    print("  RESULT: Only Riemann zeros tightened.")
    print("  zeta_phi detects zeta zeros specifically, not Dedekind zeta.")
    print("  Dedekind bridge hypothesis WEAKENED.")
elif z_z < 3 and z_l > 5:
    print("  RESULT: Only L(chi_5) zeros tightened.")
    print("  zeta_phi detects Legendre correction, not Riemann.")
    print("  Paper 150 v1.3 framing WRONG in unexpected direction.")
else:
    print(f"  RESULT: Mixed/ambiguous (zeta z={z_zeta:.2f}, L z={z_L:.2f})")

print(f"\n{'=' * 70}")
print("TEST 4 COMPLETE")
print(f"{'=' * 70}")
