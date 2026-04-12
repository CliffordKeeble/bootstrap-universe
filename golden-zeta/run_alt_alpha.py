# run_alt_alpha.py — Test 2: Alternative equidistributed sequences
# Pre-registration SHA: c2b9d0d
# Bootstrap Universe Programme
# Written by Mr Code, April 12 2026

"""
Standalone: python run_alt_alpha.py
Test 2: repeat N=1000 primary with alpha = sqrt(2), e, pi.
Tests phi-specificity of the z=-8.14 signal.
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

ALPHAS = {
    'phi':   PHI,
    'sqrt2': math.sqrt(2),
    'e':     math.e,
    'pi':    math.pi,
}

print("=" * 70)
print("TEST 2: ALTERNATIVE EQUIDISTRIBUTED SEQUENCES")
print(f"Pre-registration SHA: c2b9d0d")
print("=" * 70)

# ── Fetch 1000 Riemann zeros (reuse from cache if possible) ──────────────────

print(f"\nFetching 1000 Riemann zeros...")
t0 = time.time()
riemann_zeros = []
for n in range(1, 1001):
    z = float(mpmath.zetazero(n).imag)
    riemann_zeros.append(z)
    if n % 200 == 0:
        print(f"  {n}/1000 (t_{n} = {z:.2f})")

riemann_zeros = np.array(riemann_zeros)
T_MAX = riemann_zeros[-1] + 5
t_array = np.arange(T_MIN, T_MAX, DT)
print(f"  Range: [{riemann_zeros[0]:.2f}, {riemann_zeros[-1]:.2f}]")
print(f"  Fetched in {time.time()-t0:.0f}s")

# ── Run all four alphas ──────────────────────────────────────────────────────

results = {}
for name, alpha in ALPHAS.items():
    print(f"\n{'─' * 50}")
    print(f"Running alpha = {name} ({alpha:.8f})")
    print(f"{'─' * 50}")

    t0 = time.time()
    r = run_test(riemann_zeros, t_array, N_terms=N_TERMS, alpha=alpha,
                 norm_type='golden', n_mc=N_MC, rng_seed=RNG_SEED, label=name)
    print(f"  Completed in {time.time()-t0:.1f}s")
    print_result(r)
    results[name] = r

    # Write per-zero CSV
    if name != 'phi':
        csv_file = f'golden-zeta/per_zero_{name}.csv'
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['n', 't_n', 'nearest_minimum_t', 'delta_n', 'delta_n_null'])
            for i, m in enumerate(r['matches']):
                writer.writerow([m['n'], f"{m['zero']:.6f}", f"{m['min_t']:.6f}",
                                f"{m['delta']:.6f}", f"{r['null_per_zero'][i]:.6f}"])
        print(f"  Written to {csv_file}")

# ── Summary ──────────────────────────────────────────────────────────────────

print(f"\n{'=' * 70}")
print("TEST 2 SUMMARY")
print(f"{'=' * 70}")

print(f"\n  {'Alpha':>8s}  {'z_overall':>10s}  {'Sig_mD':>8s}  {'Null_mD':>8s}")
print("  " + "-" * 40)

for name in ['phi', 'sqrt2', 'e', 'pi']:
    r = results[name]
    print(f"  {name:>8s}  {r['z_overall']:10.4f}  "
          f"{r['signal_mean_delta']:8.4f}  {r['null_mean_delta']:8.4f}")

# Decision per protocol
z_phi = abs(results['phi']['z_overall'])
z_others = [abs(results[k]['z_overall']) for k in ['sqrt2', 'e', 'pi']]

all_strong = all(z > 5 for z in z_others)
all_weak = all(z < 3 for z in z_others)

print()
if all_strong:
    print("  RESULT: All alphas show |z| > 5.")
    print("  Effect is GENERIC to equidistribution.")
    print("  phi-specificity claim UNSUPPORTED.")
elif all_weak:
    print("  RESULT: All alternative alphas show |z| < 3.")
    print("  phi is SPECIFIC. Paper's core claim STRENGTHENED.")
else:
    print("  RESULT: Mixed — some alphas strong, some weak.")
    for name in ['sqrt2', 'e', 'pi']:
        z = results[name]['z_overall']
        print(f"    {name}: z = {z:.2f} ({'strong' if abs(z) > 5 else 'suggestive' if abs(z) > 3 else 'weak'})")
    print("  Partial phi-specificity. Requires separate investigation.")

print(f"\n{'=' * 70}")
print("TEST 2 COMPLETE")
print(f"{'=' * 70}")
