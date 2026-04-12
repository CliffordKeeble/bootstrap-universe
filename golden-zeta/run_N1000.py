# run_N1000.py
# Pre-registered N=1000 golden zeta extension
# Executes all statistics from pre_registration_N1000.md
# Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 12 2026
#
# Pre-registration SHA: 21f74ed

"""
Standalone: python run_N1000.py
Produces:
  - All pre-registered statistics to stdout
  - golden-zeta/per_zero_1000.csv
  - Version info for reproducibility
"""

import sys
import csv
import time
import numpy as np
import mpmath
import platform

# Record versions for reproducibility
print("=" * 70)
print("ENVIRONMENT")
print("=" * 70)
print(f"Python:   {sys.version}")
print(f"numpy:    {np.__version__}")
print(f"mpmath:   {mpmath.__version__}")
print(f"Platform: {platform.platform()}")
print(f"Pre-registration SHA: 21f74ed")

# ── Import from zeta_phi.py ──────────────────────────────────────────────────

from zeta_phi import compute_zeta_phi, find_minima, match_zeros

# ── Parameters (per protocol) ────────────────────────────────────────────────

N_TERMS_PRIMARY = 5000
N_TERMS_ROBUST = 10000
T_MIN, T_MAX, DT = 1.0, 1420.0, 0.008
WINDOW_PRIMARY = 1.0
WINDOWS_SECONDARY = [0.5, 0.3, 0.2]
N_MC_TRIALS = 1000
RNG_SEED = 42
BINS = [(14, 300), (300, 650), (650, 1000), (1000, 1420)]

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 1: Fetch 1000 Riemann zeros
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PHASE 1: RIEMANN ZEROS")
print("=" * 70)

t0 = time.time()
print(f"\nFetching 1000 Riemann zeros from mpmath.zetazero...")
print(f"(This takes a while at high n — patience)")

riemann_zeros = []
for n in range(1, 1001):
    z = float(mpmath.zetazero(n).imag)
    riemann_zeros.append(z)
    if n % 100 == 0:
        print(f"  {n}/1000 done (t_{n} = {z:.4f})")

riemann_zeros = np.array(riemann_zeros)
t_fetch = time.time() - t0
print(f"\nFetched in {t_fetch:.1f}s")
print(f"Range: [{riemann_zeros[0]:.4f}, {riemann_zeros[-1]:.4f}]")

# Adjust T_MAX to cover all zeros
actual_t_max = max(T_MAX, riemann_zeros[-1] + 5)
t_array = np.arange(T_MIN, actual_t_max, DT)
print(f"Scan range: [1, {actual_t_max:.1f}], {len(t_array)} grid points")

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 2: Compute ζ_φ and find minima (N=5000)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PHASE 2: ZETA_PHI (N=5000)")
print("=" * 70)

t0 = time.time()
print(f"\nComputing zeta_phi over {len(t_array)} points, N={N_TERMS_PRIMARY}...")
re, im, golden_norm = compute_zeta_phi(t_array, N=N_TERMS_PRIMARY)
t_compute = time.time() - t0
print(f"  Computed in {t_compute:.1f}s")

minima_t, minima_v = find_minima(t_array, golden_norm, percentile=5, dedup_window=0.3)
minima_density = len(minima_t) / (actual_t_max - T_MIN)
print(f"  Minima found: {len(minima_t)}, density: {minima_density:.3f}/unit")

# Per-bin density
for lo, hi in BINS:
    count = np.sum((minima_t >= lo) & (minima_t < hi))
    density = count / (hi - lo)
    print(f"  [{lo:4.0f},{hi:4.0f}): {count} minima, density {density:.3f}/unit")

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 3: MATCH SIGNAL (W=1.0)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PHASE 3: SIGNAL MATCHING (W=1.0)")
print("=" * 70)

matches = match_zeros(minima_t, riemann_zeros, window=WINDOW_PRIMARY)
signal_deltas = np.array([m['delta'] for m in matches])
signal_mean_delta = np.mean(signal_deltas)
signal_matched = np.sum([m['matched'] for m in matches])

print(f"\n  Matched: {signal_matched}/{len(matches)} at W={WINDOW_PRIMARY}")
print(f"  Mean delta: {signal_mean_delta:.6f}")
print(f"  Median delta: {np.median(signal_deltas):.6f}")
print(f"  Max delta: {np.max(signal_deltas):.6f}")

# ══════════════════════════════════════════════════════════════════════════════
# PHASE 4: MONTE CARLO NULL (1000 trials)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print(f"PHASE 4: MONTE CARLO NULL ({N_MC_TRIALS} trials)")
print("=" * 70)

rng = np.random.default_rng(RNG_SEED)

# Storage: per-trial overall mean delta, and per-bin mean deltas
null_mean_deltas_overall = []
null_mean_deltas_binned = {b: [] for b in BINS}
null_match_rates_W03 = []  # for secondary stat 3

# Also track per-zero null deltas for CSV
null_deltas_per_zero = np.zeros(len(riemann_zeros))

t0 = time.time()
for trial in range(N_MC_TRIALS):
    n_fake = len(minima_t)
    fake_minima = np.sort(rng.uniform(T_MIN, actual_t_max, size=n_fake))

    fake_matches = match_zeros(fake_minima, riemann_zeros, window=WINDOW_PRIMARY)
    fake_deltas = np.array([m['delta'] for m in fake_matches])

    null_mean_deltas_overall.append(np.mean(fake_deltas))

    # Per-zero accumulation for CSV
    null_deltas_per_zero += fake_deltas

    # Per-bin
    for lo, hi in BINS:
        bin_deltas = [m['delta'] for m in fake_matches if lo <= m['zero'] < hi]
        if bin_deltas:
            null_mean_deltas_binned[(lo, hi)].append(np.mean(bin_deltas))

    # Match rate at W=0.3 for secondary stat 3
    fake_matches_03 = match_zeros(fake_minima, riemann_zeros, window=0.3)
    rate_03 = np.mean([m['matched'] for m in fake_matches_03])
    null_match_rates_W03.append(rate_03)

    if (trial + 1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{N_MC_TRIALS} ({elapsed:.0f}s)")

null_deltas_per_zero /= N_MC_TRIALS

null_overall = np.array(null_mean_deltas_overall)
null_mean = np.mean(null_overall)
null_se = np.std(null_overall)

# ══════════════════════════════════════════════════════════════════════════════
# RESULTS — PRIMARY STATISTIC
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PRIMARY STATISTIC")
print("=" * 70)

z_overall = (signal_mean_delta - null_mean) / null_se

print(f"\n  Signal mean delta:  {signal_mean_delta:.6f}")
print(f"  Null mean delta:    {null_mean:.6f}")
print(f"  Null SE:            {null_se:.6f}")
print(f"  z_overall:          {z_overall:.4f}")

# Empirical p-value
p_empirical = np.mean(null_overall <= signal_mean_delta)
print(f"  Empirical p-value:  {p_empirical:.4f}")

# Decision
if abs(z_overall) < 2:
    decision = "NOT SIGNIFICANT (|z| < 2)"
elif abs(z_overall) < 3:
    decision = "SUGGESTIVE (2 <= |z| < 3)"
elif abs(z_overall) < 5:
    decision = "EVIDENCE (3 <= |z| < 5)"
else:
    decision = "DISCOVERY (|z| >= 5)"

print(f"  Decision:           {decision}")

if z_overall > 0:
    print(f"  WARNING: positive z — signal is WORSE than null!")

# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY 1: BIN-RESOLVED Z-SCORES
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SECONDARY 1: BIN-RESOLVED Z-SCORES")
print("=" * 70)

print(f"\n  {'Bin':>16s}  {'N_zeros':>7s}  {'Sig_mD':>8s}  {'Null_mD':>8s}  "
      f"{'Null_SE':>8s}  {'z':>8s}")
print("  " + "-" * 65)

bin_z_scores = {}
for lo, hi in BINS:
    bin_deltas_sig = [m['delta'] for m in matches if lo <= m['zero'] < hi]
    bin_null = np.array(null_mean_deltas_binned[(lo, hi)])

    if len(bin_deltas_sig) == 0 or len(bin_null) == 0:
        continue

    sig_md = np.mean(bin_deltas_sig)
    null_md = np.mean(bin_null)
    null_sd = np.std(bin_null)

    z = (sig_md - null_md) / null_sd if null_sd > 0 else 0
    bin_z_scores[(lo, hi)] = z

    print(f"  [{lo:4.0f},{hi:4.0f})  {len(bin_deltas_sig):7d}  {sig_md:8.4f}  "
          f"{null_md:8.4f}  {null_sd:8.4f}  {z:8.4f}")

# Check uniformity
signs = [1 if z < 0 else -1 for z in bin_z_scores.values()]
uniform = all(s == signs[0] for s in signs)
print(f"\n  All bins same sign? {uniform}")
if not uniform:
    print("  WARNING: non-uniform signal. Discovery framing unavailable.")

# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY 2: WINDOW SENSITIVITY
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SECONDARY 2: WINDOW SENSITIVITY")
print("=" * 70)

print(f"\n  {'W':>5s}  {'Sig_mD':>8s}  {'Match%':>7s}  {'z_approx':>9s}")
print("  " + "-" * 35)

# At W=1.0 we already have z_overall
print(f"  {1.0:5.1f}  {signal_mean_delta:8.4f}  "
      f"{signal_matched/len(matches):7.1%}  {z_overall:9.4f}")

for W in WINDOWS_SECONDARY:
    matches_w = match_zeros(minima_t, riemann_zeros, window=W)
    sig_deltas_w = np.array([m['delta'] for m in matches_w])
    sig_matched_w = np.sum([m['matched'] for m in matches_w])

    # Quick null at this window (reuse same fake sets via seeded rng)
    rng2 = np.random.default_rng(RNG_SEED)
    null_mds_w = []
    for trial in range(N_MC_TRIALS):
        fake = np.sort(rng2.uniform(T_MIN, actual_t_max, size=len(minima_t)))
        fm = match_zeros(fake, riemann_zeros, window=W)
        null_mds_w.append(np.mean([m['delta'] for m in fm]))

    null_arr_w = np.array(null_mds_w)
    z_w = (np.mean(sig_deltas_w) - np.mean(null_arr_w)) / np.std(null_arr_w) if np.std(null_arr_w) > 0 else 0

    print(f"  {W:5.1f}  {np.mean(sig_deltas_w):8.4f}  "
          f"{sig_matched_w/len(matches_w):7.1%}  {z_w:9.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY 3: MATCH RATE AT W=0.3
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SECONDARY 3: MATCH RATE AT W=0.3")
print("=" * 70)

matches_03 = match_zeros(minima_t, riemann_zeros, window=0.3)
sig_rate_03 = np.mean([m['matched'] for m in matches_03])
null_rate_03_mean = np.mean(null_match_rates_W03)
null_rate_03_sd = np.std(null_match_rates_W03)

print(f"\n  Signal match rate (W=0.3):  {sig_rate_03:.1%}")
print(f"  Null match rate (W=0.3):   {null_rate_03_mean:.1%} +/- {null_rate_03_sd:.1%}")
z_rate = (sig_rate_03 - null_rate_03_mean) / null_rate_03_sd if null_rate_03_sd > 0 else 0
print(f"  z-score (match rate):      {z_rate:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# SECONDARY 4: N_TERMS SENSITIVITY (N=10000)
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SECONDARY 4: N_TERMS = 10000")
print("=" * 70)

print(f"\nComputing zeta_phi at N={N_TERMS_ROBUST}...")
t0 = time.time()
re2, im2, gn2 = compute_zeta_phi(t_array, N=N_TERMS_ROBUST)
print(f"  Computed in {time.time()-t0:.1f}s")

minima_t2, minima_v2 = find_minima(t_array, gn2, percentile=5, dedup_window=0.3)
print(f"  Minima found: {len(minima_t2)}, density: {len(minima_t2)/(actual_t_max-T_MIN):.3f}/unit")

matches2 = match_zeros(minima_t2, riemann_zeros, window=WINDOW_PRIMARY)
sig_deltas2 = np.array([m['delta'] for m in matches2])
sig_md2 = np.mean(sig_deltas2)

# Use same null (density-matched)
rng3 = np.random.default_rng(RNG_SEED)
null_mds2 = []
for trial in range(N_MC_TRIALS):
    fake = np.sort(rng3.uniform(T_MIN, actual_t_max, size=len(minima_t2)))
    fm = match_zeros(fake, riemann_zeros, window=WINDOW_PRIMARY)
    null_mds2.append(np.mean([m['delta'] for m in fm]))

null_arr2 = np.array(null_mds2)
z_N10k = (sig_md2 - np.mean(null_arr2)) / np.std(null_arr2) if np.std(null_arr2) > 0 else 0

print(f"\n  Signal mean delta (N=10000): {sig_md2:.6f}")
print(f"  Null mean delta:             {np.mean(null_arr2):.6f}")
print(f"  z_overall (N=10000):         {z_N10k:.4f}")
print(f"  z_overall (N=5000):          {z_overall:.4f}")
print(f"  Material shift? {'YES' if abs(z_N10k - z_overall) > 1 else 'NO'}")

# ══════════════════════════════════════════════════════════════════════════════
# PER-ZERO CSV
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("WRITING PER-ZERO CSV")
print("=" * 70)

CSV_FILE = 'golden-zeta/per_zero_1000.csv'

# Determine bin for each zero
def get_bin(t):
    for lo, hi in BINS:
        if lo <= t < hi:
            return f"[{lo},{hi})"
    return "out"

rows = []
for i, m in enumerate(matches):
    rows.append({
        'n': m['n'],
        't_n': f"{m['zero']:.6f}",
        'nearest_minimum_t': f"{m['min_t']:.6f}",
        'delta_n': f"{m['delta']:.6f}",
        'nearest_null_minimum_t': '',  # not stored per-trial
        'delta_n_null': f"{null_deltas_per_zero[i]:.6f}",
        'bin': get_bin(m['zero']),
    })

fieldnames = ['n', 't_n', 'nearest_minimum_t', 'delta_n',
              'nearest_null_minimum_t', 'delta_n_null', 'bin']

with open(CSV_FILE, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"  {len(rows)} rows written to {CSV_FILE}")

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY — PRE-REGISTERED RESULTS")
print("=" * 70)

print(f"""
  Pre-registration SHA: 21f74ed
  Riemann zeros:        {len(riemann_zeros)} (range [{riemann_zeros[0]:.2f}, {riemann_zeros[-1]:.2f}])
  zeta_phi minima:      {len(minima_t)} (N={N_TERMS_PRIMARY})

  PRIMARY:
    z_overall =         {z_overall:.4f}
    Decision:           {decision}

  SECONDARY 1 (bin z-scores):""")
for (lo, hi), z in bin_z_scores.items():
    print(f"    [{lo:4.0f},{hi:4.0f}): z = {z:.4f}")
print(f"    Uniform sign:   {uniform}")

print(f"""
  SECONDARY 2 (window sensitivity):
    z at W=1.0:         {z_overall:.4f}
    (see table above for W=0.5, 0.3, 0.2)

  SECONDARY 3 (match rate W=0.3):
    Signal:             {sig_rate_03:.1%}
    Null:               {null_rate_03_mean:.1%}
    z:                  {z_rate:.4f}

  SECONDARY 4 (N_terms sensitivity):
    z at N=5000:        {z_overall:.4f}
    z at N=10000:       {z_N10k:.4f}
""")

print("=" * 70)
print("RUN COMPLETE")
print("=" * 70)
