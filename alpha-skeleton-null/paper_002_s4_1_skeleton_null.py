"""
Skeleton-fixed null distribution for α⁻¹ ≈ e^5 - 6√3 - 1 + 1/66.

Pre-registered per Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md.
Seed: 42. N: 10^6. Search space as specified in brief.

Serves Paper 2 §4.1 (broader-family null) and Paper 127 §2 (retroactive
skeleton-fixed support for DERIVED-conditional).

Run: python3 paper_002_s4_1_skeleton_null.py
"""
import numpy as np

# Pre-registered constants
SEED = 42
N = 1_000_000
TARGET = 137.036
TOL_PPM = 50e-6   # 50 ppm relative match window (primary)
TOL_PPB = 50e-9   # 50 ppb relative match window (secondary, tighter)

# Pre-registered search space bounds (inclusive at both ends)
A_RANGE = (1, 10)
B_RANGE = (1, 30)
C_RANGE = (1, 30)
D_RANGE = (1, 30)
N_RANGE = (2, 200)     # n=1 excluded (collapses 1/n into d term)

rng = np.random.default_rng(SEED)

# Uniform random sample over the search space
a = rng.integers(A_RANGE[0], A_RANGE[1] + 1, size=N)
b = rng.integers(B_RANGE[0], B_RANGE[1] + 1, size=N)
c = rng.integers(C_RANGE[0], C_RANGE[1] + 1, size=N)
d = rng.integers(D_RANGE[0], D_RANGE[1] + 1, size=N)
n = rng.integers(N_RANGE[0], N_RANGE[1] + 1, size=N)

# Compute F = exp(a) - b*sqrt(c) - d + 1/n
F = np.exp(a.astype(np.float64)) - b * np.sqrt(c.astype(np.float64)) - d + 1.0 / n.astype(np.float64)

# Relative deviation from target
rel_dev = np.abs(F - TARGET) / TARGET

# Primary: matches at 50 ppm
match_ppm = rel_dev < TOL_PPM
n_match_ppm = int(np.sum(match_ppm))

# Secondary: matches at 50 ppb (tighter)
match_ppb = rel_dev < TOL_PPB
n_match_ppb = int(np.sum(match_ppb))

# Conditioned on a = 5
a5_mask = a == 5
n_a5 = int(np.sum(a5_mask))
n_match_ppm_a5 = int(np.sum(match_ppm & a5_mask))

# Report
print(f"=== Skeleton-Fixed Null Result ===")
print(f"Brief: Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md")
print(f"Seed: {SEED}; N: {N}; Target: {TARGET}; Tolerance: {TOL_PPM*1e6:.0f} ppm")
print(f"Search space: 10 x 30 x 30 x 30 x 199 = {10*30*30*30*199:,} tuples")
print(f"")
print(f"--- PRIMARY (50 ppm) ---")
print(f"Matches: {n_match_ppm} / {N}")
print(f"Match fraction: {n_match_ppm / N:.2e}")
if n_match_ppm == 0:
    print(f"95% upper bound (rule of three): {3.0/N:.2e}")
else:
    p_hat = n_match_ppm / N
    z = 1.96
    denom = 1 + z**2/N
    centre = p_hat + z**2/(2*N)
    half = z * np.sqrt(p_hat*(1-p_hat)/N + z**2/(4*N**2))
    upper = (centre + half) / denom
    print(f"95% upper bound (Wilson): {upper:.2e}")
print(f"Pre-registered threshold: 1.00e-02 (1%)")
print(f"")
print(f"--- SECONDARY: conditioned on a = 5 ---")
print(f"Samples with a = 5: {n_a5}")
print(f"Matches within 50 ppm: {n_match_ppm_a5}")
if n_a5 > 0:
    print(f"Conditional match fraction: {n_match_ppm_a5 / n_a5:.2e}")
print(f"")
print(f"--- SECONDARY: tighter 50 ppb ---")
print(f"Matches: {n_match_ppb} / {N}")
print(f"Match fraction: {n_match_ppb / N:.2e}")
print(f"")
print(f"--- SECONDARY: match list (if any, up to 50 shown) ---")
if n_match_ppm > 0:
    match_idx = np.where(match_ppm)[0]
    for i in match_idx[:50]:
        ai, bi, ci, di, ni = int(a[i]), int(b[i]), int(c[i]), int(d[i]), int(n[i])
        Fi = float(F[i])
        rdi = float(rel_dev[i])
        print(f"  (a, b, c, d, n) = ({ai}, {bi}, {ci}, {di}, {ni}) -> F = {Fi:.6f}, rel dev = {rdi:.2e}")
    if len(match_idx) > 50:
        print(f"  ... and {len(match_idx) - 50} more")
else:
    print("  No matches at 50 ppm tolerance.")
print(f"")
print(f"=== Decision ===")
if n_match_ppm / N < 0.01:
    print("DECISION: p < 1% -- claim structurally distinguished within skeleton.")
    print("-> Paper 2 §4.1 broader-family null discharged at skeleton-fixed level.")
    print("-> Paper 127 §2 DERIVED-conditional gains additional skeleton-fixed support.")
else:
    print("DECISION: p >= 1% -- claim is GENERIC at skeleton-fixed level.")
    print("-> Paper 2 §4 and Paper 127 §2 status need reconsideration.")
    print("-> DERIVED-conditional does not survive this audit.")
