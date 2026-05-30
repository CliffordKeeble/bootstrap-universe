"""
Exhaustive skeleton-fixed enumeration for alpha^-1 ~ e^5 - 6*sqrt(3) - 1 + 1/66.

Confirmatory secondary to paper_002_s4_1_skeleton_null.py. The pre-registered
Monte-Carlo headline (N=10^6, seed 42) stands as committed and unchanged; this
script computes the EXACT match count over the entire pre-registered search
space (53,730,000 tuples), removing seed dependence and making the tight 50 ppb
window meaningful (at 10^6 samples the 50 ppb count is undersampling-limited).

Same search space, same target, same tolerances as the brief. Writes the full
50 ppm match list to matches_50ppm.csv alongside this script, with both the
brief TARGET and the actual CODATA 2018 reference, and a distinct-value id that
collapses the b*sqrt(c) surd-representation degeneracy (6*sqrt(3) = 3*sqrt(12)
= 2*sqrt(27) are one value written three ways).

Run: python3 paper_002_s4_1_skeleton_null_exhaustive.py
"""
import csv
import os
import numpy as np

TARGET = 137.036                # brief target (50 ppm window centred here)
CODATA_2018 = 137.035999084     # actual reference value
TOL_PPM = 50e-6
TOL_PPB = 50e-9

A_RANGE = range(1, 11)          # 1..10
N_LO, N_HI = 2, 200            # n in 2..200

b = np.arange(1, 31, dtype=np.float64)
c = np.arange(1, 31, dtype=np.float64)
d = np.arange(1, 31, dtype=np.float64)
n = np.arange(N_LO, N_HI + 1, dtype=np.float64)
inv_n = 1.0 / n

bb, cc = np.meshgrid(b, c, indexing="ij")          # (30,30)
bsqrtc = (bb * np.sqrt(cc)).ravel()                 # (900,)
b_idx = bb.ravel().astype(int)
c_idx = cc.ravel().astype(int)

total_space = len(list(A_RANGE)) * 30 * 30 * 30 * len(n)

n_ppm = 0
n_ppb = 0
rows = []                       # (a,b,c,d,n,bsqrtc,F,rel_t,rel_cod)
per_a_ppm = {}

for a in A_RANGE:
    ea = np.exp(float(a))
    F = ea - bsqrtc[:, None, None] - d[None, :, None] + inv_n[None, None, :]
    rel_t = np.abs(F - TARGET) / TARGET
    m_ppm = rel_t < TOL_PPM
    cnt = int(m_ppm.sum())
    per_a_ppm[a] = cnt
    n_ppm += cnt
    n_ppb += int((np.abs(F - TARGET) / TARGET < TOL_PPB).sum())
    if cnt:
        ii, jj, kk = np.where(m_ppm)
        for i, j, k in zip(ii, jj, kk):
            Fi = float(F[i, j, k])
            rows.append((a, int(b_idx[i]), int(c_idx[i]), int(d[j]), int(n[k]),
                         float(bsqrtc[i]), Fi,
                         abs(Fi - TARGET) / TARGET,
                         abs(Fi - CODATA_2018) / CODATA_2018))

# Distinct-value id: collapse surd-representation degeneracy by rounded F.
rows.sort(key=lambda r: r[7])                        # by rel-to-target
fkey = {}
for r in rows:
    k = round(r[6], 7)
    if k not in fkey:
        fkey[k] = len(fkey) + 1
n_distinct_ppm = len(fkey)

here = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(here, "matches_50ppm.csv")
with open(csv_path, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["a", "b", "c", "d", "n", "b_sqrt_c", "F",
                "rel_dev_target", "rel_dev_codata",
                "within_50ppb", "distinct_value_id"])
    for r in rows:
        a, bi, ci, di, ni, bs, Fi, rt, rc = r
        w.writerow([a, bi, ci, di, ni, f"{bs:.9f}", f"{Fi:.9f}",
                    f"{rt:.3e}", f"{rc:.3e}",
                    int(rt < TOL_PPB), fkey[round(Fi, 7)]])

# Canonical tuple, for the record
F0 = np.exp(5) - 6 * np.sqrt(3) - 1 + 1.0 / 66
print("=== Exhaustive Skeleton-Fixed Enumeration ===")
print(f"Search space (exact, full grid): {total_space:,} tuples")
print(f"Target: {TARGET}   CODATA 2018: {CODATA_2018}")
print("")
print("--- PRIMARY (50 ppm), EXACT ---")
print(f"Raw tuples: {n_ppm} / {total_space}  (fraction {n_ppm/total_space:.3e})")
print(f"Distinct values: {n_distinct_ppm}  (fraction {n_distinct_ppm/total_space:.3e})")
print(f"Pre-registered threshold: 1.00e-02 (1%)")
print("")
print("--- Matches by a (effective dimensionality) ---")
print("  " + ", ".join(f"a={a}:{per_a_ppm[a]}" for a in A_RANGE))
print("")
print("--- TIGHTER (50 ppb), EXACT ---")
ppb_rows = [r for r in rows if r[7] < TOL_PPB]
ppb_distinct = sorted({round(r[6], 8) for r in ppb_rows})
print(f"Raw tuples: {n_ppb}   Distinct values: {len(ppb_distinct)}")
for r in ppb_rows:
    a, bi, ci, di, ni, bs, Fi, rt, rc = r
    print(f"  (a,b,c,d,n)=({a},{bi},{ci},{di},{ni})  F={Fi:.9f}  "
          f"rel_target={rt:.2e}  rel_codata={rc:.2e}")
print("")
print(f"Canonical (5,6,3,1,66): F={F0:.9f}  "
      f"rel_target={abs(F0-TARGET)/TARGET:.2e}  rel_codata={abs(F0-CODATA_2018)/CODATA_2018:.2e}")
print(f"CSV written: {csv_path}")
