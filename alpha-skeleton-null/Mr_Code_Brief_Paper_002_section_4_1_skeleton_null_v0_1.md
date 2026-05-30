# Mr Code Brief — Skeleton-Fixed Null Distribution for α⁻¹

**Brief ID:** Skeleton-Fixed Null v0.1
**Date pre-registered:** 30 May 2026
**For:** Paper 2 §4.1 (front-door foundations); also retroactive support for Paper 127 §2 (alpha unification)
**Status:** Pre-registered before computation per Pattern 75 discipline
**Locked seed:** 42 (continuity with Paper 127 §2 null brief)

---

## Background and motivation

Mr Adversary's cycle-1 review of Paper 2 v8.0 (the front-door foundations paper) raised this concern:

> "The null you actually need is over the family e^a − b√c − d + 1/n with icosahedral integers — the family your own formula lives in… A null that fixes the skeleton and varies the integers is a far smaller, far less impressive space, and it is the one that bears on your claim."

Paper 127's already-discharged §2 null samples over the broader **primitive-space** {e, π, √2, √3, √5, integers 1–30} with bounded combinatorial depth. That null answers "does any small closed-form expression in this primitive space land near 137.036?". It is not the same shape as the audit Mr A is asking for, which holds the **functional skeleton** `e^a − b√c − d + 1/n` fixed and varies the five integer parameters (a, b, c, d, n) within bounded ranges.

This brief specifies the skeleton-fixed null.

---

## The claim under audit

The base α formula:

```
α⁻¹ ≈ e⁵ − 6√3 − 1 + 1/66 = 137.036006…
```

agrees with CODATA 2018 (137.035999084) to **~50 ppb relative** (≈ 0.05 ppm). The canonical integer choice is (a, b, c, d, n) = (5, 6, 3, 1, 66), each of which has an icosahedral interpretation:

- a = 5 = E/D! (edges-per-closure-term ratio, §3.4 structure-ratio exponent)
- b = 6 = D! (icosahedral closure quantity)
- c = 3 = D (spatial dimension)
- d = 1 (unity, closure)
- n = 66 = D!·(V−1) (harmonic correction denominator)

**Question:** Within the family `e^a − b√c − d + 1/n` over uniformly random integer parameters in bounded ranges, what fraction of expressions land within 50 ppm of 137.036? If the answer is decisively below 1%, the canonical (5, 6, 3, 1, 66) integer choice is structurally distinguished within the skeleton, not generically achievable.

---

## Pre-registered search space

Each parameter drawn independently from the following closed integer ranges:

| Parameter | Range | Reasoning |
|---|---|---|
| a | {1, 2, ..., 10} | e^a needs to plausibly approach 137 from above; a ≤ 6 is the practical ceiling (e^6 ≈ 403; max negative contribution from {b, c, d} ∈ {1,…,30} is ~194, leaving e^a residue ≥ 209 > 137 for a ≥ 6); a ≤ 10 gives margin and includes the canonical a = 5 |
| b | {1, 2, ..., 30} | Same closure as Paper 127 §4 generic-family null; includes the canonical b = 6 |
| c | {1, 2, ..., 30} | Same closure; includes the canonical c = 3 |
| d | {1, 2, ..., 30} | Same closure; includes the canonical d = 1 |
| n | {2, 3, ..., 200} | n ≥ 2 (excluding n = 1 since 1/1 = 1 collapses to the d term and creates a redundant subspace); upper bound 200 includes the canonical n = 66 with margin |

**Total search space size:** 10 × 30 × 30 × 30 × 199 = **53,730,000 tuples**.

**Sampling:** Uniform random sampling without replacement constraint (with replacement is acceptable given the sample size relative to space size). **N = 10⁶ samples.** Random seed: **42** (locked, pre-registered).

**Match function:**
```
F(a, b, c, d, n) = exp(a) − b·sqrt(c) − d + 1/n
match = |F − 137.036| / 137.036 < 50e-6   (50 ppm tolerance)
```

The target value 137.036 matches Paper 127 §2 null tolerance and is approximately CODATA 2018 (137.035999084) to 1 ppm — close enough that 50 ppm window centred at 137.036 contains CODATA cleanly.

---

## Pre-registered decision rule

Headline result: **observed match fraction p** with **95% binomial upper bound** computed by rule of three (if zero matches) or Wilson interval (if non-zero):

- **p < 0.01 (1%)** → **claim structurally distinguished.** The canonical (5, 6, 3, 1, 66) integer tuple is not generically achievable within the skeleton; Paper 2 §4.1 broader-null open work discharged at skeleton-fixed level. Paper 127 §2 DERIVED-conditional gains additional support beyond the primitive-space null.
- **p ≥ 0.01 (1%)** → **claim is generic at the skeleton-fixed level.** Many random integer choices within the skeleton land near 137.036; the canonical match is not structurally distinguished. Both Paper 2 §4 and Paper 127 §2 status would need reconsideration; the formula's DERIVED status would not survive this audit.

The 1% threshold is the same threshold Paper 127 §2 null pre-registered. The audit is pre-registered as binary at this threshold; intermediate values report the upper bound and let Cliff and Mr Adversary judge.

---

## Secondary reports (also pre-registered)

In addition to the headline unconditional match fraction, Mr Code reports:

1. **Match fraction conditioned on a = 5.** If only a ≈ 5 can plausibly give values near 137, this conditional fraction is the more discriminating null. Sample with a = 5 fixed and (b, c, d, n) varying; report match fraction.

2. **Tighter 50 ppb audit (relative precision 5 × 10⁻⁸).** For matches found at 50 ppm, identify whether any land within 50 ppb of 137.036. The canonical case (5, 6, 3, 1, 66) lands at ~50 ppb; we want to know if any other random tuple does. This is the tighter audit Mr Adversary's cycle-2 review on Paper 127 asked for as future work.

3. **Match list.** For any tuple matching at 50 ppm tolerance, report the tuple, F value, relative precision, and a structural interpretation if the integers have an icosahedral reading. (This is for inspection only; the headline result is the count.)

---

## Reproduction script

The Python script is the canonical implementation. The script is committed alongside this brief BEFORE Mr Code runs it. Filename: `paper_002_s4_1_skeleton_null.py` (also serves Paper 127 retroactively as `paper_127_s2_skeleton_null.py` — same script, same result; the two filenames point to the same artefact).

```python
"""
Skeleton-fixed null distribution for α⁻¹ ≈ e^5 - 6√3 - 1 + 1/66.

Pre-registered per Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md.
Seed: 42. N: 10^6. Search space as specified in brief.

Serves Paper 2 §4.1 (broader-family null) and Paper 127 §2 (retroactive
skeleton-fixed support for DERIVED-conditional).
"""
import numpy as np
import math

# Pre-registered constants
SEED = 42
N = 1_000_000
TARGET = 137.036
TOL_PPM = 50e-6   # 50 ppm relative match window
TOL_PPB = 50e-9   # 50 ppb relative match window (secondary, tighter)

# Pre-registered search space bounds
A_RANGE = (1, 10)      # inclusive at both ends
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
# Use numpy for vectorised computation
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
print(f"Search space: 10 × 30 × 30 × 30 × 199 = {10*30*30*30*199:,} tuples")
print(f"")
print(f"--- PRIMARY (50 ppm) ---")
print(f"Matches: {n_match_ppm} / {N}")
print(f"Match fraction: {n_match_ppm / N:.2e}")
if n_match_ppm == 0:
    print(f"95% upper bound (rule of three): {3.0/N:.2e}")
else:
    # Wilson 95% upper
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
print(f"--- SECONDARY: match list (if any) ---")
if n_match_ppm > 0:
    match_idx = np.where(match_ppm)[0]
    for i in match_idx[:50]:
        ai, bi, ci, di, ni = a[i], b[i], c[i], d[i], n[i]
        Fi = F[i]
        rdi = rel_dev[i]
        print(f"  (a, b, c, d, n) = ({ai}, {bi}, {ci}, {di}, {ni}) → F = {Fi:.6f}, rel dev = {rdi:.2e}")
    if len(match_idx) > 50:
        print(f"  ... and {len(match_idx) - 50} more")
else:
    print("  No matches at 50 ppm tolerance.")
print(f"")
print(f"=== Decision ===")
if n_match_ppm / N < 0.01:
    print("DECISION: p < 1% — claim structurally distinguished within skeleton.")
    print("→ Paper 2 §4.1 broader-family null discharged at skeleton-fixed level.")
    print("→ Paper 127 §2 DERIVED-conditional gains additional skeleton-fixed support.")
else:
    print("DECISION: p ≥ 1% — claim is GENERIC at skeleton-fixed level.")
    print("→ Paper 2 §4 and Paper 127 §2 status need reconsideration.")
    print("→ DERIVED-conditional does not survive this audit.")
```

---

## Commitment checklist (pre-registration discipline)

Before Mr Code runs:

- [x] Brief written with full search space, match criterion, decision rule
- [x] Locked seed: 42
- [x] N pre-registered: 10⁶
- [x] Match tolerance pre-registered: 50 ppm (with secondary 50 ppb)
- [x] Decision threshold pre-registered: p < 1% → structurally distinguished
- [x] Secondary analyses pre-registered: conditional on a=5; tighter 50 ppb; match list
- [x] Script committed alongside brief in same conversation context, time-stamped
- [ ] Mr Code runs script in fresh context
- [ ] Result reported and incorporated into Paper 2 §4.1 (and Paper 127 §2 v3.0)

---

## Notes on use

1. **One script serves both papers.** When Mr Code runs this, the result discharges:
   - Paper 2 v8.1 §4.1's "broader e^a − b√c − d + 1/n family null" open-work claim
   - Paper 127 v3.0 §2's skeleton-fixed audit (alongside the already-discharged primitive-space audit)

   The result will appear in Paper 2 v8.1 §4.1 (if cycle-2 awaits this) or v8.2 (if v8.1 ships with brief committed but result pending).

2. **Pre-registration is unconditional.** This brief and the script must not be edited after commitment, regardless of how the result turns out. If the audit fails (p ≥ 1%), Paper 2 §4 and Paper 127 §2 are reined back — Pattern 100 retirement of the DERIVED-conditional applies, and Paper 14's tension-theory retirement is the model.

3. **The secondary analyses are also pre-registered.** Mr Code should report all three (primary 50 ppm, conditional on a=5, tighter 50 ppb) even if the headline primary result is decisive. Cliff judges what each result implies; the brief commits Mr Code to producing the full report.

🐕☕⬡
