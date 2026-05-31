# Sub E findings — mechanism verification for Paper 150 v2.1

**Date**: 30 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `cfd7599` BEFORE any probe evaluation.
**Verdict**: **DECORATION-REFUTED** (mean off-diagonal sharing = 29.83%, below 40% threshold).

Per brief's stop-on-fail clause: **pause Paper 150 v2.1 drafting and flag for CinC.** The "d is decoration" mechanism reading from Sub C-1 Flag 1/2 was an overclaim. The actual mechanism is more nuanced.

---

## Methodology

Per pre-registration:

- **Probe**: `P_d(t) = sqrt(|Re²(Z_φ(t)) − d·Im²(Z_φ(t))|)`, Paper 150 v2.0 Fejér-weighted golden-angle Dirichlet sum, N_terms = 5000.
- **Grid**: t ∈ (1, 1000] at Δt = 0.01 → 99,900 points.
- **Local minimum**: three-point strict (`P[i-1] > P[i] < P[i+1]`).
- **Top-N selection**: greedy lowest-first, with non-overlap window 2.0 t-units between selected minima.
- **Sharing tolerance**: ±0.02 (= 2·Δt) per brief.

Re(Z_φ) and Im(Z_φ) computed *once* (since they don't depend on d); P_d(t) formed for each d in {3, 5, 7, 13}.

## Results

Per-d minima statistics:

| d | total local minima | top-100 lowest value | top-100 highest value |
|---|---|---|---|
| 3 | 2437 | 0.00456 | 0.06583 |
| 5 | 2449 | 0.00256 | 0.07094 |
| 7 | 2452 | 0.00592 | 0.07311 |
| 13 | 2454 | 0.00796 | 0.09128 |

The four probes have very similar total minima counts (~2450), suggesting the *density* of minima is d-independent — but the *positions* are not.

## Sharing matrix (% of row's 100 lowest minima within 0.02 of any column-d's minima)

|          | d'=3 | d'=5 | d'=7 | d'=13 |
|----------|-----:|-----:|-----:|------:|
| **d=3**  | 100.0% | 32.0% | 31.0% | 28.0% |
| **d=5**  | 32.0% | 100.0% | 32.0% | 29.0% |
| **d=7**  | 31.0% | 32.0% | 100.0% | 27.0% |
| **d=13** | 28.0% | 29.0% | 27.0% | 100.0% |

**Mean off-diagonal sharing: 29.83%.**

## Chance baseline

Computed via Monte Carlo (1000 trials of two independent uniform random selections of 100 t-values in (1, 1000] with non-overlap 2):
> Chance sharing = 0.39% ± 0.60%, 95% CI [0%, 2%], max observed 4%.

So 29.83% is **~49σ above chance**. There IS substantial shared structure — but it is far below the ≥80% threshold the brief specified for "decoration."

## Verdict per pre-registered thresholds

- Threshold for DECORATION-CONFIRMED: ≥80%. Not met.
- Threshold for DECORATION-PARTIAL: 40-80%. Not met.
- Threshold for DECORATION-REFUTED: <40%. **MET (29.83%)**.

**Verdict: DECORATION-REFUTED.**

Per brief stop-on-fail clause:
> "If DECORATION-REFUTED, pause and flag for CinC. v2.1's mechanism section needs substantial rewriting and possibly further investigation. Don't draft v2.1 until CinC has adjudicated."

## What the data actually show — honest reading

My Flag 1/2 reading in Sub C-1 was: "the probe |Re² − d·Im²| doesn't depend on d through the minima locations — only through which minima are below the 5th-percentile threshold." Sub E refutes this claim *in degree*.

The actual structure:
1. The probe minima density is essentially d-independent (~2450 minima per d in the same range).
2. ~30% of top-100 minima are **structural / d-invariant** — far above chance (0.4%), but only a minority.
3. ~70% of top-100 minima are **d-specific** — the discriminant shifts where the minima are.
4. The Sub C-1 matrix-uniform effect size (ε ≈ 0.29 across all 16 cells) therefore reflects: the L-zero distribution correlates similarly with BOTH the structural-30% and the d-specific-70% — i.e., **the L-zeros are dense enough on the t-axis that both classes of probe minima land near L-zeros at similar rates.**

This is a *real and interesting* structural fact: ~30% of the probe minima are dictated by the equidistributed Dirichlet sum independent of d, and ~70% are dictated by the d in the bilinear form. Both classes correlate with L-zero positions at similar effect sizes — which is why Sub C-1's matrix is uniform.

## Implications for Paper 150 v2.1 — what the mechanism section should say

(Pending CinC adjudication per the stop-on-fail clause.)

A more honest mechanism reading is now possible:

1. **The probe's minima come in two classes**:
   - *Structural class* (~30%): minima that are d-invariant. These arise from features of `Z_φ(t)` itself where Re² and Im² both reach near-zero simultaneously, so |Re² − d·Im²| is small for any d.
   - *D-dependent class* (~70%): minima that depend on d, occurring where Re² ≈ d·Im² (the null cone of the d-specific bilinear form). Different d's give different intersection patterns with the trajectory.
2. **Both classes correlate with L-zeros at similar rates** because the L-zero distribution on the critical line is approximately uniform-dense at the t-scale of probe minima spacing (~0.4/unit for L-zeros vs ~2.5/unit for probe minima).
3. **The uniform Sub C-1 effect size (ε ≈ 0.29)** is NOT explained by "d is decoration." It is explained by the L-zeros being dense enough to land near any reasonably-dense set of probe minima, whether structural or d-specific.

This is a more nuanced mechanism reading than my Flag 1 suggested. v2.1 (if drafted post-CinC adjudication) should reflect this.

## Honest limitations

- The sharing test uses a tight tolerance (±0.02) requiring near-exact coincidence. A looser tolerance (±0.05 or ±0.1) might give higher apparent sharing percentages — but that would test "approximately the same" rather than "the same minimum."
- The 100-lowest selection with non-overlap window 2 was the brief's specification. Different choices (200 lowest? non-overlap window 1?) might give different sharing percentages.
- The shared 30% might decompose further into "exactly d-invariant" and "structurally similar but slightly shifted by d." A finer analysis could distinguish these.
- The verdict applies to the specific (d ∈ {3, 5, 7, 13}, t ∈ (1, 1000], α = φ) configuration of Sub C-1. A different α (e.g., √2 from Paper 150 v2.0 Test 2) might give different sharing patterns.

## Pattern flags

- **Pattern 19 (adversary)**: my Flag 1 reading made a structural claim that was directly testable. Sub E ran the test and the claim failed in degree (30% not 80%). This is exactly the discipline Mr Adversary's review structure demands — making claims that can be checked. **The honest verdict is REFUTED**, not "qualified."
- **Pattern 39 (DERIVED vs OBSERVED)**: the probe construction is DERIVED. The minima locations, sharing percentages, and chance baseline are all OBSERVED (computed exactly from the probe).
- **Pattern 75 (null)**: chance baseline computed via 1000-trial Monte Carlo. The 30% sharing is far above the 0.4% chance baseline, so the result is not "no shared structure" — it's "less shared structure than my analytical reading predicted."

## What I will NOT do (per brief stop-on-fail)

- Draft Paper 150 v2.1's mechanism section. The brief explicitly says wait for CinC.
- Reinterpret the verdict to avoid the REFUTED label. The data are what they are: 29.83% < 40%.
- Run a different version of the test post-hoc (e.g., wider tolerance) to convert REFUTED → PARTIAL. That would violate pre-registration.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — pre-registered methodology
- [sub_e.py](sub_e.py) — implementation
- [sharing_matrix.csv](sharing_matrix.csv) — 4×4 results
- [minima_per_d.csv](minima_per_d.csv) — top-100 minima per d for traceability
- [findings_paper_150_sub_e.md](findings_paper_150_sub_e.md) — this report

## Compute

~7 seconds total. Far below the brief's 1-2 hour budget — single Z_φ evaluation reused across d's, minima detection is vectorised.
