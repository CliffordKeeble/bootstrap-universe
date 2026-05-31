# Sub F follow-up findings — RvM-random draw distribution

**Date**: 31 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 31 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `078cc38` BEFORE any new seed was sampled.
**Verdict**: **REJECT SUSTAIN** (K = 4.45σ ≥ 3).

Per brief stop-on-fail: **pause Paper 150 v2.3 drafting. Flag for CinC adjudication. The Sub F retraction needs reframing — the probe IS partially L-zero-specific in z-score terms, even though it isn't in effect-size terms.**

---

## Methodology

Identical to Sub F (commit `d03517b`) except the RvM-random seed varies across 10 pre-registered values:
```
[20260530, 20260531, 20260601, 20260602, 20260603,
 20260604, 20260605, 20260606, 20260607, 20260608]
```
Seed 20260530 is the Sub F original.

All other parameters frozen: probe `|Re² − 5·Im²|`, N_terms = 5000, α = φ, dt = 0.008, W = 1.0, N_null = 1000, mpmath dps = 15, N_target = 2792.

The probe (Re, Im) computed once on t ∈ [1.0, 3335.7] (416,837 grid points). Reused across all 10 RvM draws. Minima identified once (3942 total, density 1.18/unit).

## Per-seed results

| seed | n_targets_in_range | n_matched | signal_mean | null_mean | null_std | z | ε |
|---|---|---|---|---|---|---|---|
| 20260530 | 2792 | 2697 | 0.3225 | 0.4282 | 0.0104 | −10.175 | 0.2467 |
| 20260531 | 2792 | 2713 | 0.3128 | 0.4282 | 0.0108 | −10.692 | 0.2695 |
| 20260601 | 2792 | 2700 | 0.3216 | 0.4280 | 0.0104 | −10.197 | 0.2488 |
| 20260602 | 2792 | 2704 | 0.3204 | 0.4278 | 0.0102 | −10.572 | 0.2516 |
| 20260603 | 2792 | 2709 | 0.3214 | 0.4279 | 0.0106 | −10.018 | 0.2489 |
| 20260604 | 2792 | 2707 | 0.3251 | 0.4279 | 0.0105 | −9.835 | 0.2401 |
| 20260605 | 2792 | 2699 | 0.3162 | 0.4280 | 0.0105 | −10.622 | 0.2617 |
| 20260606 | 2792 | 2718 | 0.3066 | 0.4281 | 0.0105 | −11.614 | 0.2836 |
| 20260607 | 2792 | 2704 | 0.3179 | 0.4282 | 0.0108 | −10.179 | 0.2578 |
| 20260608 | 2792 | 2699 | 0.3225 | 0.4281 | 0.0105 | −10.031 | 0.2466 |

CSV: [sub_f_followup_rvm_distribution.csv](sub_f_followup_rvm_distribution.csv).

## Distribution statistics

### z-score across 10 RvM draws

| Statistic | Value |
|---|---|
| mean | **−10.39** |
| std | **0.488** |
| min | −11.61 |
| max | −9.83 |
| median | −10.19 |
| range | 1.78 |

Notably: the Sub F single-draw value (seed 20260530, z = −10.175 here, slightly different from the original z = −9.89 reported in Sub F due to a small change in the probe grid extent — see "Note on Sub F numerical reproducibility" below) **is close to the median, not particularly low**.

### Effect size ε across 10 RvM draws

| Statistic | Value |
|---|---|
| mean | **0.256** |
| std | **0.012** |
| min | 0.240 |
| max | 0.284 |
| median | 0.250 |

## Comparison to Riemann control

Riemann control from Sub F: z = **−12.56**, ε = **0.234**.

**z-space comparison:**
> K = (12.56 − 10.39) / 0.488 = **4.45**

Riemann |z| = 12.56 is **above max(RvM |z|) = 11.61**. **0 of 10** RvM draws reach Riemann's |z|.

**ε-space comparison:**
> Riemann ε = 0.234 is **below all 10 RvM ε** (min 0.240).

This is the nuance the pre-reg didn't explicitly anticipate: **the Riemann control has a SMALLER effect size than the RvM-random distribution, but a LARGER |z|-score, because the Riemann null has tighter std.**

| | Riemann | RvM mean |
|---|---|---|
| signal_mean | 0.328 | 0.320 |
| null_mean | 0.428 | 0.428 |
| null_std | **0.0079** | **0.0105** |
| ε | 0.234 | 0.256 |
| \|z\| | 12.56 | 10.39 |

**The null_std is the key difference.** Riemann's null_std = 0.0079 is ~25% tighter than RvM's 0.0105. This tightening, applied to roughly the same signal_mean − null_mean gap (~0.10), gives Riemann a larger |z| even though its ε is smaller.

## Pre-registered verdict

Per pre-registration:
> K-σ separation: K = (|z_Riemann| − mean(|z_RvM|)) / std(|z_RvM|)

K = **4.45**.

Per pre-registered thresholds:
- K < 1 → CONFIRM SUSTAIN
- 1 ≤ K < 3 → PARTIAL DISTINCTION
- K ≥ 3 → **REJECT SUSTAIN**

**Verdict: REJECT SUSTAIN.**

Per brief stop-on-fail:
> "If REJECT SUSTAIN (K ≥ 3): pause Paper 150 v2.3 drafting. The retraction may need to be reframed or partially withdrawn. Flag for CinC adjudication."

I **pause v2.3 drafting** and flag for CinC.

## What this means — the structural story behind the verdict

The Sub F single-draw RvM result (z = −9.89) was a low outlier (lower than 7 of 10 in the distribution). The distribution mean |z| is 10.39, which is meaningfully but not catastrophically below Riemann's 12.56 (K = 4.5σ).

Critical nuance: **the K-σ difference is driven by null_std, not by signal_mean**. Riemann zeros have GUE-like spacing (Montgomery's pair correlation), giving them a more uniform 1D arrangement than RvM-random (Poisson) at small scales. This GUE-like uniformity makes the random-minima null *more concentrated* (tighter std) when the target is Riemann zeros than when the target is RvM-random. The Riemann null_std = 0.0079 vs RvM null_std = 0.0105 — a 25% gap.

In **effect-size** terms, the probe matches RvM-random points *slightly better* than Riemann zeros (ε 0.256 vs 0.234). The probe is NOT detecting Riemann zeros more effectively in magnitude.

In **z-score** terms, the probe's detection of Riemann zeros is 4.5σ more significant than its detection of RvM-random — but only because the Riemann null is tighter due to GUE spacing.

**This is partial L-zero specificity**, but of a subtle kind: the probe's effect size is not L-zero-specific (RvM gives slightly larger ε), but the statistical significance via z-score IS L-zero-specific due to the GUE-like spacing of Riemann zeros tightening the null variance.

## Note on Sub F numerical reproducibility

The Sub F seed-20260530 result was z = **−9.89** (per `d03517b`). This follow-up run on the same seed gives z = **−10.18**. The 0.3-difference is due to a minor change in the probe grid: Sub F used `t_max = max(target) + 5`, where max(target) for RvM with seed 20260530 was 3329.65; this follow-up uses `t_max_grid = T_Riemann + 10 = 3335.69`, so the grid extends slightly further and includes slightly more minima.

The two numbers are within sampling noise; the qualitative result (RvM ≈ Riemann at single-draw level) is unchanged. The shift IS large enough (~3% in z) to affect the strict comparison with the Sub F-cached Riemann control, but K = 4.45 is robustly ≥ 3 either way. The mean of 10 RvM draws here is −10.39, and the Sub F single draw at −9.89 sat 1σ below this mean — slightly low but not extreme.

## Implications for Paper 150 v2.3

(Pending CinC adjudication per stop-on-fail.)

The Sub F verdict needs revision. The right framing is now:

1. **The probe detects L-zero structure at higher statistical significance than RvM-random density-matched points, by 4.45σ.** (K-σ in z-score.)
2. **The probe's effect size against L-zeros is slightly SMALLER than against RvM-random.** (ε comparison: 0.234 vs mean 0.256.)
3. **The discriminating quantity is the null variance, which is tighter for L-zeros due to their GUE-like pair correlation.** The probe doesn't "find" L-zeros better; the L-zeros' inherent spacing structure tightens the null.

This is **partial L-zero specificity**, mediated by target spacing, not by effect size. Paper 150 v2.3 should articulate this nuance clearly:

- Sub F's headline "probe detects density, not zeros" was correct in **effect-size** terms.
- Sub F follow-up's headline "Riemann distinguished from RvM at 4.5σ" is correct in **z-score** terms.
- The mechanism: L-zeros' GUE-like spacing makes their detection statistically more significant via tighter null variance, not via higher effect size.

The detection IS partially L-zero-specific, but in a way that none of the Paper 150 v2.0, v2.1, v2.2 mechanism stories captured. v2.3 must be written carefully to:
- Acknowledge the K = 4.5σ z-space distinction.
- Acknowledge the ε-space null result.
- Explain the GUE-like spacing as the mediator.
- Avoid over- or under-claiming.

## Honest pre-vs-post comparison

My pre-registered prior: "leans CONFIRM SUSTAIN... K ≈ 1.8σ guess." Reality: K = 4.45, **REJECT SUSTAIN**.

**This is the third time in this brief sequence my pre-registered prior was wrong**:
- Sub E: predicted Flag 1 "d is decoration" would be CONFIRMED at ≥80% sharing; reality was 30% (REFUTED).
- Sub F: predicted RvM COLLAPSE; reality was PARTIAL near SUSTAIN.
- Sub F follow-up: predicted CONFIRM SUSTAIN; reality is REJECT SUSTAIN.

The pattern: I keep underestimating L-zero specificity. My readings have been over-generalising the "detection is generic" framing. The empirical answer in each case has pulled back toward "there's something L-zero-specific going on" — though never in the simple way v2.0 or v2.1 proposed.

The pre-registration discipline keeps working: claims pre-committed, tests run cleanly, refuted priors stand refuted, no rescue attempted.

## Honest limitations

- Only 10 seeds. Larger ensemble (50-100 seeds) would tighten the K-σ estimate and the std estimate. With 10 seeds, std(|z|) = 0.488 has its own sampling uncertainty.
- The GUE-like-spacing-tightens-null hypothesis is an *interpretation* of why null_std differs between Riemann and RvM. A direct test would be: compute the null with target being a random shuffle of the actual Riemann zero positions (which preserves their spacing structure) vs RvM. If shuffled-Riemann gives Riemann-like null_std, the GUE-spacing hypothesis is confirmed.
- The Sub F → Sub F-follow-up shift (z=−9.89 vs z=−10.18 at seed 20260530) is a small discrepancy due to grid extent. It doesn't change qualitative reading but I should flag it.

## Pattern flags

- **Pattern 19 (adversary)**: Mr A's NULL catch fully sustained — the single-draw was indeed not a distribution. The fix produced REJECT SUSTAIN, reversing the previous strong claim. Adversary catches improving the result is what the loop is for.
- **Pattern 75 (null)**: pre-registered random null discipline. The RvM distribution IS a null over targets; the resulting distribution shows L-zeros are not in the bulk of this null at K = 4.5σ.
- **Pattern 39 (DERIVED vs OBSERVED)**: 10-seed distribution is OBSERVED. The GUE-spacing-tightens-null interpretation is FRAMEWORK-INTERPRETATION, not derived.

## What I will NOT do per brief stop-on-fail

- Draft Paper 150 v2.3. CinC adjudicates first.
- Attempt to reframe the v2.0/v2.1 mechanism stories. The CinC adjudication may produce a NEW mechanism reading (GUE-spacing-mediated detection), but I do not draft it pre-CinC.
- Re-interpret the verdict. K = 4.45, REJECT SUSTAIN, stands.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — the 10 seeds and K-σ thresholds, pre-committed
- [sub_f_followup.py](sub_f_followup.py) — implementation
- [sub_f_followup_rvm_distribution.csv](sub_f_followup_rvm_distribution.csv) — per-seed results
- [findings_paper_150_sub_f_followup_rvm_distribution.md](findings_paper_150_sub_f_followup_rvm_distribution.md) — this report

## Compute

- Probe Re/Im: ~45 s (one-shot, reused).
- Minima identification: ~1 s.
- Per RvM seed (matching + 1000-trial MC null): ~28-31 s each.
- Total: ~5 min, exactly the brief's estimate.
