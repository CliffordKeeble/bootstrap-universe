# Sub G findings — target gap-shuffle null

**Date**: 31 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 31 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `a3a04a0` BEFORE any gap-shuffle was performed.
**Verdict**: **NEAREST-NEIGHBOR SUFFICIENT (H1)** — with a clean and surprising twist (see "Findings beyond verdict" below).

This is the **last** empirical step before Paper 150 v2.3 per the brief's stop-on-clean-verdict clause. v2.3 will be drafted from this finding.

---

## Methodology

Identical to Sub F follow-up except the target type:

- **Target**: gap-shuffled Riemann zeros. For each pre-registered seed, the empirical gaps `g_i = t_{i+1} - t_i` of the first 2792 Riemann zeros are randomly permuted (numpy `default_rng(seed).permutation(2791)`), then the target positions are reconstructed as `t_1, t_1 + Σ g_{σ(1..i-1)}, ...`. The result preserves:
  - First point (= first Riemann zero)
  - Last point (= last Riemann zero)
  - Multiset of nearest-neighbor gaps
  - Overall density
  while destroying:
  - Order of gaps (long-range correlations)
  - All n-point correlations beyond pair
  - Any L-function-specific positional information

The gap-shuffle preservation was verified before running (sorted gap arrays of original and shuffled are identical; first/last points unchanged).

All other parameters frozen (probe Re/Im computed once, reused; MC null seed 42, 1000 trials).

10 pre-registered seeds (same as Sub F follow-up):
`[20260530, 20260531, 20260601, 20260602, 20260603, 20260604, 20260605, 20260606, 20260607, 20260608]`

## Per-seed results

| seed | n_matched | signal_mean | null_mean | null_std | z | ε |
|---|---|---|---|---|---|---|
| 20260530 | 2710 | 0.3155 | 0.4275 | 0.00793 | −14.120 | 0.2637 |
| 20260531 | 2702 | 0.3171 | 0.4283 | 0.00816 | −13.623 | 0.2602 |
| 20260601 | 2711 | 0.3117 | 0.4288 | 0.00822 | −14.245 | 0.2732 |
| 20260602 | 2708 | 0.3162 | 0.4290 | 0.00811 | −13.877 | 0.2635 |
| 20260603 | 2692 | 0.3179 | 0.4284 | 0.00798 | −13.833 | 0.2579 |
| 20260604 | 2697 | 0.3171 | 0.4283 | 0.00805 | −13.824 | 0.2613 |
| 20260605 | 2701 | 0.3198 | 0.4292 | 0.00829 | −13.208 | 0.2552 |
| 20260606 | 2698 | 0.3174 | 0.4292 | 0.00803 | −13.917 | 0.2597 |
| 20260607 | 2708 | 0.3127 | 0.4288 | 0.00823 | −14.108 | 0.2706 |
| 20260608 | 2704 | 0.3173 | 0.4286 | 0.00817 | −13.607 | 0.2599 |

CSV: [sub_g_target_shuffle_results.csv](sub_g_target_shuffle_results.csv).

## Aggregate statistics

### z distribution (10 shuffles)

| Statistic | Value |
|---|---|
| mean | **−13.84** |
| std | **0.287** |
| min | −14.25 |
| max | −13.21 |
| median | −13.85 |
| range | 1.04 |

### Effect size ε distribution

| Statistic | Value |
|---|---|
| mean | **0.263** |
| std | 0.005 |
| min | 0.255 |
| max | 0.273 |

### Null_std distribution

| Statistic | Value |
|---|---|
| mean | **0.00800** |
| std | 0.00012 |
| min | 0.00793 |
| max | 0.00829 |

## Three-way comparison

| Quantity | Riemann (Sub F) | Gap-shuffle (Sub G, 10 seeds) | RvM-Poisson (Sub F follow-up, 10 seeds) |
|---|---:|---:|---:|
| signal_mean | 0.328 | 0.316 mean | 0.319 mean |
| null_mean | 0.428 | 0.4286 mean | 0.428 mean |
| **null_std** | **0.00793** | **0.00800 mean** | **0.0105 mean** |
| ε | 0.234 | 0.263 mean | 0.256 mean |
| **\|z\|** | **12.56** | **13.84 mean** | **10.39 mean** |

**Key observations**:
1. Gap-shuffle null_std (0.0080) is essentially identical to Riemann's (0.0079) — within 1%.
2. RvM null_std (0.0105) is 31% wider than Riemann's.
3. **Gap-shuffle |z| (13.84) is HIGHER than Riemann |z| (12.56)** by ~10%.
4. **Gap-shuffle ε (0.263) is HIGHER than Riemann ε (0.234)** by ~12%.

## K_shuffle

K_shuffle = (|z_Riemann| − mean(|z_shuffle|)) / std(|z_shuffle|)
         = (12.56 − 13.84) / 0.287
         = **−4.44**

The K is **negative**: Riemann |z| is 4.44σ *below* the gap-shuffle distribution. Gap-shuffle produces *stronger* detection than actual Riemann zeros.

## Pre-registered verdict

Verdict criteria (from PRE_REGISTRATION.md):

| Criterion | Met? |
|---|---|
| K_shuffle < 1.5 | **TRUE** (K = −4.44) |
| null_std mean within 15% of Riemann 0.0079 | **TRUE** (mean = 0.0080, 1% off) |
| K_shuffle > 3 | False |
| null_std within 30% of RvM 0.0105 | True (loose criterion) |

H1 fires: both K_shuffle < 1.5 AND null_std within 15% of Riemann.

**Verdict: NEAREST-NEIGHBOR SUFFICIENT (H1).**

The K = 4.5σ distinction between Riemann and RvM-Poisson (Sub F follow-up) reduces entirely to a 2-point statistic — specifically, the nearest-neighbor gap distribution. Long-range correlations beyond pair are *not* responsible for the null_std tightening.

## Findings beyond verdict — the surprising structural twist

The pre-registered verdict is clean (H1), but the data reveal something more interesting than the pre-reg anticipated. The gap-shuffle gives **better detection than actual Riemann zeros** — both in |z| (13.84 vs 12.56) and in ε (0.263 vs 0.234).

This decomposes the K = 4.5σ Riemann-vs-RvM gap into two distinct components:

| Component | Comparison | Δε | Δ\|z\| |
|---|---|---|---|
| **A** (null_std tightening from GUE spacing) | gap-shuffle vs RvM | +0.007 | +3.45 |
| **B** (signal degradation from long-range correlations) | Riemann vs gap-shuffle | −0.029 | −1.28 |
| **Net** (A − B) | Riemann vs RvM | −0.022 | +2.17 |

(All values comparing means; signs from the perspective of moving Riemann → other target.)

Component A is what Mr Code's Sub F follow-up reading captured ("GUE-spacing tightens null"). It is entirely a 2-point property — gap-shuffle preserves the nearest-neighbor structure and reproduces the null_std tightening.

Component B is NEW from Sub G: **actual Riemann zeros are HARDER to detect than gap-shuffled targets.** The Riemann zeros' long-range correlations *reduce* probe-target alignment. Gap-shuffling, by destroying these correlations, IMPROVES detection.

The net K = 4.5σ Riemann-vs-RvM advantage is dominated by Component A (null tightening) and partially offset by Component B (signal degradation). The Riemann |z| is bigger than RvM, but smaller than gap-shuffle's.

### Mechanistic reading (the new picture)

The probe's minima have a distribution along the t-axis that:
1. Has approximately GUE-like local spacing (level repulsion) and similar one-point density to the L-zero distribution.
2. Lacks the LONG-range correlations of L-zeros.

When matched against:
- **Random Poisson points**: weak detection (looser null due to no spacing structure; probe minima don't share Poisson clustering).
- **Gap-shuffled L-zeros**: STRONGEST detection (probe minima share local spacing with target; no long-range conflict).
- **Actual L-zeros**: intermediate (same local spacing as gap-shuffle, but long-range correlations of L-zeros that the probe's minima distribution doesn't share → ~10% degradation in alignment).

So the probe minima themselves behave more like a GUE-random-eigenvalue spectrum (locally repulsive, no long-range L-zero arithmetic structure) than like actual Riemann zeros. **The probe-vs-target correlation is best when the target also lacks long-range L-zero correlations.**

## Implications for Paper 150 v2.3

(Pending CinC adjudication. The brief specifies v2.3 lands from this verdict, so I lay out what the v2.3 mechanism section should now say.)

The honest mechanism summary:

> **The probe (golden-angle Dirichlet sum with indefinite norm) has minima that distribute along the critical line with approximately GUE-like local spacing. The detection of any target on the critical line decomposes into two effects:**
>
> **(A) Null variance tightening**: when the target itself has GUE-like nearest-neighbor spacing, the random-minima MC null distribution is tighter. This produces higher |z| for L-zero targets vs Poisson targets. This effect is a 2-point statistic (Sub G confirmed: gap-shuffle preserves it).
>
> **(B) Signal alignment**: when the target's long-range correlations match the probe's minima distribution, alignment is best. For actual L-zeros, long-range correlations are DIFFERENT from the probe's minima distribution, REDUCING alignment compared to gap-shuffled L-zeros.
>
> **The K = 4.5σ Riemann vs RvM-Poisson distinction is dominated by effect (A); effect (B) partially offsets it but does not erase it. The detection is therefore L-zero-related but in a particular technical sense: the L-zeros' nearest-neighbor structure (shared with the probe) drives the null tightening, while their longer-range structure (NOT shared with the probe) degrades signal alignment.**

This is a substantially different mechanism story from any of:
- v2.0: "probe detects L-zeros specifically" — wrong (gap-shuffle is detected better than L-zeros).
- v2.1: "dense intersection" — incomplete (doesn't predict null_std tightening or signal degradation).
- v2.2: "probe detects density not zeros" — Sub F's reading, refuted by Sub F follow-up.
- Sub F follow-up: "GUE-spacing tightens null" — partially correct (matches Sub G H1) but missed signal degradation.

The cleanest v2.3 reading is what Sub G empirically establishes:
- The probe is a **2-point-spacing detector**, not an L-zero detector.
- Its null tightening is sensitive to nearest-neighbor structure (any GUE-like spacing).
- Its signal alignment is reduced by long-range correlations specific to L-zeros (vs gap-shuffled GUE-spaced targets).

Net effect: the probe's correlation with L-zeros is real but not L-zero-specific in either direction — it's BETTER for some non-L targets (gap-shuffle) and WORSE for some (RvM-Poisson), with both effects on the same null/signal scale.

## Honest pre-vs-post

My pre-registered prior: weak LONGER-RANGE MATTERS (H2). I expected longer-range correlations to be needed for the null_std tightening — specifically anticipating null_std_shuffle to be intermediate between 0.0079 (Riemann) and 0.0105 (RvM), around 0.009.

Reality: null_std_shuffle = 0.0080, essentially identical to Riemann's. **My prior was wrong — fourth in a row**. The null_std tightening is fully captured by the nearest-neighbor distribution, exactly as CinC predicted.

But the SECOND finding (gap-shuffle outperforming Riemann in detection) was anticipated by neither of us. My Sub F follow-up reading had the GUE-spacing-mediated story but didn't predict that gap-shuffle would produce *higher* |z| than Riemann. This is the productive surprise of Sub G.

## What I will NOT do (per brief: no further empirical iteration)

- Run further mechanism tests (H2 vs H3 distinction would need GUE-eigenvalue targets; out of scope).
- Reframe the verdict thresholds.
- Add seeds.

Per brief: "Whatever Sub G returns, v2.3 lands from there. No further empirical iteration on Paper 150 after this."

## Honest limitations

- Only 10 shuffles. A larger ensemble (50-100) would tighten std(|z_shuffle|) = 0.287 estimate, possibly affecting K with greater precision. But the qualitative result (gap-shuffle |z| > Riemann |z| by ~1.3 with std 0.29) is robust.
- Component B (signal degradation from long-range correlations) is a NEW empirical observation; I do not have a theoretical explanation for it beyond noting it. The probe's own minima distribution must lack the specific long-range structure of L-zeros; what *kind* of long-range structure the probe's minima have is open work.
- The H2 vs H3 distinction (longer-range correlations vs L-zero arithmetic specificity) was acknowledged in the brief as not directly addressable by Sub G. A future test against GUE-random-eigenvalue targets would distinguish H2 from H3. Out of scope here per brief.

## Pattern flags

- **Pattern 19 (adversary)**: Mr A's STRUCTURAL catch ("what about pair correlation?") fully addressed. Sub G isolates pair correlation as the source of the K = 4.5σ effect. The cleaner mechanism story for v2.3 is more honest than any previous version.
- **Pattern 75 (null)**: gap-shuffle null is itself a more sophisticated null than RvM-Poisson, designed to isolate the pair-correlation contribution. The pre-registration discipline is now operating on second-order nulls.
- **Pattern 39 (DERIVED vs OBSERVED)**: gap-shuffle preservation properties (first/last point, gap multiset) are DERIVED. The 10-seed |z| distribution, the K = −4.44, and the gap-shuffle-beats-Riemann finding are all OBSERVED.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
- [sub_g.py](sub_g.py)
- [sub_g_target_shuffle_results.csv](sub_g_target_shuffle_results.csv)
- [findings_paper_150_sub_g_target_shuffle.md](findings_paper_150_sub_g_target_shuffle.md) — this report

## Compute

- Probe Re/Im (one-shot): ~33 s.
- 10 gap-shuffle seeds (match + 1000-trial MC null each): ~18 s each, ~3 min total.
- Total: ~4 min. Inside brief's 5-min estimate.
