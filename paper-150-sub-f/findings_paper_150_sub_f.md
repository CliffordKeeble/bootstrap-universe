# Sub F findings — density-matched non-L target null

**Date**: 31 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `88fc4d3` BEFORE any non-Riemann target was generated and BEFORE the probe was run against any non-L set.
**Verdict**: **SUSTAIN** (aggregate). RvM-random alone gives PARTIAL just below SUSTAIN; Gram triggers SUSTAIN but is interleaving-uninformative per the anti-circularity check. Either reading triggers stop-on-fail.

Per brief: **pause Paper 150 v2.1 publication. Do not upload to Zenodo. Flag for CinC adjudication.**

---

## Methodology

Identical to Sub C-1 except for the target set:

| Parameter | Value |
|---|---|
| Probe | `sqrt(|Re²(Z_φ(t)) − 5·Im²(Z_φ(t))|)`, N_terms = 5000, α = φ |
| Grid | t ∈ [1, T_max + 5], dt = 0.008 |
| Minima | 3-point, 5th percentile, dedup 0.3 |
| Match window | W = 1.0 |
| MC null | 1000 trials, seed = 42 |
| N_target | 2792 (matches Sub C-ext-1 for direct comparison) |
| mpmath precision | dps = 15 |

The probe `(Re, Im)` was computed *once*; the same minima set was matched against each of the three targets. Only the *target set* changes between runs.

## Targets

1. **Riemann zeros** (control): first 2792 non-trivial Riemann zeros via `mpmath.zetazero(k)`. T_max = t-coordinate of 2792nd Riemann zero ≈ 3326.55.

2. **Gram points**: g_n satisfying `siegeltheta(g_n) = (n−1)·π`. 2792 computed by `mpmath.findroot` with local-bracket initialisation using the Gram-spacing asymptotic `2π / log(t/(2π))`. Range [17.85, 3326.58].

3. **RvM-density-matched random**: 2792 i.i.d. samples from `ρ(t) = log(t)/(2π)` on (1, T_max + 5] via inverse-CDF (bisection). NumPy RNG seed = 20260530 (per pre-registration).

Range matching of all three targets is confirmed: all sit in the same t-interval at the same density.

## Anti-circularity: Gram–Riemann distance distribution

Per pre-registration the Gram-zero distance was a mandatory check before interpreting the Gram result.

| Statistic | Value |
|---|---|
| Median | 0.378 |
| 10th percentile | 0.092 |
| 90th percentile | 0.794 |
| Within W = 1.0 | 2681 / 2792 = **96.0%** |
| Within W/2 = 0.5 | 1845 / 2792 = **66.1%** |

**The Gram-point control is NOT informative.** Per pre-registration: "if median is below W = 1.0, then Gram detection ≈ Riemann detection by interleaving... the RvM-random control becomes the decisive test." Median = 0.378 < 1.0, and 96% of Gram points sit within W of a Riemann zero. The Gram probe-vs-target matching is effectively the same as the Riemann matching.

**Decisive test: RvM-random.**

## Detection results

| Target | n_minima | signal_mean | null_mean | null_std | **z** | **ε** |
|---|---|---|---|---|---:|---:|
| Riemann (control) | 3304 | 0.328 | 0.428 | 0.008 | **−12.56** | **0.234** |
| Gram | 3304 | 0.320 | 0.428 | 0.008 | **−14.38** | **0.253** |
| RvM random | 3304 | 0.324 | 0.430 | 0.011 | **−9.89** | **0.246** |

Effect sizes are **indistinguishable** across the three targets:
- Riemann ε = 0.234
- Gram ε = 0.253
- RvM ε = 0.246

The probe's effect size on RvM-random points (ε = 0.246) is essentially the same as on actual Riemann zeros (ε = 0.234) — the difference (~5%) is well within sampling noise.

## Pre-registered verdicts per target

| Target | \|z\| | Threshold | Verdict |
|---|---:|---|---|
| Gram | 14.38 | ≥ 10 → SUSTAIN | **SUSTAIN** (but interleaving-uninformative) |
| RvM random | 9.89 | 3 ≤ \|z\| < 10 → PARTIAL | **PARTIAL** (just below SUSTAIN) |

### Aggregate verdict per pre-reg conjunction rule

> "SUSTAIN on either → SUSTAIN verdict (one positive is enough to refute L-zero specificity)"

**Aggregate verdict: SUSTAIN.**

### Aggregate verdict adjusted for anti-circularity

Strict reading of the conjunction rule: SUSTAIN (Gram triggers).
Anti-circularity-adjusted reading (Gram inadmissible because interleaved with Riemann at sub-W scale): **PARTIAL** (RvM at |z| = 9.89, just below the SUSTAIN threshold of 10).

**Either reading triggers pause-and-flag.** The strict reading is SUSTAIN; the adjusted reading is "very strong PARTIAL approaching SUSTAIN".

## What the data actually show

The probe `|Re² − 5·Im²|` of the Fejér-weighted golden-angle Dirichlet sum **detects ANY sufficiently-dense critical-line-positioned point set** at effect size ε ≈ 0.25:

- On real Riemann zeros: ε = 0.234.
- On Gram points (interleaved with Riemann): ε = 0.253.
- On i.i.d. random points sampled from the RvM density: ε = 0.246.

The matching effect size on **uniformly random points** with no L-function content refutes the v2.0 and v2.1 mechanism stories at the empirical level. **Whatever the probe is detecting, it is not L-function-zero structure specifically; it is the density of the target on the critical line.**

Mr Adversary's gating catch is **fully sustained**:
> "Two Poisson processes of comparable rate are uncorrelated and would give z ≈ 0, not z = −22.91."

Indeed: two genuinely *independent* Poisson processes would give z ≈ 0. The probe's minima set, however, is **not independent of the target set's density** — it tracks the t-axis density of *whatever* is there. The mechanism that makes this happen is not L-function-specific; it is something about the probe's minima distribution that correlates with any density-matched point set.

## Implications for Paper 150 v2.1

Per brief stop-on-fail:
> "Pause Paper 150 v2.1 publication. Do not upload to Zenodo. Flag for CinC. Do NOT attempt to 'rescue' v2.1's mechanism story by post-hoc rationalisation."

I **pause** v2.1 drafting and publication. v2.1's §7.2 dense-intersection mechanism is empirically refuted by RvM-random at ε ≈ 0.246. The L-function-specific framing of v2.0 is doubly refuted.

What the result *positively* establishes (without claiming a mechanism):
- The probe minima form a point set on (1, T] at density ~2.5/unit.
- This minima set has a fixed positive correlation with ANY point set of comparable density on the same interval.
- The correlation is **ε ≈ 0.25** — 25% tightening of mean distance vs. uniform random.

The next paper would need to ask: **what is it about the probe's minima distribution that produces a ε ≈ 0.25 correlation with arbitrary same-density targets?** Mr A's framing: "Two Poisson processes of comparable rate are uncorrelated and would give z ≈ 0, not z = −22.91." So the probe minima are *not* Poisson-distributed; they have some structure that makes them auto-correlate with density-matched targets.

This is genuinely puzzling and would deserve a separate investigation. Sub F establishes the negative result; the positive mechanism question is open.

## Honest limitations

- N = 2792 was the pre-registered N. Larger N might give slightly different z values; the trend would not change.
- The RvM-random target uses a single seed (20260530). A handful of additional seeds would test seed-stability of the RvM result.
- The GUE-pair-correlation-matched target (optional Target 4 per brief) was not run. It would test whether GUE-style higher-correlation matters. Given that even fully-i.i.d. RvM-random sustains detection, GUE would almost certainly also sustain.
- The pipe `compute_probe → find_minima → MC null` is the Sub C-1 standard machinery. If there's a methodological issue with this pipeline, all of Sub C and Sub F would inherit it.

## Pattern flags

- **Pattern 19 (adversary)**: Mr Adversary's gating catch fully sustained. The detection is density-driven, not L-zero-driven. Stop-on-fail triggered.
- **Pattern 75 (null)**: same null methodology as Sub C-1. The RvM-random target IS a null — and the probe "detects" it. This is the null-collapses-the-claim outcome.
- **Pattern 39 (DERIVED vs OBSERVED)**: probe construction is DERIVED. The z-scores and effect sizes are OBSERVED. The mechanism interpretation is now demoted to "open question."

## What I will NOT do per brief stop-on-fail

- Attempt to "rescue" v2.0 or v2.1's mechanism story. The pre-registered SUSTAIN verdict stands.
- Draft a v2.2 reframing. CinC must adjudicate first.
- Reinterpret the verdict to avoid the SUSTAIN label. The data are what they are: the probe detects density.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — pre-registered methodology, seed, thresholds
- [sub_f.py](sub_f.py) — implementation (with Unicode bug-fix after the first run; the original output captures the pre-fix run, which produced identical numerics)
- [sub_f_results.csv](sub_f_results.csv) — raw z, ε per target
- [gram_points.csv](gram_points.csv) — 2792 Gram points (cached for reproducibility)
- (Riemann zeros cached at `paper-203-sub-c/riemann_zeros_10000.csv`)
- [findings_paper_150_sub_f.md](findings_paper_150_sub_f.md) — this report

## Honest expectation vs. result

My pre-registered prior:
> "weak prior leans COLLAPSE on RvM-random and PARTIAL on Gram (because Gram-Riemann interleaving may make the test less discriminating)."

Reality:
- Gram: SUSTAIN (interleaving-uninformative as suspected)
- RvM-random: PARTIAL very near SUSTAIN

My prior was wrong on RvM. I predicted COLLAPSE; reality gave PARTIAL/SUSTAIN. Mr A's catch was correct; my "weak prior" leaned toward L-zero specificity, which is now refuted.

This is the second time in this brief sequence my structural reading has been overturned by direct test (after Sub E's "d is decoration" refutation). The pre-registration discipline keeps working: claims are pre-committed, tests run cleanly, refuted claims stand refuted.

## Compute

- Riemann extension (1900 → 2792): ~30 min (mpmath.zetazero at dps=15, slowing at high t).
- Gram points (2792 via siegeltheta + findroot): ~2 min (asymptotic-bracket init made this fast).
- RvM-random sampling: <1 s.
- Probe Re/Im over 416,706 grid points: ~30 s.
- Three target matchings + MC nulls: ~1 min total.

Total: ~35 min. Inside the brief's 1-2 hour budget.
