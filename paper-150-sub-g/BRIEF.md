# Mr Code brief — Sub G: target-shuffle null (gap-shuffle test)

**Date**: 31 May 2026
**Brief author**: CinC
**Scope**: One bounded follow-up to Sub F follow-up. Tests Mr Code's GUE-spacing-mediated hypothesis directly via gap-shuffle of the Riemann zero target. Outcome decides Paper 150 v2.3's mechanism section. **This is explicitly the LAST mechanism-clarifying empirical step before v2.3** — if it gives a clean result either way, v2.3 is drafted from there.

**Reads required**:
- This brief
- Sub F follow-up pre-registration (`078cc38`) and findings (`cd2ab8a`)
- Mr Adversary review of Paper 150 v2.2 (NULL catch and STRUCTURAL catch)
- Sub F findings (`d03517b`) for null_std numbers

**Reads NOT required** (Pattern 97 — scope protection):
- Paper 150 v2.2 / v2.3 drafts — this sub doesn't depend on the prose
- Sub C, Sub D, Sub E — independent scope
- Other programme papers — independent scope

---

## Question

Sub F follow-up established K = 4.45σ distinguishing the probe's |z| against Riemann zeros from its |z| against RvM-Poisson density-matched random points. Mr Code interpreted this as GUE-spacing-mediated: the Riemann zeros' GUE-like pair correlation tightens the MC null variance (Riemann null_std = 0.0079 vs RvM null_std = 0.0105, 25% gap), raising |z| relative to Poisson targets.

**Is this GUE-spacing-mediated interpretation empirically supported?** Specifically: does the nearest-neighbor gap distribution of Riemann zeros (the dominant 2-point structure) account for the K = 4.5σ distinction, or is something beyond nearest-neighbor (long-range correlations, n-point correlations beyond pair, L-function-specific structure) also at play?

## Background

Three hypotheses are in play:

(H1) **Nearest-neighbor sufficient**: K = 4.5σ comes entirely from the GUE-like nearest-neighbor gap distribution. Any target with the same nearest-neighbor gap distribution and same overall density would give the same null_std.

(H2) **Longer-range correlations matter**: K = 4.5σ depends on n-point correlations beyond nearest-neighbor (number variance, three-point functions). Nearest-neighbor alone is insufficient.

(H3) **L-zero specific beyond pair correlation**: There's something about the actual Riemann zero positions — beyond their statistical pair-correlation structure — that drives the null_std distinction. (E.g., arithmetic content of the zeros, not captured by pair correlation statistics.)

Sub G distinguishes these via a gap-shuffle target.

## Test design

**Same probe, same null methodology, same N_target = 2792, same Riemann control**:
- Probe: `sqrt(|Re²(Z_φ(t)) − 5·Im²(Z_φ(t))|)`, N_terms = 5000, α = φ
- Grid: t ∈ [1, T_max + 5], dt = 0.008 (Sub F follow-up's grid)
- Minima: 3-point, 5th percentile, dedup 0.3, W = 1.0
- MC null: 1000 trials, seed = 42
- mpmath dps = 15
- All other parameters frozen

**What changes**: a new target type — **gap-shuffled Riemann zeros**.

**Gap-shuffle procedure**:
1. Take the 2792 Riemann zero positions t_1 < t_2 < ... < t_{2792} (cached from Sub F).
2. Compute the 2791 gaps g_i = t_{i+1} - t_i.
3. Randomly permute the gap order: g_{σ(1)}, g_{σ(2)}, ..., g_{σ(2791)} for permutation σ.
4. Reconstruct shuffled positions: t'_1 = t_1, t'_{i+1} = t'_i + g_{σ(i)}.
5. Result: 2792 points covering [t_1, t_{2792}] (same start, same end since sum-of-gaps preserved), with the same multiset of nearest-neighbor gaps, but randomized order.

**Properties of gap-shuffled Riemann**:
- ✓ Same overall density (same count in same range)
- ✓ Same nearest-neighbor gap distribution (multiset of gaps preserved)
- ✗ Long-range correlations destroyed (gap order randomized)
- ✗ Three-point and higher-n correlations destroyed
- ✗ Any L-function-specific structure beyond gap distribution destroyed

**10 independent gap-shuffles** with pre-registered seeds:
```
[20260530, 20260531, 20260601, 20260602, 20260603,
 20260604, 20260605, 20260606, 20260607, 20260608]
```
(Same seeds as Sub F follow-up for cross-comparison; each seed parameterises a numpy permutation.)

For each shuffle:
- Compute z, ε, null_mean, null_std vs the gap-shuffled target.

Plus the existing Riemann control from Sub F (z = −12.56) for comparison.

## Pre-registered analysis

Across the 10 gap-shuffle realizations:
- mean(z_shuffle), std(z_shuffle), min, max
- mean(ε_shuffle), std(ε_shuffle)
- mean(null_std_shuffle), std(null_std_shuffle)

Compare to Sub F numbers:
- Riemann: z = −12.56, null_std = 0.0079
- RvM-random (Sub F follow-up mean): z = −10.39, null_std = 0.0105

Define K_shuffle = (|z_Riemann| − mean(|z_shuffle|)) / std(|z_shuffle|).

## Pre-registered verdicts

Three pre-registered readings:

- **NEAREST-NEIGHBOR SUFFICIENT (H1 confirmed)**: K_shuffle < 1.5 AND mean(null_std_shuffle) is within 15% of Riemann null_std = 0.0079. The gap-shuffle reproduces Riemann's null_std. **Conclusion: K = 4.5σ is entirely about the nearest-neighbor gap distribution; nothing beyond pair correlation matters. Paper 150 v2.3's mechanism section: the K-σ distinction is a generic property of GUE-like spacing distributions, not L-function-specific.**

- **LONGER-RANGE MATTERS (H2 confirmed)**: K_shuffle > 3 AND mean(null_std_shuffle) is within 30% of RvM null_std = 0.0105. The gap-shuffle does NOT reproduce Riemann's null_std; it reverts to roughly RvM-like behavior despite preserving nearest-neighbor gaps. **Conclusion: K = 4.5σ depends on longer-range correlations beyond nearest-neighbor; nearest-neighbor distribution alone is not the driver. Paper 150 v2.3: mechanism is GUE-correlations-beyond-pair, narrower than "any GUE-spaced target."**

- **PARTIAL (H1 + H2 mixed)**: 1.5 ≤ K_shuffle ≤ 3, OR mean(null_std_shuffle) is intermediate. Both nearest-neighbor and longer-range correlations contribute. **Conclusion: nuanced mechanism; v2.3 acknowledges both components.**

- **H3 NOT TESTABLE BY THIS SUB**: distinguishing H2 from H3 requires comparing GUE-random-eigenvalue targets to gap-shuffled Riemann, which is a separate investigation. Sub G can refute H1 (in favour of H2 or H3) but cannot distinguish H2 from H3 alone.

## Stop-on-fail (or rather, stop-on-clean-verdict)

This is explicitly the **last** mechanism-clarifying empirical step before v2.3:
- If verdict is clean (NEAREST-NEIGHBOR SUFFICIENT or LONGER-RANGE MATTERS), proceed to v2.3 drafting from the corresponding mechanism reading.
- If verdict is PARTIAL, proceed to v2.3 drafting acknowledging the nuance; do NOT request further empirical investigation. We commit to landing v2.3 from whatever Sub G returns.

This protects against "just one more test" iteration creep. The empirical work for Paper 150 v2.3 is bounded by Sub G's conclusion.

## Anti-circularity

- Probe, Riemann cache, MC null methodology, N_target — all unchanged from Sub F follow-up. No re-analysis.
- Gap-shuffle procedure is target-side only.
- 10 seeds pre-registered (same as Sub F follow-up).
- All 10 results reported regardless of outcome.
- The gap distribution preserved by shuffle is the EMPIRICAL distribution of Riemann gaps (not a fit to a theoretical GUE). This is exact in the sense that the shuffle preserves what's actually there.

## Implementation notes

- Reuse `sub_f_followup.py` machinery with new target type.
- Gap-shuffle: load cached Riemann zeros, compute gaps via np.diff, np.random.permutation with seed, cumulative sum starting from t_1.
- Per-seed: ~30 sec (same as Sub F follow-up per-seed cost).
- Total: ~5 min.

## Honest CinC prior (recorded BEFORE the experiment)

I weakly lean **NEAREST-NEIGHBOR SUFFICIENT (H1)**. The reason: the dominant effect of GUE-like statistics on null variance should come from level repulsion at short scales, which IS captured by the nearest-neighbor gap distribution. Longer-range correlations (number variance scaling as log N for GUE vs N for Poisson) might contribute, but my intuition is that the K = 4.5σ effect is mostly a 2-point statistic.

But I have been wrong about Sub E, Sub F, and Sub F follow-up priors. And Mr Code has been wrong about all three of his priors. The empirical answer is the only route.

## Pre-registration

Commit to repository:
- BRIEF.md (this file)
- PRE_REGISTRATION.md (the 10 seeds, the K_shuffle thresholds, the null_std thresholds)
- sub_g.py (skeleton with gap-shuffle implementation, BEFORE results)

Tag the commit. Then run.

## Deliverable

Markdown report `findings_paper_150_sub_g_target_shuffle.md` with:
- Methodology delta from Sub F follow-up (just the gap-shuffle target)
- Per-seed table (z, ε, null_std, n_matched)
- Aggregate statistics across 10 shuffles
- Comparison: Riemann (Sub F cache) vs gap-shuffle (Sub G) vs RvM-Poisson (Sub F follow-up)
- K_shuffle and null_std comparison
- Verdict per pre-registered thresholds
- Brief mechanism implications for Paper 150 v2.3

## Compute budget

~5 minutes total. Bounded.

## Order

Single bounded sub-task. Run. Report. The empirical work for Paper 150 v2.3 ends here.

---

⌨️ over to you.

The third bounded test in the Paper 150 chain (after Sub F and its follow-up). The mechanism story has progressively refined from "golden norm hears Dedekind zeta" → "indefinite norm hears L-zeros generically" → "probe detects density not zeros" → "probe distinguishes targets via null variance, mediated by pair correlation." Sub G tests whether that pair correlation reduces to nearest-neighbor structure.

Whatever Sub G returns, v2.3 lands from there. No further empirical iteration on Paper 150 after this.

🐕☕⬡

---

## Mr Code's report

**Verdict: NEAREST-NEIGHBOR SUFFICIENT (H1)** — with a clean, productive
surprise (see "Findings beyond verdict" below).

Long form: [findings_paper_150_sub_g_target_shuffle.md](findings_paper_150_sub_g_target_shuffle.md).
Pre-registration at `a3a04a0`.

### Per-shuffle distribution (10 seeds)

| Statistic | z | ε | null_std |
|---|---:|---:|---:|
| mean | **−13.84** | 0.263 | **0.00800** |
| std | 0.287 | 0.005 | 0.00012 |
| min | −14.25 | 0.255 | 0.00793 |
| max | −13.21 | 0.273 | 0.00829 |

### Three-way comparison

| Quantity | Riemann (Sub F) | Gap-shuffle (Sub G) | RvM-Poisson (Sub F follow-up) |
|---|---:|---:|---:|
| **null_std** | **0.0079** | **0.0080** | **0.0105** |
| ε | 0.234 | 0.263 | 0.256 |
| **\|z\|** | **12.56** | **13.84** | **10.39** |

### K_shuffle

K_shuffle = (12.56 − 13.84) / 0.287 = **−4.44**.

**Riemann |z| is 4.44σ BELOW the gap-shuffle distribution mean.**
0 of 10 shuffles give |z| ≤ Riemann.

### Pre-registered verdict criteria

- K_shuffle < 1.5: **TRUE** (K = −4.44).
- null_std mean within 15% of Riemann's 0.0079: **TRUE** (0.0080, 1% off).
- → **H1 NEAREST-NEIGHBOR SUFFICIENT** confirmed.

### Findings beyond verdict — the productive surprise

The K = 4.5σ Riemann-vs-RvM distinction decomposes into:

- **Component A (null tightening)**: gap-shuffle preserves it. Riemann
  nearest-neighbor structure is enough. Pure 2-point statistic.

- **Component B (signal degradation)**: NEW from Sub G. **Actual Riemann
  zeros give weaker probe-target alignment than gap-shuffled targets.**
  Riemann ε = 0.234; gap-shuffle ε mean = 0.263; gap-shuffle |z| = 13.84
  vs Riemann |z| = 12.56. **The Riemann long-range correlations REDUCE
  detection by ~10%** in both ε and |z|.

The mechanism story for Paper 150 v2.3:
- Probe minima distribute on the critical line with approximately
  GUE-like local spacing (level repulsion at short scales).
- Null tightening for L-zero targets: 2-point property (H1).
- Signal alignment is BEST for gap-shuffled targets (same local spacing,
  no long-range structure conflicting with probe).
- The actual Riemann |z| being above RvM but below gap-shuffle reflects
  competing effects: GUE-spacing tightens null (advantage) but
  long-range correlations degrade alignment (disadvantage). The net is
  positive (K = 4.5σ Riemann vs RvM) but smaller than gap-shuffle's
  K vs RvM (about 6σ).

### Pre-vs-post

My weak prior: LONGER-RANGE MATTERS (H2). Reality: H1, null fully
captured by nearest-neighbor. **Fourth consecutive wrong prior** (Sub E,
Sub F, Sub F follow-up, Sub G). CinC's prior (H1) was correct.

The productive surprise (gap-shuffle outperforming Riemann) was
unanticipated by either of us. Sub G earns its keep beyond just
distinguishing H1 from H2 — it surfaces Component B, which sharpens the
v2.3 mechanism reading substantially.

### v2.3 directive (per brief)

> "Whatever Sub G returns, v2.3 lands from there. No further empirical
> iteration on Paper 150 after this."

Per brief: I am ready to assist with v2.3 drafting. The mechanism
section should now describe:
1. Probe = 2-point-spacing detector with GUE-like minima distribution.
2. Detection has TWO components: null tightening (2-point of target) +
   signal alignment (long-range structure of target).
3. L-zeros benefit from (1) but suffer slight degradation in (2);
   net effect is K = 4.5σ Riemann-vs-RvM.
4. Gap-shuffled L-zeros maximise BOTH components → strongest detection
   yet observed in the chain.

The Sub F-line retraction stands at full strength in **effect-size**
terms: Riemann ε < RvM ε < gap-shuffle ε. The probe does not
"specifically detect L-zeros" by alignment strength.
In **z-score** terms, Riemann is distinguished from RvM by null
tightening (4.5σ), but is DOMINATED by gap-shuffle (4.44σ in the other
direction). So z-space "L-zero specificity" is also refuted at the
"gap-shuffle vs Riemann" level.

The mechanism is **GUE-spacing-with-Riemann-long-range-anti-correlation**.
That's the cleanest, most honest statement.

### Compute

~4 min total, inside brief's 5-min estimate.
