# Pre-registration — Sub G (target gap-shuffle null)

**Date**: 31 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 31 May 2026)
**Binding**: committed to git **before any gap-shuffle is performed**.

---

## Honest expectation, recorded BEFORE the experiment

I've been wrong three times in a row about priors (Sub E Flag 1, Sub F
RvM COLLAPSE, Sub F follow-up CONFIRM SUSTAIN). I'm going to record a
prior, but with low confidence:

**Weak prior leans LONGER-RANGE MATTERS (H2)**.

Reasoning: my analytical reading in Sub F follow-up attributed the
K=4.5σ to GUE-spacing-mediated tightening of null_std. But "GUE
spacing" includes both nearest-neighbor repulsion AND longer-range
number-variance suppression. The latter is the more distinctive GUE
feature vs Poisson. Pure nearest-neighbor repulsion (preserved by
gap-shuffle) still allows the gaps to occur in any order — which
destroys the smooth-distribution-along-the-axis property that comes
from long-range correlations. I'd expect that to noticeably affect
the null_std.

Specifically I'd guess null_std_shuffle ≈ 0.009 (between Riemann's
0.0079 and RvM's 0.0105), giving K_shuffle around 2-3. That'd be
PARTIAL.

But CinC leans H1 NEAREST-NEIGHBOR SUFFICIENT. So we differ.

Either way, the empirical answer decides.

## Construction commitments

Identical to Sub F follow-up except the target is gap-shuffled Riemann
(not RvM-Poisson):

- Probe machinery, MC null, all thresholds, seeds: same.
- 10 pre-registered shuffle seeds:
  `[20260530, 20260531, 20260601, 20260602, 20260603,
    20260604, 20260605, 20260606, 20260607, 20260608]`
- Each seed parameterises `numpy.random.default_rng(seed).permutation(2791)`
  on the gap array.

## Gap-shuffle procedure (pre-registered, exact)

```
zeros = load_riemann_first_n(2792)         # cached from Sub F
gaps  = np.diff(zeros)                     # 2791 gaps
rng   = np.random.default_rng(seed)
perm  = rng.permutation(len(gaps))
shuffled_gaps = gaps[perm]
shuffled_zeros = np.concatenate([[zeros[0]], zeros[0] + np.cumsum(shuffled_gaps)])
```

The shuffled targets:
- Have the same first point (zeros[0]).
- Have the same last point (zeros[0] + sum(gaps) = zeros[-1]).
- Have the same multiset of nearest-neighbor gaps as Riemann.
- Have randomised gap ORDER — destroying long-range structure.

## Pre-registered analysis

For each of 10 seeds:
- Compute z, ε, null_mean, null_std, n_matched.

Aggregate:
- mean and std of z, ε, null_std across 10 seeds.
- K_shuffle = (|z_Riemann| − mean(|z_shuffle|)) / std(|z_shuffle|),
  using Sub F's z_Riemann = −12.56.
- mean(null_std_shuffle) vs Riemann null_std = 0.0079 (Sub F).
- mean(null_std_shuffle) vs RvM null_std mean = 0.0105 (Sub F follow-up).

## Pre-registered verdicts (re-stated)

- **NEAREST-NEIGHBOR SUFFICIENT (H1)**:
  K_shuffle < 1.5 AND mean(null_std_shuffle) within 15% of 0.0079
  (i.e., in [0.00672, 0.00909]).
- **LONGER-RANGE MATTERS (H2)**:
  K_shuffle > 3 AND mean(null_std_shuffle) within 30% of 0.0105
  (i.e., in [0.00735, 0.01365]).
- **PARTIAL**: 1.5 ≤ K_shuffle ≤ 3, OR null_std intermediate
  (not satisfying either H1 or H2 numerical criteria).

Note ambiguity: the null_std bands of H1 and H2 overlap
([0.00735, 0.00909]). If mean(null_std_shuffle) falls in that overlap,
the K_shuffle threshold decides the verdict.

## Anti-circularity

- Gap-shuffle is target-side only.
- 10 seeds pre-registered, all reported.
- The shuffle preserves the empirical gap multiset (not a fit to GUE).
- Probe / null machinery / Riemann cache unchanged.

## What I will NOT do

- Re-run a "bad" seed.
- Adjust K_shuffle or null_std thresholds after seeing data.
- Add or remove seeds.
- Request further empirical iteration regardless of outcome. v2.3 lands
  from Sub G's verdict.

## Files

- `sub_g.py` — gap-shuffle implementation, reuses Sub F machinery
- `sub_g_target_shuffle_results.csv` — per-seed results
- `findings_paper_150_sub_g_target_shuffle.md` — long-form report
