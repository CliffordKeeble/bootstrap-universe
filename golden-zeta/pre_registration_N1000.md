# Pre-Registration Protocol — Golden Zeta Extension to N=1000

**Pre-registered:** 12 April 2026, before any extended run.
**Purpose:** establish whether ζ_φ's systematic improvement over same-density
random minima, observed at N=100 (mean Δ ~30% tighter than null, z ≈ −1.8
overall), survives extension and at what sigma level.

This protocol is committed to the repository *before* the extension is
computed. Any deviation from the protocol must be recorded in a follow-up
commit and justified explicitly in the findings. The point of pre-registration
is that the result we report is an honest test, not a post-hoc selection.

## Hypothesis

**H₀ (null).** The minima of ζ_φ(t) match Riemann zero locations no better
than same-density uniform random minima would, once classification window
width is controlled.

**H₁ (signal).** ζ_φ's minima match Riemann zero locations systematically
tighter than same-density random minima, with the effect robust across the
range tested.

## Data

- **Riemann zeros:** first 1000 non-trivial Riemann zeros from mpmath's
  `zetazero(n)` or equivalent LMFDB-sourced high-precision values. Range
  approximately t ∈ [14.13, 1419.4].
- **ζ_φ construction:** identical to Paper 150 v1.3 §2.1–2.2, N_terms = 5000
  (primary), with N_terms = 10000 as robustness check. Fejér smoothing.
  Golden-angle phase θ_n = 2π{nφ}. Golden norm √|Re² − 5·Im²|.
- **Scan:** Δt = 0.008 over the full range [1, 1420]. Minima identified by
  local three-point comparison below the 5th percentile threshold.
  Deduplication window 0.3.

Reproducibility: all RNG seeds explicit. Python version, numpy/mpmath
versions recorded in findings.

## Test statistics — pre-registered

For each Riemann zero t_n, define Δ_n = distance to nearest ζ_φ minimum
within window W.

### Primary statistic (pre-committed)

**z_overall = (mean(Δ_signal) − mean(Δ_null)) / SE(Δ_null)**

where mean(Δ_null) and SE(Δ_null) are estimated from 1000 Monte Carlo trials
of same-density random minima matched to the same 1000 zeros. W = 1.0
(matching the v1.3 classification window; note that the null itself uses
this window).

### Secondary statistics (also pre-committed, reported alongside)

1. **Bin-resolved z-scores.** Partition zeros into four t-bins:
   [14, 300], [300, 650], [650, 1000], [1000, 1420]. Compute z separately.
   Report all four. If the effect is genuinely uniform, bins should be
   consistent to within their individual SEs. If the effect is concentrated
   in one bin, this is informative.
2. **Window sensitivity.** Recompute z_overall at W = 0.5, 0.3, 0.2.
   Report all four. If the z-score is robust to window tightening, the
   signal is window-independent. If it collapses, the signal is
   window-coupled.
3. **Match rate at W = 0.3.** Fraction of zeros with a ζ_φ minimum within
   0.3 of them, compared to the same rate under null. Report both.
4. **N_terms sensitivity.** Rerun at N_terms = 10000. Report whether z_overall
   shifts materially.

## Decision thresholds — pre-committed

| z_overall | Interpretation | Paper 150 v2.0 framing |
|---|---|---|
| \|z\| < 2 | Not significant | "Modest improvement over chance, within 2σ of null" |
| 2 ≤ \|z\| < 3 | Suggestive | "Systematic improvement over chance at 2σ; not a detection" |
| 3 ≤ \|z\| < 5 | Evidence | "3σ evidence for systematic ζ_φ–zero correspondence" |
| \|z\| ≥ 5 | Discovery | "5σ discovery of systematic ζ_φ–zero correspondence" |

**Direction matters.** A positive z (signal mean Δ *larger* than null) would
be anti-evidence — ζ_φ fits zeros worse than random. We report whichever
direction the z takes, not only the favourable one.

**Uniformity matters.** If z_overall is 5σ but one of the four bin z-scores
has the opposite sign, the "discovery" framing is unavailable until the
inhomogeneity is understood. In that case, the v2.0 framing becomes:
"Aggregate significance 5σ, but non-uniform across t — concentrated in
bins X." This is a real result and reported as such, not dressed up as
a universal law.

## What would count as refutation

- z_overall ∈ [−1, +1] at N = 1000 with the same ~−1.8 seen at N = 100 would
  indicate the earlier signal was statistical fluctuation — the effect does
  not scale as √N and is not present.
- Window-collapse: if z_overall is −3 at W = 1.0 but approaches 0 at W = 0.3,
  the result is a window artefact, not a signal. Paper 150 v2.0 would then
  remove the ζ_φ-matching claim entirely and retain only the Dedekind/Zeckendorf
  results.

## Non-deviations permitted without protocol amendment

- Replacing mpmath with a faster backend if numerical agreement is verified
  on the first 100 zeros to high precision.
- Adjusting Δt downward (finer scan) if coarse grid misses obvious minima.
- Extending the Monte Carlo null to 10000 trials if 1000 gives noisy SE.

## Deviations that would require protocol amendment (committed before run)

- Changing the construction of ζ_φ (different weights, different phase).
- Changing the classification procedure for minima.
- Changing the null model (e.g., from uniform to clustered random).
- Changing the definition of the test statistic.

If any of these need to change, the amendment is committed first, in a
separate commit, before running the extended analysis.

## Reporting

Findings go in `golden-zeta/findings_N1000.md`. The narrative there cites
this protocol file by git SHA (so the pre-registration is verifiable post
hoc). Key numbers reproduced in the commit message for the findings commit.

Per-zero CSV `per_zero_1000.csv` captures all 1000 zeros with columns:
n, t_n, nearest_minimum_t, Δ_n, nearest_null_minimum_t, Δ_n_null, bin.

## Authorship and review

Drafted by CinC. Reviewed by Cliff before run. Executed by Mr Code.
Adversary pass by CinC on findings before any paper incorporation.
No shortcut through that loop.

---

*Committed 12 April 2026. CinC drafted, Cliff approved.*

🐕☕⬡
