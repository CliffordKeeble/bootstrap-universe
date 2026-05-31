# Mr Code brief — Sub F follow-up: RvM-random draw distribution

**Date**: 31 May 2026
**Brief author**: CinC
**Scope**: One small bounded follow-up to Sub F. Addresses Mr Adversary's NULL catch on Paper 150 v2.2: the RvM-random result (z = −9.89) is a single draw of a random object; we need the distribution across multiple seeds to convert "≈ same effect size as Riemann" from eyeball comparison to a quantitative claim. Outcome feeds Paper 150 v2.3 §6.

**Reads required**:
- This brief
- Sub F pre-registration (`88fc4d3`) and findings (`d03517b`)
- Mr Adversary review of Paper 150 v2.2 (NULL catch specifically)

**Reads NOT required** (Pattern 97 — scope protection):
- Paper 150 v2.2 draft itself — this sub doesn't depend on the prose
- Other programme papers — independent scope

---

## Question

The Sub F RvM-random control used a single i.i.d. realisation (seed = 20260530, z = −9.89). What is the *distribution* of z across multiple independent RvM-random realisations, holding the probe, the null methodology, the t-range, and N_target = 2792 fixed?

## Background — Mr Adversary's NULL catch on v2.2

> "Sub F's decisive number, the RvM-random z = −9.89, is one realisation of an i.i.d. random target. A density-matched random target is a random object; one draw gives one z. Did Mr Code draw it once, or many times and report a distribution? If once, the honest statement is 'a single density-matched draw gave z = −9.89'; if many, give me the mean and spread, because that is what turns '≈ same effect size as Riemann' from an eyeball into a result."

Mr A is right. The Sub F single-draw result is one data point. To claim quantitatively that "the probe gives essentially the same effect size against RvM-random as against actual Riemann zeros," we need the distribution of z over multiple RvM draws and a comparison of the Riemann z to that distribution.

This is the cleanest way to land Paper 150 v2.2 (→ v2.3) at ★★★★★ from Mr A's current ★★★★ review. The catch is precise; the fix is small.

## Test design

**Identical to Sub F except for the seed**:

| Parameter | Value |
|---|---|
| Probe | `sqrt(|Re²(Z_φ(t)) − 5·Im²(Z_φ(t))|)`, N_terms = 5000, α = φ (same as Sub F) |
| Grid | t ∈ [1, T_max + 5], dt = 0.008 (same as Sub F) |
| Minima | 3-point, 5th percentile, dedup 0.3 (same as Sub F) |
| Match window | W = 1.0 (same as Sub F) |
| MC null | 1000 trials, seed = 42 (same as Sub F) |
| N_target | 2792 (same as Sub F) |
| mpmath precision | dps = 15 (same as Sub F) |
| Riemann target | first 2792 zeros (cached from Sub F) |
| RvM-random seeds | **10 independent seeds**: 20260530 (Sub F original), 20260531, 20260601, 20260602, 20260603, 20260604, 20260605, 20260606, 20260607, 20260608 (sequential dates from Sub F + 9 follow-ups) |

The probe Re/Im is computed once and reused across all 10 seeds (since the probe doesn't depend on the target). The matching step and effect-size computation runs once per seed.

For each seed, record:
- z_RvM(seed)
- ε_RvM(seed)
- Number of matched minima (sanity check; should be ~3304 as in Sub F)

## Pre-registered analysis

Before running:

1. Compute the empirical distribution of z over the 10 RvM-random draws: mean(z), std(z), min(z), max(z), median(z).
2. Compute the empirical distribution of ε similarly.
3. Compare the Riemann control z (= −12.56 per Sub F) to the RvM-random distribution. Report:
   - K-σ separation: K = (|z_Riemann| − mean(|z_RvM|)) / std(|z_RvM|)
   - Within / outside the RvM-random range
   - Percentile of |z_Riemann| relative to the RvM-random distribution

## Pre-registered verdicts

Three pre-registered readings of the comparison:

- **CONFIRM SUSTAIN** (K < 1, i.e., Riemann z within 1σ of RvM-random mean): the probe truly does not distinguish L-zeros from random density-matched points. Sub F's verdict stands; "essentially the same effect size" is empirically established as a quantitative claim. Paper 150 v2.2 (→ v2.3) retraction holds at full strength.

- **PARTIAL DISTINCTION** (1 ≤ K < 3): the Riemann control sits 1-3σ above the RvM-random distribution. Some weak L-zero-specific signal survives Sub F's null, though it is much smaller than v2.0 claimed. Paper 150 v2.3 needs to acknowledge this nuance: the retraction holds but the framing should say "the probe gives an effect size against L-zeros that is at most marginally distinguishable from arbitrary density-matched targets at K-σ."

- **REJECT SUSTAIN** (K ≥ 3): the Riemann control is significantly outside the RvM-random distribution. The original Sub F verdict was a Type-I error driven by a particularly low single seed; the probe does distinguish L-zeros from random points after all. Paper 150 v2.3 must rewrite the retraction substantively. CinC adjudication required.

## Stop-on-fail

If **REJECT SUSTAIN** (K ≥ 3): pause Paper 150 v2.3 drafting. The retraction may need to be reframed or partially withdrawn. Flag for CinC adjudication.

If **PARTIAL DISTINCTION** (1 ≤ K < 3): proceed with v2.3 drafting but adjust the abstract and §6/§7/§12 to acknowledge the partial-distinction nuance.

If **CONFIRM SUSTAIN** (K < 1): proceed with v2.3 drafting; the retraction claim is now quantitatively grounded.

CinC's honest prior, recorded before the experiment: I expect CONFIRM SUSTAIN. The Sub F single-draw result (z = −9.89) is close enough to the Riemann control (z = −12.56) that 9 additional draws should produce a distribution centred near z ≈ −10 to −11 with std ~1-2, putting Riemann at K ≈ 1-2σ. But I have been wrong about my priors twice in a row (Sub E and Sub F), so the empirical answer is the only route.

## Anti-circularity

- The probe and Riemann cache are unchanged from Sub F. No reanalysis of Riemann.
- The RvM-density sampling uses inverse-CDF on log(t)/(2π), same as Sub F's seed-20260530 draw.
- Seeds are pre-registered (10 specific integers); no post-hoc seed selection.
- All 10 seed results are reported, regardless of outcome. No cherry-picking.

## Implementation notes

- Reuse `sub_f.py` machinery with seed parameterised.
- Loop over the 10 seeds, capturing z and ε for each.
- Output a CSV (`sub_f_followup_rvm_distribution.csv`) with one row per seed.
- ~30 seconds per seed (matching step + MC null). Total: ~5 minutes.
- The 10 seeds include the original Sub F seed (20260530) so the result is self-consistent with the Sub F single-draw finding.

## Pre-registration

Commit to repository:

- BRIEF.md (this file)
- PRE_REGISTRATION.md (the 10 seeds, the three pre-registered verdicts, the K-σ threshold)
- sub_f_followup.py (script skeleton, BEFORE results)

Tag the commit. Then run.

## Deliverable

Markdown report `findings_paper_150_sub_f_followup_rvm_distribution.md` with:

- Methodology (delta from Sub F: just the seeds)
- Per-seed table (z, ε, n_matched)
- Empirical distribution statistics (mean, std, range, median for both z and ε)
- Comparison to Riemann control: K-σ separation, percentile
- Verdict per pre-registered thresholds
- Brief note on implications for Paper 150 v2.3

## Compute budget

~5 minutes total. Bounded.

## Order

Single bounded sub-task. Run. Report. No parallelism.

---

⌨️ over to you.

Small follow-up that converts Sub F's central admissible-control number from "an eyeball comparison" to "a quantitative claim about distribution overlap." The result either confirms the retraction at its full strength (CONFIRM SUSTAIN), nuances it (PARTIAL DISTINCTION), or forces a rewrite (REJECT SUSTAIN). My prior leans toward CONFIRM but the empirical answer is what matters.

🐕☕⬡

---

## Mr Code's report

*[To be filled in.]*
