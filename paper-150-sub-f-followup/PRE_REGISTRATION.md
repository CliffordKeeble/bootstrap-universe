# Pre-registration — Sub F follow-up (RvM-random draw distribution)

**Date**: 31 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 31 May 2026)
**Binding**: this document is committed to git **before any new seed is run**.

---

## Honest expectation, recorded BEFORE the experiment

My pre-registered prior leans **CONFIRM SUSTAIN**. Reasoning: the single
draw z = −9.89 vs Riemann z = −12.56 is a gap of ~2.7 in |z|. If the
RvM distribution has std(|z|) ~ 1.5 (a rough guess based on typical
single-seed sampling noise in this setup), the gap is K ≈ 1.8σ — i.e.,
PARTIAL DISTINCTION territory.

But my priors have been wrong twice in a row (Sub E Flag 1, Sub F
expectation of RvM COLLAPSE), so I am not confident. The empirical
answer is what matters.

## Construction commitments

Identical to Sub F except for the seed:

- Probe: `sqrt(|Re² − 5·Im²|)`, N_terms = 5000, α = φ.
- Probe grid, dt, t-padding, minima detection, match window W = 1.0,
  MC null (1000 trials, seed = 42 internal), N_target = 2792,
  mpmath dps = 15: ALL identical to Sub F.
- The probe (Re, Im) is computed **once** and reused across the 10
  RvM-random seeds.

## The 10 RvM-random seeds (frozen at pre-registration)

```
[20260530, 20260531, 20260601, 20260602, 20260603,
 20260604, 20260605, 20260606, 20260607, 20260608]
```

Seed 20260530 is the Sub F original; the other 9 are sequential.

## Pre-registered analysis

- Mean, std, min, max, median of z over the 10 draws.
- Mean, std, min, max, median of ε over the 10 draws.
- K-σ comparison: `K = (|z_Riemann| − mean(|z_RvM|)) / std(|z_RvM|)`
  where z_Riemann = −12.56 from Sub F's cached Riemann control.
- Percentile of |z_Riemann| relative to the 10-draw distribution.

## Pre-registered verdicts (re-stated)

- **CONFIRM SUSTAIN**: K < 1 — Riemann within 1σ of RvM mean.
- **PARTIAL DISTINCTION**: 1 ≤ K < 3 — Riemann 1-3σ above RvM.
- **REJECT SUSTAIN**: K ≥ 3 — Riemann ≥ 3σ above RvM.

## What I will NOT do

- Add or remove seeds after the run starts.
- Cherry-pick which seeds to report.
- Reframe verdict thresholds.
- Rerun a "bad" seed if the verdict is unfavourable.

## Files

- `sub_f_followup.py` — implementation
- `sub_f_followup_rvm_distribution.csv` — per-seed results
- `findings_paper_150_sub_f_followup_rvm_distribution.md` — long-form report
