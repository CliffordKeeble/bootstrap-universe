# Pre-Registration — Test 1b: Parallel-Convention Density Test

**Frozen at commit time. Never edited. Use git's commit hash as the
timestamp authority.**

Per Brief 4 May 2026 (CinC), Mr Adversary follow-up to Paper 197 v0.1.
Independent test of the second heuristic convention; neither a re-run nor
a continuation of Test 1.

---

## Hypothesis

**FRAMEWORK-IMPLIED, NOT DERIVED.** The digital boundaries at 10⁶¹ and 10⁸⁰
register as density-transition signatures in the joint counting function

> Δ(T) = N_ζ(T) − N_{L(χ₅)}(T)

at heights determined by the **resonance-period convention**:

> γ_d = 2π · d / ln(10)

| Boundary | d  | Predicted γ                          |
|---|---|---|
| 1        | 61 | 2π · 61 / ln(10) ≈ **166.4524**     |
| 2        | 80 | 2π · 80 / ln(10) ≈ **218.2982**     |

This is the second of two defensible conventions linking digit-gap scale
10ᵈ to zero-height test windows. Test 1 (Paper 197 v0.1) tested the first
convention γ = d · ln(10) and returned null; this test asks whether the
second convention shows the predicted signature. **Neither convention is
rigorously derived from the explicit formula** — both are heuristics with
reasonable intuitions behind them.

If both conventions return null, the density-transition framing is
strengthened-null and the framework's empirical anchor moves decisively to
the special-value framing. If this convention returns positive, family-wise
α correction (Bonferroni-style; with two conventions, per-test threshold
should be 0.025 to maintain α = 0.05 family-wise) applies in v1.0 review.

---

## Locked parameters

| Parameter | Value | Source |
|---|---|---|
| Boundary 1 prediction (γ_61)        | 166.4524                          | 2π · 61 / ln(10) |
| Boundary 2 prediction (γ_80)        | 218.2982                          | 2π · 80 / ln(10) |
| Boundary 1 window                   | `[161.4524, 171.4524]`            | ±5 about γ_61 |
| Boundary 2 window                   | `[213.2982, 223.2982]`            | ±5 about γ_80 |
| Window half-width                   | 5                                 | matches Test 1 |
| Test statistic                      | RMS of R(T) within window         | matches Test 1 |
| R(T) definition                     | Δ(T) − Δ_smooth(T)                | matches Test 1 |
| Δ_smooth(T) (Selberg null)          | −(T/2π) · log(5) + 1              | DERIVED, audited |
| T-grid step                         | 0.05                              | matches Test 1 |
| Permutation null                    | 1000 random windows of width 10   | matches Test 1 |
| Window-fits-inside-control-piece    | mandatory (no straddling)         | matches Test 1 |
| RNG seed                            | **316** (= 2 · 137 + 42)          | brief seed scheme, trial 2 |
| p-value convention                  | one-sided, (count + 1)/(N + 1)    | matches Test 1 |
| **Pass criterion**                  | p_emp(1) < 0.05 AND p_emp(2) < 0.05 | conjunctive |
| **Stop-on-fail**                    | yes; do not run further tests     | matches Test 1 |

### Control regions

Carved to avoid the four boundary windows (Test 1 + Test 1b) and the
no-coverage region above the L(χ₅) sweep (γ_max = 244.36). Each window
must fit entirely inside one piece — no straddling.

| Region | Interval        | Length | Notes |
|---|---|---|---|
| C1 | [55, 130]       | 75     | between origin and Test 1 boundary 1 |
| C2 | [150, 161.4524] | 11.45  | between Test 1 boundary 1 and Test 1b boundary 1 — fits one width-10 window |
| C3 | [189, 213.2982] | 24.30  | between Test 1 boundary 2 and Test 1b boundary 2 — fits two non-overlapping width-10 windows |
| C4 | [223.2982, 239.3594] | 16.06 | above Test 1b boundary 2; upper bound = γ_max − 5 = 244.3594 − 5; truncation rule did not trigger (γ_max ≥ 233) |

Pre-registered fallback rule, locked at calibration time:
- If γ_max < 233, drop C4 entirely.
- Otherwise, C4 = [223.2982, γ_max − 5].

After Stage 1 coverage extension: γ_max = 244.3594 → C4 = [223.2982, 239.3594].

---

## Pre-registered outcome interpretations

The interpretation of any result is fixed before the result is seen.

### Outcome A — Both boundaries fail
Conditions: p_emp(1) ≥ 0.05 AND p_emp(2) ≥ 0.05.

The density-transition framing is null at *both* tested conventions
(Test 1 and Test 1b). Paper 197 v1.0 reports both nulls; the framework's
empirical anchor moves decisively to the special-value framing in
[Paper 198]. This is the strengthened-null outcome.

### Outcome B — Both boundaries pass
Conditions: p_emp(1) < 0.05 AND p_emp(2) < 0.05.

The resonance-period convention shows the predicted signature; the
framework has empirical support for the density-transition framing — at
this convention, not the natural-log one. Paper 197's title and structure
change radically; this becomes a positive-result paper.

**Mr Adversary review of v1.0 must scrutinise the family-wise
false-positive question:** with two conventions tested, the per-test
threshold should be 0.025 (not 0.05) to maintain family-wise α = 0.05.
Whether the observed p-values clear that bar is part of the v1.0 story.

### Outcome C — Asymmetric (one passes, one fails)
The most ambiguous outcome. Paper 197 v1.0 reports honestly and notes
the asymmetry. Possible readings:

(a) The convention works for the first boundary but the second is at the
edge of coverage and may suffer from edge effects.

(b) The two boundaries actually have different physical mechanisms at
play, consistent with the digit-gaps paper's distinction between
staggered (10⁶¹) and coincident (10⁸⁰) gaps.

Either reading deserves its own follow-up. **No third-convention test
launched in response to this outcome (Pattern 41 — parking discipline).**

### Outcome D — Boundary 2 marginal band
Conditions: any boundary lands in 0.05 ≤ p < 0.15.

Treat as null per pre-registration. Note the marginal status honestly in
the paper. Do not run further tests.

---

## Status flag summary

- **DERIVED.** The Selberg null formula Δ_smooth(T) = −(T/2π) log(5) + 1
  (calibrated through T=230 in Stage 1's 9-point audit, drift indicator
  −0.0006, all bounds satisfied).
- **FRAMEWORK-IMPLIED.** The resonance-period convention γ = 2π·d/ln(10).
  Defensible heuristic; not rigorously derived from the explicit formula.
- **OBSERVED.** Whatever Test 1b returns. Reported plainly per Pattern 28.

---

## Audit gap

This file lives in the **pre-registration commit** of the three-commit
chain for Test 1b:

1. coverage extension + 9-point calibration audit (commit `08aee65`)
2. **pre-registration commit (this file + locked test script)** — the
   parameters above are fixed at this commit's timestamp; the test has
   not yet run; the result does not yet exist
3. result commit (Stage 3) — full numerical output, README, plot

The audit gap is the time between commit (2) and commit (3). The
pre-registration discipline is what makes the audit gap real. Any post-hoc
modification of the locked parameters above invalidates the test.

🐕☕⬡

*Frozen pre-registration. Mr Code, 4 May 2026.*
