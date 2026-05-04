# test1b_parallel_convention — result summary

*Written after the result. Never edited to retrofit. Plain reporting per
Pattern 28.*

---

## Result — Test 1b: FAIL (Outcome A)

**The density-transition framing is null at the resonance-period
convention as well. Combined with Test 1's null, this is now null at
*both* tested conventions. The strengthened-null outcome.**

| | observed RMS | p-value | verdict |
|---|---|---|---|
| Boundary 1 [161.45, 171.45] (γ ≈ 166.45, x = 10⁶¹) | 0.6510 | **0.4286** | fail (well within null) |
| Boundary 2 [213.30, 223.30] (γ ≈ 218.30, x = 10⁸⁰) | 0.6438 | **0.4585** | fail (below null mean) |
| Null distribution (n=1000) | mean 0.6275, median 0.6372, std 0.0625, max 0.7892 | — | — |

Pass criterion `p₁ < 0.05 AND p₂ < 0.05`: not met. Both observed RMS values
sit slightly above the null *median* but well inside the bulk of the null
distribution. Boundary 2 is at the 46th percentile of the null —
essentially the null mean.

**No marginal flag** for either boundary. Neither lands in the [0.05, 0.15)
inconclusive band. The null is decisive.

## Pre-registered outcome interpretation: **A**

From `pre_registration.md` §"Outcome A": *"The density-transition framing
is null at both tested conventions (Test 1 and Test 1b). Paper 197 v1.0
reports both nulls; the framework's empirical anchor moves decisively to
the special-value framing in [Paper 198]. This is the strengthened-null
outcome."*

The interpretation was fixed before the result was seen; the result selects
which interpretation activates.

## Combined picture across both tests

| Test | Convention | γ at d=61 | γ at d=80 | p₁ | p₂ | Outcome |
|---|---|---|---|---|---|---|
| 1   | γ = d · ln(10)            | 140.48 | 184.21 | 0.353 | 0.092 | FAIL (boundary 2 marginal) |
| **1b** | γ = 2π · d / ln(10)    | 166.45 | 218.30 | **0.429** | **0.459** | **FAIL (outcome A — both decisive)** |

Test 1's boundary 2 had been the only suggestive signal across either
test (p = 0.092, in the inconclusive band); under the resonance-period
convention, that signal disappears entirely (p = 0.459).

The two conventions disagree about the predicted γ at each boundary by
factor ~1.19 — and neither location shows windowed RMS anomaly above the
permutation null. The framework prediction at the level of *zero-density
transition* does not survive either heuristic linkage.

## What this means

DERIVED — survives unchanged:
- Selberg null `Δ̄(T) = −(T/2π) log(5) + 1` (calibrated through T = 230
  in Stage 1's 9-point audit; mean drift −0.0006).
- The Dedekind factorisation `ζ_{ℚ(√5)}(s) = ζ(s) · L(s, χ₅)`.
- The L(χ₅) zero data (integrity-checked at 25/25 vs LMFDB to 4
  decimals).

OBSERVED — null at both conventions:
- Empirical residuals at γ_61 and γ_80 (under either heuristic) are
  indistinguishable from residuals at non-boundary control heights.

FRAMEWORK-IMPLICATION — strengthened-null:
- The unification claim *in the density-transition framing* is null at
  both defensible conventions. Paper 197's empirical anchor moves
  decisively to special-value transitions per the next brief.

## What this does NOT address

- The deeper question of whether *localised γ-windowing* is the right
  way to detect a *localised x-feature* (Mr Adversary's concern (a) in
  v0.1 review). That's a broadband-vs-localised question, properly
  addressed by a cumulative-residual or partial-sum observable. Out of
  scope for this brief.
- Special-value transitions (Paper 198 territory).
- Pair-correlation or any other observable.

Per pre-registered Pattern 41 discipline: **no third-convention test.**
The parallel-convention arc parks here.

## Files in this folder

| File | Role |
|---|---|
| `pre_registration.md` | Frozen pre-registration; locked at commit `015f043` |
| `extend_lchi5_sweep.py` | Stage 1 coverage extension; integrity check 25/25 vs LMFDB |
| `extended_calibration.py` | Stage 1 9-point Selberg null audit; PASS |
| `extended_calibration_results.json` | Stage 1 audit numerics |
| `test1b_density_parallel.py` | Stage 2 locked test script; constants frozen |
| `test1b_density_parallel_results.json` | Stage 3 result |
| `test1b_density_parallel.png` | Stage 3 4-panel diagnostic |
| `README.md` | This file |

## Commit chain

```
(stage 3) result(rdgt-1b): Test 1b FAIL — outcome A, both conventions null
(stage 2) preregister(rdgt-1b): pre_registration.md + locked test script
(stage 1) feat(rdgt-1b): coverage extension + 9-point calibration audit
```

The audit gap was between stages 2 and 3.

## What did not happen between Stage 2 and Stage 3

- No edit to `pre_registration.md`.
- No edit to the locked-constants block.
- No re-run with a tweaked window or seed.
- No statistic-shopping.
- No third-convention spec drafted.

Pattern 75 + 41 + 28 in operation.

🐕☕⬡

*Mr Code, 4 May 2026 — sub-arc closed.*
