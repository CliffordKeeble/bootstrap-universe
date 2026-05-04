# WINDUP — riemann_digit_gap_test arc terminus

*4 May 2026. Final note on this arc, before merge.*

This arc tested the **density-transition framing** of the framework's
digital-boundaries claim: that χ₅-feedback structure visible in the digit
gaps at decimal scales 10⁶¹ and 10⁸⁰ should also be visible as anomalies
in the joint ζ / L(χ₅) zero distribution at corresponding heights γ ≈ 140
and γ ≈ 184.

**Result: null at registered thresholds.**

Test 1 (the cleanest, density-transition): boundary 1 sat at the null mean
(p = 0.353); boundary 2 was marginal but did not pass (p = 0.092).
Pre-registered pass criterion `p₁ < 0.05 AND p₂ < 0.05` was not met.
Calibration audit on the Selberg null had passed before the test ran, so
the negative result is informative — not an artefact of null
mis-specification.

Per pre-registered discipline, Tests 2 and 3 were not run. No rescue
commits, no statistic-shopping, no window-narrowing.

## What survives

- **DERIVED** — Dedekind factorisation `ζ_{ℚ(√5)}(s) = ζ(s) · L(s, χ₅)`,
  the Selberg null formula, the calibration audit at non-boundary heights.
- **OBSERVED** — the spectral floor (Paper 158), Pattern 56, the digit-gaps
  observation itself (Paper v3.9). Each on its own evidence.
- **The unification claim** in *this specific framing* (digit-gap boundaries
  visible as zero-density fingerprint) is not confirmed.

## Where the empirical anchor moves

Per follow-up brief: **special-value transitions**, not density transitions.
Different test, different statistics, fresh pre-registration. This arc
ends here; the next arc starts elsewhere with a clean slate.

## Pattern citations

- **Pattern 74** — work is only real when merged; this includes negative
  results. The branch carries a clean record of what was tried and what
  was ruled out.
- **Pattern 75** — numerical match is a claim about separation from null,
  and this one separated less than chance. Discipline held under live
  pressure.
- **Pattern 28** — show, don't tell. The plot, the JSON, and the README
  speak for themselves; this WINDUP.md just marks the closing gate.
- **Pattern 14** — the negative result is itself a lemma. The framework
  architecture updates (anchor moves to special-value transitions); the
  individual programme components stand.

---

## Addendum — Test 1b (parallel-convention) result, 4 May 2026

Mr Adversary follow-up to Paper 197 v0.1 raised: only one of two
defensible heuristic conventions for linking digit-gap scale 10ᵈ to
zero-height γ had been tested. Test 1b ran the second convention
γ = 2π · d / ln(10) under independent pre-registration (separate seed,
separate windows, separate control regions, separate folder).

**Result: FAIL. Outcome A — both boundaries null at decisive p.**

| | window | observed RMS | p_emp |
|---|---|---|---|
| Boundary 1 (γ ≈ 166.45) | [161.45, 171.45] | 0.6510 | 0.429 |
| Boundary 2 (γ ≈ 218.30) | [213.30, 223.30] | 0.6438 | 0.459 |

Neither marginal. The signal that *had* been suggestive in Test 1
boundary 2 (p = 0.092) disappears entirely under the resonance-period
convention. The two conventions disagree about predicted γ by factor
~1.19, and neither location shows windowed RMS anomaly above permutation
null.

The density-transition framing is now null at **both** defensible
conventions. The original arc message stands and is strengthened: the
framework's empirical anchor moves decisively to special-value
transitions for the next paper. No third-convention test launched
(Pattern 41 — parking discipline).

Calibration upgraded to 9 points (Mr Adversary recommendation),
calibrated cleanly through T = 230 (mean drift −0.0006).

Sub-arc files: `test1b_parallel_convention/` (pre-registration,
calibration audit, locked test, result, plot, README).

🐕☕⬡

— Mr Code, signing off this arc (now with parallel-convention closure).
