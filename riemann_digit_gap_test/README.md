# riemann_digit_gap_test

Spike to test the framework prediction that χ₅-feedback structure visible in
the digit-gaps boundaries (decimal scales 10⁶¹ and 10⁸⁰) should also be
visible in the existing Riemann ζ and L(s, χ₅) zero data at corresponding
heights γ ≈ 140 and γ ≈ 184.

Per Mr Code brief, 4 May 2026 (CinC). Three tests in increasing order of
structural ambition; stop on first fail.

---

## Result — Test 1: FAIL

**The framework prediction at γ = 140 is not borne out by zero-count data.
γ = 184 is marginal. Tests 2 and 3 not run, per pre-registered stop-on-fail
discipline.**

| | observed RMS | p-value | verdict |
|---|---|---|---|
| Boundary 1 [135, 145] (γ≈140, x=10⁶¹) | 0.643 | **0.353** | fail |
| Boundary 2 [179, 189] (γ≈184, x=10⁸⁰) | 0.689 | **0.092** | fail (marginal) |

Pass criterion was `p₁ < 0.05 AND p₂ < 0.05`; both fail. Boundary 1 sits at
the null mean (no detectable signal). Boundary 2 sits ~1.3σ above the null
mean — within CinC's pre-flagged "inconclusive band" (0.05–0.15) but not
significant at the registered threshold.

Calibration audit on Selberg null (run before the test, committed in
`zeros_data.py` self-test): residuals at non-boundary heights T = 60, 80,
100, 120, 160 sit at +0.37, +0.49, −0.39, −0.26, −0.02. Mean +0.04, range
[−0.39, +0.49], all well within ±1.5. Linear-plus-constant null is correctly
calibrated, so the windowed boundary test was meaningful.

## What this means

DERIVED: the Selberg null formula `Δ̄(T) = −(T/2π) log(5) + 1`. Constants
derived from the standard Riemann–von Mangoldt formula for ζ (+7/8) and the
even-character analogue for L(s, χ₅) (conductor 5, parity δ_a = −1/8).

OBSERVED: the empirical residuals at γ ≈ 140 and γ ≈ 184 are *not*
distinguishable from residuals at non-boundary control heights, given the
fluctuation scale set by S(T).

FRAMEWORK-IMPLICATION (FAIL): the prediction that digit-gap boundaries
correspond to *visible* structural transitions in the joint ζ / L(χ₅) zero
distribution is not confirmed by the simplest test — counting zeros up to
height T and looking for excess fluctuation in a pre-registered window.

## What survives

The negative result here does not invalidate:

- The Dedekind factorisation `ζ_{ℚ(√5)}(s) = ζ(s) · L(s, χ₅)` (DERIVED).
- The spectral floor (Paper 158, separate evidence).
- The Pattern 56 framework (separate evidence).
- The digit-gaps observation itself (Paper v3.9, separate evidence).

Each of those stands or falls on its own merits. **The unification claim —
that the digit-gap boundaries are a direct numerical fingerprint of the
χ₅-feedback structure visible in the zero distribution — weakens.**

## What might be tried, with appropriate caution

If CinC and Cliff judge the marginal boundary-2 result worth pursuing
(p = 0.092 is suggestive but not significant), candidate follow-ups:

1. **Extend the L(χ₅) sweep** beyond γ = 220 to admit control windows up to
   T = 250 (would marginally increase test power; CinC said upfront this is
   only worthwhile if Test 1 lands in the marginal band, which boundary 2
   does).
2. **Try a different test statistic** that uses the *signed* residual or
   the local autocorrelation rather than RMS. But: Pattern 75 says we don't
   tune until we get a signal. Any statistic-shopping after seeing Test 1's
   data is post-hoc and the result loses meaning.
3. **Investigate boundary 2 separately** with a finer-grained statistic at
   that specific height, but only if pre-registered as a *new* test — i.e.
   a fresh brief with a fresh seed and fresh windows.

The honest interpretation right now is: the simplest framework-implied
test does not see the predicted signal at γ ≈ 140 and is borderline at
γ ≈ 184. If the missing paper (provisional Paper 197) exists, it will need
either a different empirical anchor or a sharpened test.

## Files

- `zeros_data.py` — shared zero loader, Selberg null formulas, calibration
  audit. Self-test prints the audit table and constants on each run.
- `zeros_cache/` — cached zero CSVs for reproducible runs.
  - `riemann_zeros_dps30_n100.csv` (100 zeta zeros, mpmath.zetazero, dps=30)
  - `lchi5_zeros_dps30_n150.csv` (146 L(χ₅) zeros, max γ ≈ 219.25)
- `test1_density.py` — Test 1, pre-registered constants at top of file.
- `test1_density_results.json` — full numerics + verdict.
- `test1_density.png` — 4-panel diagnostic (Δ vs Selberg null; residual
  R(T); permutation distribution vs observed; survival function log-y).

Tests 2 and 3 not run. Their scripts (`test2_partial_sum.py`,
`test3_correlation.py`) are not present in this commit.

## Commit chain

```
1. feat(rdgt): zero-data infrastructure + Selberg null + calibration audit
2. feat(rdgt): test1_density scaffolding (pre-registered, not yet run)
3. result(rdgt): Test 1 FAIL - this commit
```

Audit happened in the gap between commits 2 and 3: the pre-registration was
locked before the result existed.

🐕☕⬡

*Pattern 28: show, don't tell. Pattern 75: numerical match is a claim
about separation from null, and this one separated less than chance.
Pattern 14: the negative result is itself a lemma — the unification claim
weakens, and that's a real finding for the missing-paper architecture.*
