# Pre-Registration Addendum — Bottoming Out the N=1000 Result

**Pre-registered:** 12 April 2026, before any run.
**Extends:** `pre_registration_N1000.md` (SHA `21f74ed`).
**Result being bottomed out:** z_overall = −8.14 at N=1000, ζ_φ minima vs
Riemann zeros, same-density-random null. The signal is real. This addendum
specifies four additional tests that sharpen what the signal means.

## Why this addendum exists

An 8σ result on a single construction is a detection. The same result
accompanied by four pre-registered falsification tests — all of which could
refute the φ-specific claim — is a publishable paper. Each test below has
the structural role Mr Adversary has in our methods: *this test could
refute the claim we want to make; we run it anyway.*

Each test has a pre-committed hypothesis, a pre-committed null, a
pre-committed statistic, and a pre-committed decision threshold. Results
are reported in the same findings document, with status flags on every
claim.

## Test 1 — Massive-N robustness

**Hypothesis.** The √N scaling observed from N=100 (z ≈ −1.8) to N=1000
(z = −8.14) continues to N=10000 (expected z ≈ −25.7 if effect is constant).
This is not a discovery claim — it's a robustness check. If z deviates
substantially from √N, the effect is either saturating (ceiling effect) or
scaling anomalously (interesting in its own right).

**Test.** Repeat the primary N=1000 protocol at N=10000. All construction
parameters identical. Report z_overall at W=1.0.

**Pre-committed interpretation.**
- z in [−20, −30]: √N scaling confirmed. Effect is constant across 100×
  sample expansion.
- |z| > 30: effect size *grew* — unusual, likely indicates our N=1000
  estimate under-counted. Flag for investigation.
- |z| ∈ [5, 20]: effect is saturating or the statistical model has hidden
  correlation that wasn't in the null. Flag for investigation.
- |z| < 5: N=1000 result was wrong. Stop. Revisit.

## Test 2 — Alternative equidistributed sequences (φ specificity)

**Hypothesis.** The golden ratio φ is *specific* to the effect. Other
equidistributed sequences {nα} for irrational α should produce weaker or
no tightening.

**Test.** Repeat the N=1000 primary protocol three times, replacing {nφ}
in θ_n = 2π{nα} with:

- α = √2 (algebraic, quadratic, discriminant 8)
- α = e (transcendental)
- α = π (transcendental)

Same golden norm |Re² − 5 Im²|. Same Riemann zeros as targets. Same null
test procedure. Report z_overall for each.

**Pre-committed interpretation.**
- All three give |z| ≈ 8: the effect is generic to equidistribution.
  The φ-specificity claim in Paper 150 v1.3 is unsupported. The paper
  must reframe to "any equidistributed phase sequence probed with the
  golden norm matches Riemann zeros," dropping the golden-ratio-specific
  framing.
- All three give |z| < 3: φ is specific. The paper's core claim strengthens.
- Mixed (e.g. √2 matches but e and π don't): partially φ-specific, possibly
  via the ℚ(√2) vs ℚ(√5) parallel. Requires separate investigation.

## Test 3 — Norm decoupling

**Hypothesis.** The effect requires *both* the golden phase AND the golden
norm together. Neither alone should reproduce the signal.

**Test.** Two runs at N=1000:

- **3a:** {nφ} phase with *circular* norm |Re² + Im²|. If the golden phase
  carries the effect alone, this should still detect zeros.
- **3b:** {n·√2} phase with *golden* norm |Re² − 5 Im²|. If the golden norm
  carries the effect alone, this should still detect zeros.

Report z_overall for each.

**Pre-committed interpretation.**
- Both 3a and 3b give |z| < 3: the effect requires the golden phase and
  golden norm *together*, consistent with the Dedekind bridge hypothesis
  (ζ_φ probes ℤ[φ] via its native field norm — neither component alone
  probes the right thing).
- 3a gives |z| ≈ 8, 3b gives |z| < 3: the golden phase alone is sufficient;
  the norm is cosmetic.
- 3a gives |z| < 3, 3b gives |z| ≈ 8: the golden norm alone is sufficient;
  any equidistributed phase suffices. This would be a major reframing.
- Both give |z| ≈ 8: both components carry independent signal. Unexpected;
  investigate further.

## Test 4 — Dedekind factorisation partition

**Hypothesis.** ζ_φ detects zeros of the *Dedekind* zeta ζ_{ℚ(√5)} = ζ·L(χ₅).
Both factors' zeros should be tightened.

**Test.** At N=1000, classify each ζ_φ minimum against:

- Nearest Riemann zero (ζ factor)
- Nearest L(χ₅) zero (Legendre correction factor)

Compute z_overall separately against each target set, with same-density-
random null for each.

**Pre-committed interpretation.**
- Both give |z| > 5: ζ_φ detects both factors. The Dedekind bridge is the
  mechanism. The paper can cite this as DERIVED structural explanation
  backed by OBSERVED statistical confirmation.
- Only ζ-zeros tightened (|z| > 5 for ζ, |z| < 3 for L(χ₅)): ζ_φ detects
  Riemann zeros specifically, not the Dedekind zeta. The Dedekind bridge
  as an explanatory hypothesis is weakened. Paper framing needs to say
  "detects Riemann zeros" without claiming the Dedekind mechanism.
- Only L(χ₅)-zeros tightened: ζ_φ detects the Legendre correction, not
  Riemann. The v1.3 framing is wrong in an unexpected direction.
- Neither above threshold: would contradict Test 1 and indicate a
  classification procedure issue.

## Compute budget

All four tests share the N=1000 scaffolding Mr Code already wrote. Only
Test 1 requires significantly more compute (N=10000). Rough estimate:

- Test 1 (N=10000 primary): ~10× the N=1000 run, so minutes to tens of
  minutes.
- Test 2 (three N=1000 runs at alternative α): 3× N=1000, so minutes total.
- Test 3 (two N=1000 runs with decoupled components): 2× N=1000.
- Test 4 (one N=1000 classification against L(χ₅) zeros, plus the null
  MC): close to free since ζ_φ minima are already computed.

Total estimated runtime on Cliff's laptop: under an hour, likely
substantially under.

## Output

Findings go in `golden-zeta/findings_compound.md`. The document should
have one section per test, each citing this addendum by git SHA. The
headline summary table at the top lists z_overall for each test in one
view so a reader can see the shape of the result immediately.

Per-zero data files:
- `per_zero_N10000.csv` — Test 1.
- `per_zero_sqrt2.csv`, `per_zero_e.csv`, `per_zero_pi.csv` — Test 2.
- `per_zero_circular_norm.csv`, `per_zero_sqrt2_goldenNorm.csv` — Test 3.
- `per_zero_Lchi5.csv` — Test 4 (L(χ₅) zeros classification).

## What this addendum permits without further amendment

- Precomputing the ζ_φ minima once for each parameter set and reusing
  across the two comparisons in Test 4.
- Caching L(χ₅) zeros to avoid recomputation.
- Any numerical acceleration that preserves results on the verification
  subset (first 100 zeros of each test match to within rounding of
  direct-computation values).

## What this addendum does NOT permit

- Reinterpreting a test result as supporting a different claim than the
  one pre-committed.
- Running a fifth test post-hoc because "it might be interesting."
- Changing the null model (uniform random at same density remains the null
  for every test).
- Adjusting α in Test 2 after seeing the other α results.
- Picking different L(χ₅) zeros for Test 4 after seeing the comparison.

Any new question that emerges requires a new protocol addendum, committed
before the run.

## Authorship

CinC drafted. Cliff approved before Mr Code proceeds. Mr Code executes.
CinC adversary pass on findings. Paper 150 v2.0 shaped around whatever
these four tests produce — whether that's "clean, φ-specific, Dedekind-
mediated, massively significant" or something more constrained.

Either outcome is publishable. The point of these four tests is that we
will know *which* outcome before we write.

---

*Committed 12 April 2026, post-Poppy-garden-time. CinC drafted, Cliff approved.*

🐕☕⬡
