# Mr Code Brief: Paper 48 §5 Audit

**Pre-registered Pattern 75 audit on the m_e c² icosahedral relation**
**Brief: v0.1**
**Date: 25 May 2026**
**Author: Dr. Clifford G. Keeble (ORCID: 0009-0003-6828-2155)**
**For: Mr Code computational session**

---

## Pre-registration commitments

This brief is committed to the working directory by git commit BEFORE the computational script is run, per Pattern 75 discipline. No fishing; single run; locked seed; output reported regardless of result. The reproduction script (`paper_48_s5_audit.py`) is run once with the brief's locked specifications. Mr Code does not iterate on tolerance, search space, sample size, or seed after the brief is committed.

## Background and motivation

Paper 48 v4.4 §5 (*The Hydrogen Bond: Polynomial-Expansion Architecture from Icosahedral Geometry*) records a numerical observation at OBSERVED status:

> Within the eV system of units, the numerical value M = m_e c² / eV satisfies the relation 2(α⁻¹ − 1)α⁻²/(V − χ) to 0.0008%.

The "0.0008%" precision figure was carried from v4.1 through v4.2–v4.4 without recomputation. Hand-arithmetic in the v4.2 colophon turned up a discrepancy (computing 2 × 136.036 × 137.036² / 10 gave ~510,824 against the v4.1 claim of 511,003 against measured M = 510,999). CinC's analysis of the §5 relation also notes algebraic equivalence with the §3.2 fossil claim (the (α⁻¹ − 1)/(V − χ) ≈ R∞hc/eV observation at ~154 ppm precision), suggesting the §5 relation should match at the same ~154 ppm precision rather than at 0.0008%.

This audit settles the §5 precision question definitively and runs the Pattern 75 null distribution that the v4.4 audit schedule requires before any future programme paper claims m_e derivation.

## Part A: Arithmetic verification

**Task A1:** Compute the right-hand side of the §5 relation to 12 significant figures.

Specification:
- α⁻¹ = 137.035999084 (CODATA 2018, exact value used for comparison)
- V = 12, χ = 2 (icosahedral primitives), V − χ = 10
- Formula: `RHS = 2 * (alpha_inv - 1) * alpha_inv**2 / (V - chi)`

Report:
- The computed value of RHS to 12 significant figures
- The measured value of M = m_e c² / eV. Use CODATA 2018: m_e c² = 0.51099895000 MeV = 510998.95000 eV
- The absolute residual |M − RHS| in eV
- The relative residual (M − RHS) / M expressed in: parts-per-million (ppm), parts-per-billion (ppb), parts-per-trillion (ppt)
- The relative residual expressed as a percentage to two significant figures

This is straightforward arithmetic; no sampling, no iteration. One number out, no decisions to make.

**Task A2:** Confirm or refute the algebraic equivalence claim.

CinC's analysis: the §5 relation should be algebraically equivalent to the §3.2 fossil observation. Specifically:
- §3.2 observation: E₀_v3 = (α⁻¹ − 1)/(V − χ) matches R∞hc in eV at ~154 ppm
- §5 relation: M = m_e c² / eV ≈ 2(α⁻¹ − 1)α⁻²/(V − χ)
- These should be algebraically equivalent up to the standard relation R∞hc = m_e c² × α² / 2

Specification:
- Compute |E₀_v3 − R∞hc| / R∞hc (the §3.2 fossil residual) where:
  - E₀_v3 = (α⁻¹ − 1)/(V − χ) (a pure number, dimensionless)
  - R∞hc taken as the measured Rydberg energy in eV: R∞hc = 13.605693122994 eV (CODATA 2018)
  
- Compute |M − RHS| / M (the §5 residual from Task A1)
- Report whether these two residuals match to within a few parts per billion (which they should if the algebraic-equivalence claim holds)
- If they differ by orders of magnitude, the algebraic-equivalence claim is wrong and CinC needs to reconsider

If the two residuals agree to within ~10 ppb, the algebraic equivalence holds and §5 should be rewritten with the corrected precision figure (probably ~154 ppm or 0.015%, not 0.0008%).

## Part B: Pattern 75 null distribution

**Task B1:** Measure the search-space density of icosahedral-primitive constructions that match M to 100 ppm tolerance.

Following Paper 127's §2 audit template directly:

### Locked search space specification

**Primitive set:** S = {e, π, √2, √3, √5} ∪ {integers 1, 2, …, 30}. Cardinality |S| = 35. This is the same primitive set used in Paper 127's §2 audit, chosen for fair comparison.

**Expression form:** Each candidate expression has 1 to 4 terms. Each term is a product or quotient of 1 to 3 distinct primitives drawn from S (without repetition within a term), with each primitive raised to an integer exponent in {−5, −4, −3, −2, −1, 1, 2, 3, 4, 5} (exponent zero excluded). Terms are added with sign + or −.

**Target:** M = 510998.95 (the numerical value of m_e c² in eV from CODATA 2018).

**Tolerance:** 100 ppm absolute = 51.1 eV. A candidate expression "matches" the target if its computed value falls within [M − 51.1, M + 51.1] = [510947.85, 511050.05].

**Sample size:** N = 10⁶ uniform-random samples from the expression space, drawn with random seed = 48 (chosen as the paper number; locked here).

**Decision threshold:** Pre-registered match-rate threshold = 1% (10⁻²). If the random match-rate at 100 ppm tolerance is < 1%, the §5 observation lies in a sparse region of the icosahedral-primitive search space and the OBSERVED status carries non-trivial structural weight. If ≥ 1%, the §5 observation is consistent with generic match-density in this primitive space and the OBSERVED status reflects coincidence rather than structural connection.

### Sanity check on the search space

The §5 relation as printed, 2(α⁻¹ − 1)α⁻²/(V − χ), is itself constructible in the locked search space if α⁻¹ is interpreted via the substitution rule α⁻¹ → 137 (the nearest integer in the primitive set) or via an exact icosahedral construction. For the audit, treat α⁻¹ as an external constant and not part of the primitive set — the audit asks "how many constructions in S match M independently of any prior α-knowledge?" The §5 relation itself uses α⁻¹ explicitly, which is acceptable; what the audit measures is the density of *other* constructions matching M in the locked space.

If the audit cannot resolve this — i.e., if including or excluding α⁻¹ from the primitive set changes the result by orders of magnitude — that itself is informative and should be reported.

### Output

Report exactly:
- N samples drawn
- Number of matches within 100 ppm tolerance
- Match fraction (matches / N)
- 95% upper bound on match fraction via rule of three: 3/N if 0 matches; otherwise binomial confidence interval
- Comparison to pre-registered 1% threshold
- Verdict: § 5 OBSERVED status earns non-trivial standing (match-rate ≪ 1%) or does not (≥ 1%)

**Task B2 (optional, if N = 10⁶ is unhelpful):** Tighter audit at 10 ppm tolerance.

If the 100 ppm audit returns a high match-rate (≫ 1%), a tighter audit may discriminate. Repeat Task B1 with tolerance 10 ppm (5.11 eV) and N = 10⁷ (random seed = 48 still). Report the same metrics.

This is the "if Task B1 result is inconclusive" branch. If Task B1 cleanly returns ≪ 1% or ≫ 1%, B2 is not needed.

## Output format

Mr Code produces:
1. `paper_48_s5_audit.py` — reproduction script
2. `paper_48_s5_audit_results.md` — markdown results document with the items above
3. `paper_48_s5_audit_brief_v0_1.md` (this file, committed before execution)

The results document should be self-contained and citable from Paper 48 v4.5 (the version that will absorb the audit outcome).

## Sign-off conditions

The audit is complete when:
- Task A1: numerical residual reported to 12 sig figs
- Task A2: algebraic equivalence claim confirmed or refuted
- Task B1: match-rate reported with 95% CI, decision threshold compared, verdict given
- (B2 if A2 inconclusive)
- Git commit of all three files with pre-registration brief committed BEFORE the script runs

CinC will read the results and absorb them into Paper 48 v4.5. Likely revision: §5 precision figure corrected from "0.0008%" to the actual residual; §5 status retained at OBSERVED with the audit outcome cited; §8 audit schedule updated to mark §5 audit DISCHARGED with Paper 127-style verdict line; the §3/§4 audit (the larger one) advances to next computational priority.

## Notes for Mr Code

- This is Paper 127 §2 audit template patterned to Paper 48 §5. Same shape, different target.
- Pre-registration is mandatory: git commit of brief before script execution. No tolerance/seed/threshold adjustments after the brief is in.
- If Task A1 returns the actual residual (~154 ppm expected), please flag this explicitly so CinC can correct v4.4 §5 quickly. If Task A1 returns the v4.1 figure (~8 ppm = 0.0008%), the algebraic-equivalence claim of CinC's analysis is wrong and needs reconsideration before B1 runs.
- The random seed 48 is the paper number, chosen for traceability.
- Sample size N = 10⁶ matches Paper 127 §2 audit for cross-comparison.
- 100 ppm tolerance is chosen to over-support the audit: even if §5 actually matches at 154 ppm, the 100 ppm tolerance still tests "matches close enough to be interesting." Tighter tolerances (Task B2) provide additional discrimination if needed.

---

🐕☕⬡

*Brief committed before execution. Single run. Locked specifications.*
*Mr Code: when ready, git commit this brief, then run the script.*
