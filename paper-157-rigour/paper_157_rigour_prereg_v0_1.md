# Mr Code Brief — Paper 157 Rigour Tests
## Pre-registration v0.1 — the two tests that move 157 from OBSERVED to DERIVED

**Status: PRE-REGISTRATION. To be committed (timestamped, frozen) BEFORE any computation runs.**
**Discipline: Pattern 75 — null and decision thresholds locked before the data is touched. No fishing. Single run per locked spec. Seed fixed.**

**Author:** CinC, for Cliff, 2I Universe Programme
**Date drafted:** 1 June 2026
**Target paper:** Paper 157 (The Self-Reference Tower), currently α⁻¹ = (|A₅|−1)/L(1,χ₅) = 59√5/(2 log φ) = 137.06 at STATUS: OBSERVED
**Goal:** determine whether 157's two unproven steps can be closed, converting OBSERVED → DERIVED at 157's honest 0.02% precision (NOT higher — see scope guard §0).

---

## §0. Scope guard (read first — this is the anti-numerology fence)

This brief does **not** seek higher precision. Paper 157 sits at 0.02% and that is the target ceiling. The programme's position (Pattern 98: the discipline IS the pitch — topological consistency, not super-precision) is that 157's modest precision with fully-forced terms beats Paper 90's parts-per-trillion with soft degrees of freedom. These tests aim to convert 157's two *assumptions* into *consequences* — they do not aim to add decimal places.

If either test is passed, the win is **rigour, not precision**: 157's formula becomes derived rather than observed, at unchanged precision. If a test tempts a precision-improvement sub-result, that sub-result is OUT OF SCOPE and flagged for separate pre-registration. The back-EMF / proton-g−2 precision direction is **pocketed as Phase 2, not activated** (see §4).

---

## §1. What 157 already has (not under test)

These are taken as given; the brief does not re-test them:

- **|A₅| − 1 = 59** — the augmentation ideal dimension of the regular representation of A₅. Forced by representation theory (the regular rep contains the trivial rep exactly once; the complement has dimension |G|−1). NOT under test.
- **L(1, χ₅) = 2 log(φ)/√5 = 0.43041…** — the Dirichlet L-function at s=1 for the conductor-5 character, via the analytic class number formula for ℚ(√5) with h=1, w=2. A theorem. NOT under test.
- **The Dedekind factorisation** ζ_{ℚ(√5)}(s) = ζ(s)·L(s,χ₅). A theorem. NOT under test.

The formula α⁻¹ = 59/L(1,χ₅) = 59√5/(2 log φ) follows arithmetically from these PLUS one assumption (equipartition, §2). The arithmetic is not in question. The assumption is.

---

## §2. TEST A — Does the spectral gap force equipartition?

### The gap in 157 (its own §7(a), quoted)
157 concedes: *"The step from 'the augmentation ideal is 59-dimensional' to 'the coupling divides equally across 59 directions' requires a physical assumption (equipartition across interaction modes) that is natural but not proved."*

Equipartition is the single unproven step between the forced ingredients and the numerical claim. If the coupling does NOT divide equally across the 59 augmentation modes, α ≠ L(1,χ₅)/59 and the formula is coincidence. So: **is equipartition forced, or assumed?**

### The hypothesis under test
**H_A:** The spectral gap of the relevant operator on the A₅ augmentation ideal is tight enough (the non-trivial spectrum is sufficiently degenerate) that no mode is energetically preferred, forcing uniform distribution of the coupling across all 59 directions — i.e. equipartition is a *consequence* of spectral rigidity, not an independent assumption.

The candidate operator and the candidate gap (to be locked before computation):
- The relevant structure is the **600-cell graph spectrum / the A₅ (or 2I) action on the augmentation ideal**, with λ₁ = 168 named in 157 §6 as "the quadratic floor between vertex and rotation levels" (168 = |PSL(2,7)|, the Klein quartic automorphism order).
- **Tightness criterion (LOCK BEFORE RUN):** equipartition is forced iff the augmentation-ideal spectrum, under the A₅-invariant operator, has its non-trivial eigenvalues equal (degenerate) up to a tolerance that leaves no mode able to carry more than (1/59)(1 + ε) of the coupling, for ε to be fixed below.

### Null hypothesis
**H_A0 (null):** The augmentation-ideal spectrum is NOT uniformly degenerate; the modes carry unequal spectral weight; a generic A₅-invariant operator on the 59-dim ideal distributes weight non-uniformly, so equipartition is an external assumption with no spectral forcing.

### Pre-registered decision rule (LOCK BEFORE RUN)
1. **Compute** the eigenvalue multiplicities of the A₅ regular-representation augmentation ideal under the natural Laplacian / Casimir operator (the operator must be named and fixed before running — proposed: the graph Laplacian of the A₅ Cayley graph on the standard generating set, OR the Casimir of the regular rep; **Cliff/Mr Code lock the choice before run**).
2. **Decompose** the 59-dim augmentation ideal into A₅ irreducibles (expected: 3+3+4+5 = 15 per the irrep dims {1,3,3,4,5}, with the regular rep giving each irrep with multiplicity = its dimension; 59 = 60−1 = (1·1 + 3·3 + 3·3 + 4·4 + 5·5) − 1 = 1+9+9+16+25−1 = 59 ✓). Confirm this decomposition first as a sanity check.
3. **Equipartition forced (H_A passes)** iff: the coupling, distributed by the operator's spectral weights, lands within ε = 1% of uniform (1/59 per direction) across all irrep blocks — i.e. the irrep structure does NOT concentrate weight on any block beyond the 1% tolerance. *(ε = 1% chosen to match 157's own precision ceiling — equipartition need only hold to the precision 157 claims.)*
4. **Equipartition NOT forced (H_A0 stands)** iff any irrep block carries weight deviating from uniform by more than ε = 1%. In that case equipartition remains an assumption and 157 stays OBSERVED.
5. **Single run.** The operator choice and ε are locked before execution. No adjusting ε after seeing the spread.

### What passing buys
If H_A passes: the equipartition step is forced by spectral rigidity, 157's only physical assumption becomes a consequence, and α⁻¹ = 59√5/(2 log φ) is DERIVED at 0.02%. This is the load-bearing test.

### What failing costs
If H_A0 stands: 157 stays OBSERVED. No loss — the formula is still a convergence result; it simply isn't promoted. Honest null.

---

## §3. TEST B — Does higher-tower augmentation close the 0.02% residual WITHOUT a fitted coefficient?

### The gap in 157 (its own §4.2)
157 says the 0.02% discrepancy (137.06 vs 137.036) "is consistent with corrections from higher levels of the tower" — the 2I cover (order 120, augmentation 119) and E₈ (240 roots, augmentation 239) — and defers this to companion Paper 158. **Is that a real prediction or a hope?**

### The hypothesis under test
**H_B:** The next augmentation level(s) of the tower (2I: 120→119, and/or E₈: 240→239) contribute a correction to α⁻¹ of the right SIGN and the right ORDER OF MAGNITUDE to close the 0.02% gap, with the correction's coefficient FORCED by the augmentation dimension alone (no fitted parameter).

### Null hypothesis
**H_B0 (null):** Closing the 0.02% gap from higher tower levels requires a fitted coefficient; the augmentation dimensions 119 / 239 do not by themselves produce a correction of the observed size and sign; any apparent closure is tuning.

### Pre-registered decision rule (LOCK BEFORE RUN)
1. **State the correction form BEFORE computing.** The augmentation principle gives each level a coupling L(1,χ₅)/(|G|−1). The candidate correction to α⁻¹ from the next level is a FORCED function of {59, 119, 239} and L(1,χ₅) with NO free coefficient. The exact functional form must be written down and frozen in this brief before the number is computed. *(Proposed forms to choose ONE from, then lock: (a) additive next-level term L(1,χ₅)/239; (b) ratio correction (1 + 59/239); (c) a forced series Σ over levels. Cliff/Mr Code lock ONE before run.)*
2. **Predict sign first.** Does the chosen forced correction move 137.06 toward 137.036 (down) or away (up)? Lock the predicted sign BEFORE computing magnitude. The 0.02% gap is negative (157 overshoots: 137.06 > 137.036), so a passing correction must be NEGATIVE.
3. **H_B passes** iff: the forced correction (no fitted coefficient) has the correct sign AND brings α⁻¹ within a pre-registered factor of 3 of CODATA's 137.035999 (i.e. residual shrinks from 0.02% to < 0.007%) — order-of-magnitude closure, not precision-matching.
4. **H_B0 stands** iff: the forced correction has wrong sign, OR requires a fitted coefficient to close, OR overshoots/undershoots by more than the factor-of-3 window.
5. **Anti-fishing guard:** only ONE correction form is tested (the one locked in step 1). Testing multiple forms and reporting the best is a Pattern 75 violation and voids the result. If the locked form fails, that is the result.

### What passing buys
A forced, signed, order-of-magnitude closure of 157's residual from the tower's own structure — turning "consistent with higher corrections" (hope) into "predicted by the next augmentation level" (result). Combined with Test A, this would make 157 DERIVED with a structurally-predicted residual.

### What failing costs
The residual stays an open "higher corrections" gesture. 157 stays at 0.02% OBSERVED on the residual question even if Test A passes the equipartition question. Honest — and notably, a wrong-sign or fitted-coefficient result here is exactly the kind of finding that should be *welcomed*, because it stops a numerology trap before it ships.

---

## §4. Phase 2 — POCKETED, NOT ACTIVATED (back-EMF / proton g−2)

Recorded for completeness; **no trigger, do not run under this brief.**

The back-EMF hypothesis (the 2I→S³ self-reference leak returns via the quintic wave / A₅ non-solvability; the interface carries a reactive cost that opposes the leak, Lenz-style; this reactance IS the residual and shows up as the proton/electron g−2 offset ~10⁻⁸) is a candidate Phase 2 test. Its pre-registration, IF activated later, must lock:
- the SIGN of the reactive correction (back-EMF opposes the leak → effective coupling slightly weaker than bare augmentation value) BEFORE comparing to the g−2 offset;
- the scale, order-of-magnitude only;
- explicit acknowledgement that proton/electron g−2 is the highest-precision arena in physics and therefore the highest numerology risk — the test must be a SIGN-and-SCALE structural prediction, never a precision match.

**Not activated. Pocketed per Cliff, 1 June 2026.** Listed so the plan exists; deliberately not triggered.

---

## §5. Execution discipline (applies to Tests A and B)

- **Pre-register, then run.** This brief is committed to the working directory (timestamped) before any computation. Operator choice (Test A), correction form (Test B), ε, sign predictions, and decision thresholds are frozen here.
- **Seed fixed**, single run per locked spec.
- **Report the result of the locked test**, pass or fail. Do not re-spec and re-run to get a pass. A null is a publishable result (cf. the γ HALT-NULL, May 2026 — z=16 analytic collapsed to p=0.15 empirical; pre-registration caught it).
- **Sanity checks first** (Test A step 2's irrep decomposition; Test B's sign-before-magnitude) before the load-bearing computation, so a wrong setup is caught before it contaminates the result.
- **Outputs:** for each test, a HALT-NULL or PASS verdict against the pre-registered threshold, the computation, and the reproduction script. No promotion of 157's status happens until CinC + Cliff review the verdicts.

---

## §6. What the three outcomes mean for the cluster audit

| Test A (equipartition) | Test B (residual) | Consequence for 157 |
|---|---|---|
| PASS | PASS | 157 DERIVED at 0.02% with structurally-predicted residual. 157 becomes the cluster spine; 90/127 base formula is the lossy coordinate image. Strong anti-numerology position. |
| PASS | NULL | 157 DERIVED at 0.02%; residual stays open "higher corrections". Still a major promotion. |
| NULL | PASS | Residual mechanism works but equipartition unforced — 157 stays OBSERVED; investigate why coupling distributes as it does. |
| NULL | NULL | 157 stays OBSERVED. Convergence argument (Pattern 9) still stands; no promotion. The honest floor. |

In every outcome the convergence argument survives (three independent routes onto the same forced object). The tests determine whether 157 is promoted from OBSERVED to DERIVED, nothing weaker.

---

*Mr Code Brief — Paper 157 rigour tests, pre-registration v0.1*
*To be frozen before run. Pattern 75 discipline. 🐕☕⬡*
