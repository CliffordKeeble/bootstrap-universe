# Mr Code Brief — Paper 157 Rigour Test
## Pre-registration v0.2 — the block-weight test (reframed after v0.1 scoping)

**Status: PRE-REGISTRATION. To be committed (timestamped, frozen) BEFORE the locked computation runs.**
**Discipline: Pattern 75 — null and decision thresholds locked before the data is touched. No fishing. Single run per locked spec. Seed fixed.**

**Author:** CinC, for Cliff, 2I Universe Programme
**Date drafted:** 1 June 2026
**Supersedes:** Mr_Code_Brief_Paper_157_rigour_v0_1.md (the two-test version). v0.1's scoping run (commit 0c79bc6) is the input to this reframe and is summarised in §1.
**Target paper:** Paper 157 (The Self-Reference Tower)

---

## §0. Why v0.2 exists — what the v0.1 scoping established

v0.1's lock-free scoping run (pure A₅ rep theory, no locked test executed) produced three results that force a reframe. All three are taken as established input here, not re-litigated:

1. **ARITHMETIC FLAG (Pattern 39, against Paper 157 itself).** 59√5/(2 log φ) = **137.0789**, verified at 40 dp — **not 137.06** as stated in Paper 157 v1.3 (abstract, §4 eq. 3, Table 3, conclusion). True residual against CODATA 137.035999177 is **0.031%**, not the 0.02% the paper claims. Sign unchanged (overshoot). This is a headline error in a live paper and requires a 157 rein-back independent of this test (see §5). The effective channel count that hits CODATA exactly is **58.9815** (a downward correction of 0.0185 from 59).

2. **EQUIPARTITION IS FORBIDDEN, NOT FORCED (kills v0.1 Test A as framed).** By Schur's lemma, any A₅-invariant operator is block-scalar on the augmentation ideal's four pairwise-inequivalent isotypic blocks (irreps 3, 3′, 4, 5; dims summing 9+9+16+25 = 59). Equipartition across all 59 directions would require all four block-scalars equal → a global scalar operator → **no spectral gap**. A real gap (λ₁ = 168) is by definition a *distinction* between modes — the opposite of equipartition. So v0.1's "does the gap force equipartition?" was self-undercutting: a gap forbids equipartition. v0.1 Test A is recorded as **NULL by structure** and not re-run.

3. **THE TOWER RESIDUAL ROUTE IS DEAD (kills v0.1 Test B).** Forced next-level augmentation ratios (1/119, 1/239, …) overshoot the 3.13×10⁻⁴ gap by 10×–1600×. Only squared/product terms reach the right scale, and selecting one is precisely the fitted-coefficient freedom H_B0 forbids. v0.1 Test B is recorded as **NULL by structure** and not re-run.

**What survives.** The non-trivial fact under (2): the coupling across the 59 directions is *not* uniform — Schur says it is block-scalar with four (potentially distinct) weights on the (3, 3′, 4, 5) blocks. The effective count needed for CODATA is 58.98, not 59. So the question is no longer "is the coupling uniform?" (no) but **"are the four block weights FORCED, and do they produce the 58.98 effective count?"** That is the v0.2 test.

---

## §1. The reframed hypothesis — and the trap it must avoid

### The structural picture
The 59-dim augmentation ideal of A₅ decomposes into four isotypic blocks:

| Block | Irrep dim | Multiplicity (= dim, regular rep) | Directions in block |
|---|---|---|---|
| 3 | 3 | 3 | 9 |
| 3′ | 3 | 3 | 9 |
| 4 | 4 | 4 | 16 |
| 5 | 5 | 5 | 25 |
| | | | **Σ = 59** |

By Schur, the coupling operator acts as a scalar wᵢ on each block. The total coupling is L(1,χ₅); the effective inverse coupling (α⁻¹) is L(1,χ₅)⁻¹ weighted by how the blocks carry it. Uniform weights (all wᵢ = 1) give the flat count 59 → 137.0789. The CODATA-matching effective count 58.9815 corresponds to a *specific small non-uniformity* in the {wᵢ}.

### The hypothesis
**H:** The four block weights {w₃, w₃′, w₄, w₅} are FORCED by an independently-specified structure (named and locked in §2 BEFORE computation), and the resulting effective channel count matches 58.98 — i.e. the forced non-equipartition of the coupling across the four A₅ irrep blocks closes the 0.031% residual with NO fitted parameter.

### The trap (stated loudly, because this reframe is seductive)
Replacing one count (59) with four weights {w₃, w₃′, w₄, w₅} is **four knobs instead of one** unless the weights are forced. If Mr Code is free to choose four weights that sum to 58.98, the test is *worse* numerology than the original, not better — it has more freedom, dressed in rep-theory costume. 

**Therefore the test passes ONLY IF the weighting rule is stated and justified BEFORE the effective count is computed, and the rule is forced by structure external to the target number.** If the weighting rule can only be stated by reference to "what lands on 58.98," the test is NULL BY CONSTRUCTION and Mr Code halts and reports exactly that. This is the load-bearing guard of v0.2.

---

## §2. The forcing gate — LOCK BEFORE RUN (this is the whole test)

Before any number is computed, ONE forced weighting rule must be written down here and frozen. Candidate forcing structures (Cliff/CinC select ONE, or reject all → null):

**Candidate W1 — χ₅ character weighting.** L(1,χ₅) already encodes the χ₅ Legendre character, which knows the cooperative/stubborn (split/inert) prime split. If the four A₅ blocks inherit a forced weight from how χ₅ acts on each irrep — e.g. the character value χ₅ restricted to each block's conjugacy-class content — those weights are forced by the L-function that is *already in the formula*, not chosen. **Test: compute the χ₅-induced weight on each block from the character table alone, BEFORE comparing to 58.98.**

**Candidate W2 — dimension weighting.** The blocks carry weight proportional to a fixed function of irrep dimension (e.g. the trivial-rep removal already distinguishes the augmentation ideal; perhaps a second forced removal — the 5-dim block as the "vertex degree" mode — gives a forced deficit). **Test: the deficit 0.0185 ≈ ? a forced dimensional quantity. 0.0185 × 59 ≈ 1.09 — is the correction ≈ 1 mode out of 59, i.e. a SECOND augmentation (removing one more structurally-forced direction)? LOCK the identity of that direction before computing.**

**Candidate W3 — spectral-gap weighting.** The λ₁ = 168 operator assigns each block its eigenvalue; weights are the (forced) eigenvalue ratios. **Test: compute the four eigenvalues of the locked operator on the four blocks, derive weights as their forced ratios, BEFORE comparing to 58.98.** (Note: requires locking the operator — graph Laplacian of A₅ Cayley graph on a NAMED generating set, vs Casimir — and the generating set choice is itself a forcing question, not a free choice.)

**Decision at the gate:**
- If ONE candidate gives a weighting rule statable entirely from structure (character table / dimensions / locked operator spectrum) with NO reference to 58.98 → lock it, proceed to §3.
- If NO candidate can state its weights without peeking at the target → **HALT. Report "NULL BY CONSTRUCTION: no forced weighting available."** This is a legitimate, publishable outcome and the correct one if the weights aren't forced.

**The gate is the test.** §3 is just arithmetic once the gate is passed.

---

## §3. Pre-registered decision rule (LOCK BEFORE RUN)

Given a forced weighting rule from §2:

1. **Compute the four weights from the locked rule alone.** Record them before any comparison.
2. **Compute the effective channel count** N_eff = (the weighted combination, formula fixed in §2) from those weights.
3. **Predict direction first:** the residual is an overshoot (137.0789 > 137.036), so a passing weighting must give N_eff < 59 (downward correction toward 58.98).
4. **H passes** iff: N_eff lands within a pre-registered window of 58.9815 — **window = ±0.003** (i.e. N_eff ∈ [58.978, 58.985], corresponding to closing the residual from 0.031% to better than ~0.005%) AND the weights were stated from structure alone per the §2 gate.
5. **H NULL** iff: N_eff outside the window, OR wrong direction, OR (critically) the weighting rule could not be stated without reference to the target.
6. **Single run.** One weighting rule, one computation, pass or fail. No trying W1 then W2 then W3 and reporting the best — that is multiple comparisons and voids the result (Pattern 75). If the locked rule fails, that is the result. *(If Cliff wishes to test more than one candidate, EACH must be a separately pre-registered single-shot with its own frozen brief and the multiple-comparison correction stated — not a menu explored in one run.)*

---

## §4. Anti-fishing guards (v0.2-specific, beyond standard discipline)

- **The four-knobs trap:** the §2 gate exists precisely to prevent four free weights. A weighting rule that requires four independent choices is NOT forced and fails the gate, even if each choice is "icosahedral." (This is the Paper 90 lesson at the block level: "constructible from icosahedral quantities" ≠ forced.)
- **One rule, pre-stated:** the weighting rule is named in §2 and frozen before §3 runs. The character table / dimensions / operator spectrum are computed; the *rule mapping them to weights* is fixed in advance.
- **Direction before magnitude:** N_eff < 59 predicted before the value is seen.
- **Null is welcome:** "no forced weighting available" is the expected outcome if the residual is genuinely just the honest 0.031% precision ceiling of a convergence result. A null here protects 157 from acquiring a numerology trap; it does not weaken the convergence argument.

---

## §5. Independent of this test — the 157 arithmetic rein-back

**This happens regardless of the test outcome.** Paper 157 v1.3 states α⁻¹ = 137.06 in the abstract, §4 (eq. 3), Table 3, and the conclusion. The correct value is **137.0789** and the correct residual is **0.031%** (not 0.02%). This is a Pattern 39 / Pattern 75 headline error in a live paper.

Required 157 edit (MINOR rein-back, separate from this brief):
- Correct 137.06 → 137.0789 at all four occurrences.
- Correct residual 0.02% → 0.031% wherever stated.
- §4.2's "0.02% gap" framing updated to 0.031%.
- The Table 3 comparison row (59/L(1,χ₅) → 137.06, OBSERVED) corrected to 137.0789.
- Note: this makes 157 ~1.6× further from CODATA than v1.3 claimed — the OBSERVED status is unchanged and honest, but the number must be right.

This is logged as an audit finding against 157 and should be actioned at 157's own rein-back, with or without the block-weight test passing.

---

## §6. Outcomes and what each means

| §2 gate | §3 result | Consequence for 157 |
|---|---|---|
| Forced rule found | N_eff ∈ window | 157's residual is the forced block-weight non-equipartition. Equipartition was the wrong mechanism; block-weighting is the right one and it's FORCED. 157 promotes toward DERIVED at corrected precision. Strong result. |
| Forced rule found | N_eff outside window | The forced rule is wrong; residual unexplained. 157 stays OBSERVED at 0.031%. Honest null. |
| No forced rule | HALT | Residual is the honest precision ceiling of a convergence result, not a forced correction. 157 stays OBSERVED at 0.031%. The convergence argument (Pattern 9) stands; no numerology trap acquired. The expected, safe outcome. |

In every row the convergence argument survives and the arithmetic correction (§5) happens. The test only determines whether the 0.031% residual has a forced structural cause or is the honest cost of self-reference at this geometry (Cliff's "1% is the price" intuition — which, if the gate finds no forced rule, is simply *true* and recorded as such).

---

## §7. Execution discipline

- Pre-register (this brief, frozen, timestamped) before the §2 gate is evaluated against any target number.
- The §2 forcing gate is evaluated FIRST and can halt the whole test with "NULL BY CONSTRUCTION."
- Seed fixed; single run; one weighting rule.
- Report pass/null against the frozen thresholds. Do not re-spec to pass.
- The §5 arithmetic correction is reported separately as a 157 finding regardless of test outcome.
- Outputs: gate verdict, weights (if computed), N_eff (if computed), PASS/NULL against §3, reproduction script. No status promotion until CinC + Cliff review.

---

*Mr Code Brief — Paper 157 rigour test, pre-registration v0.2 (block-weight reframe)*
*Supersedes v0.1. Incorporates v0.1 scoping commit 0c79bc6: arithmetic flag, equipartition-forbidden, tower-residual-dead.*
*To be frozen before run. Pattern 75 discipline. The §2 forcing gate is the test. 🐕☕⬡*
