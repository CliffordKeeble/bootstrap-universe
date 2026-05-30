# Mr Code Sub C brief — Discriminant variation in the Galois-norm probe

**Date**: 30 May 2026
**Brief author**: CinC
**Context**: Mr Adversary reviewed Paper 150 v2.0 (the golden-norm Dedekind detection) and Paper 203 v0.3 (the closure-in-time synthesis) in the same week. Both reviews share a single load-bearing catch: *the discriminant in p² − 5q² has never been varied*. Mr A's specific test: vary the D in p² − Dq² across multiple fields, test against L-functions of multiple characters, pre-register the diagonal that the Dedekind-bridge claim predicts. If the diagonal emerges, Paper 150 earns its 5th star and Paper 203's spine has its empirical anchor; if a full matrix emerges, the field-specificity claim of both papers is decoration and the closure-in-time framework loses its grip on the icosahedron specifically.

A structural prediction has been derived in the corpus since Paper 196 (3 May 2026) but not connected to Paper 150 (12 April 2026): the golden specificity is forced by the Gauss-sum identity 2·cos(π/5) = φ. Sub C tests whether that structural prediction is empirically confirmed.

**Pre-registration document**: commit to git BEFORE any zero data is pulled and BEFORE any computation is run.

**Estimated runtime**: ~half a day of compute. Pure exact arithmetic where possible, mpmath at fixed dps for the L-function zero refinement; same machinery as Paper 150 v2.0 with the discriminant parameterised.

**Reads required, before any computation:**

- This brief (BRIEF.md)
- Paper 150 v2.0 (the existing four-test machinery; what worked and what was pre-registered)
- Paper 196 §5 (the Gauss-sum derivation; the structural prediction for the diagonal)
- Paper 164 §2 (the Dedekind factorisation for ℚ(ρ) at conductor 7; the second concrete field worked out)

**Read NOT required, by hygiene:**

- Mr Adversary's two reviews — Cliff sent them to CinC; Mr Code stays scope-protected per Pattern 97 (debugging needs context; reviewing needs independence). CinC's summary above is sufficient briefing.

---

## Sub C-1 — The canonical matrix (Mr A's specification)

### Question

Does the indefinite Galois norm |p² − Dq²| of ℚ(√d) detect zeros of L(s, χ_D') for D = D' (own field), and fail to detect zeros for D ≠ D' (cross field)?

### Structural prediction (pre-registered, from Paper 196 §5)

Paper 196 §5 derives the field-specificity from the Gauss-sum identity:

> sin(2π/q) / sin(π/q) = 2·cos(π/q)

For q = 5: 2·cos(π/5) = φ (golden).
For q = 7: 2·cos(π/7) = ρ (Steinbach heptagonal, Paper 164 §1.1).
For q = 3: 2·cos(π/3) = 1 (trivial).
For q = 13: 2·cos(π/13) ≈ 1.942 (no special name).

The Gauss sum τ(χ_q) carries this q-specific algebraic content into the functional equation of L(s, χ_q). The indefinite Galois norm |p² − Dq²| of ℚ(√d) carries the d-specific bilinear form structure. For own-field (D = D'), these contents match and detection succeeds; for cross-field (D ≠ D'), they mismatch and detection fails.

**Prediction**: clean diagonal in the (D, D') matrix below; off-diagonal cells |z| < 5 (below discovery threshold).

### Test design

For each fundamental discriminant D and character conductor D', build:

- **Probe**: ζ_φ,D(t) = |Re²(Z_φ(t)) − D · Im²(Z_φ(t))|

where Z_φ(t) is the Fejér-weighted golden-angle Dirichlet sum of Paper 150 v2.0 with phase α = φ (Paper 150 Test 2 showed phase is interchangeable; we use φ for consistency with v2.0 baseline). The probe is the same as Paper 150 v2.0 with the discriminant D parameterised in place of the fixed 5.

- **Target**: the non-trivial zeros of L(s, χ_D') from LMFDB, refined to 30 dps if necessary.

- **Classification window**: W = 1.0 (same as Paper 150 v2.0). Do NOT tune; window is fixed by precedent.

- **Null**: same-density uniform-random null. For each cell, generate N_null = 1000 null samples and compute the null distribution of mean classification distance Δ.

- **Test statistic**: z_overall against null, as in Paper 150 §3.

### Canonical matrix

Norms (rows): p² − Dq² for D ∈ {3, 5, 7, 13}

Characters (columns): zeros of L(s, χ_D') for D' ∈ {3, 5, 7, 13}

For each cell, classify probe minima against L-function zeros at matched N and matched density. Use N_probe = 1000 for first pass (Paper 150 v2.0 N=1000 reference); extend to N = 10⁴ for cells with |z| ≥ 5 to confirm √N scaling.

Note on conductors: real quadratic field ℚ(√d) has discriminant D where D = d if d ≡ 1 (mod 4), else D = 4d. The character mod D is the Kronecker symbol. For our matrix:

- d = 5: D = 5 (Kronecker = Legendre mod 5)
- d = 13: D = 13 (Legendre mod 13)
- d = 3: D = 12 (Kronecker mod 12; we use the formal conductor)
- d = 7: D = 28 (Kronecker mod 28)

To keep the test clean and match Mr A's specification, we parameterise by the squarefree d in the norm (so the probe is |p² − dq²| with d ∈ {3, 5, 7, 13}) and we use the corresponding fundamental L-function (which may be at conductor 4d for d ≡ 2, 3 mod 4). Each cell is documented explicitly with (d, conductor, character) triple before running.

### Pre-registered thresholds

For each cell:

- **PASS** (own-field detection): |z| ≥ 5 (discovery threshold from Paper 150 v2.0).
- **FAIL** (off-diagonal, cross-field): |z| < 5.
- **AMBIGUOUS**: 3 ≤ |z| < 5. Cell flagged for re-running at N = 10⁴.

Bonferroni correction across 16 new cells: α_per = 0.05/16 ≈ 0.003, which gives an adjusted z-threshold of approximately ±2.97. Discovery threshold |z| ≥ 5 is well above this and is unaffected.

### Verdict logic (pre-registered)

- **Diagonal verdict** (PASS on all diagonal cells, FAIL on all off-diagonal cells): predicted outcome from Paper 196 §5. Massive confirmation; field-specificity is structurally derived AND empirically confirmed. Paper 150 v2.1 earns 5th star; Paper 203 v0.4 spine is anchored.

- **Matrix verdict** (PASS on all cells, including off-diagonal): the indefinite-norm-on-critical-line is field-general; Paper 196's mechanism predicts the diagonal but probe detection sees something else. Paper 150 v2.1 retitles ("An Indefinite Norm Hears the Critical Line"); Paper 203 v0.4 must severely re-anatomise. Either way, finding is publishable (it would be a real and interesting structural fact).

- **Mixed verdict** (some diagonal PASS, some FAIL; some off-diagonal PASS): the most informative outcome. Identify which (d, conductor) combinations work and which don't; map onto the Gauss-sum algebraic content; refine the structural prediction. Paper 150 v2.1 retains its claim for the cells that pass; the title is qualified.

- **Diagonal-fail verdict** (some diagonal cells fail): the probe-Dedekind connection is more subtle than Paper 196 §5 predicts. Paper 196's derivation reopened; Paper 150's title bet not won.

### Stop-on-fail

If three or more diagonal cells fail (i.e., the predicted-PASS cells produce |z| < 5), STOP the matrix and report. The structural prediction has failed at the empirical level and we need to understand why before continuing. Do NOT proceed to Sub C-2 in this case.

If the (5, 5) cell — which Paper 150 v2.0 has already established — comes back at |z| < 5 in this re-run, STOP immediately. A failed reproduction of the headline result indicates a methodological problem with the discriminant-parameterised code.

---

## Sub C-2 — Effect size vs N at (5, 5)

### Question

Is the effect size at the (5, 5) cell (which Paper 150 v2.0 reports at z = −22.91 for N = 10⁴) stable across N, or does it decline as N grows?

### Context

Mr Adversary's catch on Paper 150 v2.0: from N = 100 (|z| ≈ 1.8) to N = 1000 (|z| = 8.14) the effect grew above √N prediction (~5.7 expected); from N = 1000 to N = 10⁴ (|z| = 22.91) it shrank below √N prediction (~25.7 expected). This is "rose then fell." z-scores grow as √N for any fixed nonzero effect, so √N scaling is not in itself confirmation of mechanism; what matters is whether the effect size (percentage tightening of mean Δ versus null) is flat across N. If the effect size declines, the result has a horizon problem — it may wash out at N → ∞.

### Test design

For the canonical (d, conductor) = (5, 5) cell — the Paper 150 v2.0 baseline:

- Compute the effect size: ε(N) := (mean Δ_null − mean Δ_observed) / mean Δ_null

  This is the percentage tightening of classification distance. ε ≈ 0.25 at N = 1000 per Paper 150 v2.0 reporting.

- Compute at N ∈ {100, 1000, 10⁴, 10⁵} (the last requires extending Paper 150 v2.0's machinery; if compute budget tight, drop to N = 10⁴ max).

- Report ε(N) per N. Plot ε vs log N.

### Pre-registered thresholds

- **Flat verdict** (ε(N) stable within ±20% across N): no horizon problem; the effect size is preserved as sample grows, consistent with a genuine fixed mechanism. Paper 150's headline z = −22.91 is then a real and stable result. Mr A's concern dissolves.

- **Declining verdict** (ε(N) monotonically decreasing with N, especially if ε(10⁵) < 0.5 · ε(10³)): horizon problem real. The effect washes out as N grows; mechanism is finite-N artefact. Paper 150 v2.1 must report this honestly; the title becomes a bet on the relevant N regime only.

- **Stable-then-declining**: noted but interpretation deferred until the asymptotic behaviour is clearer.

### Stop-on-fail

If Sub C-1 already stopped due to diagonal failure, Sub C-2 does not run.

If Sub C-1 confirmed the diagonal but Sub C-2 shows a declining effect size, both subs are reported: the diagonal is correct (structural) but the effect-size horizon is a separate concern, and Paper 150 v2.1 must address it.

---

## Implementation notes

- **Build on golden-zeta/ infrastructure**. Paper 150 v2.0 has run_N1000.py, run_alt_alpha.py, run_norm_decoupling.py, run_dedekind_partition.py, run_N10000.py. Sub C extends with a parameterised norm coefficient D and parameterised target L-function χ_D'. Add: run_discriminant_matrix.py and run_effect_size_vs_N.py.

- **Zero data**: LMFDB for L(s, χ) zeros at the specified conductors. Use mpmath dps = 30 for refinement. All zeros external; no programme-generated zeros.

- **Reproducibility**: fixed random seed for null generation; fixed dps for refinement; version-pin Sympy, mpmath, and the zero-data version.

- **Reporting**: per cell, log (d, conductor, N_probe, N_null, observed mean Δ, null mean Δ, null std Δ, z-score, PASS/FAIL/AMBIGUOUS verdict).

---

## Hygiene

- **Pre-registration is binding**. Commit BRIEF.md + PRE_REGISTRATION.md before any zero data is pulled. The matrix specification, classification window, thresholds, and verdict logic are fixed at pre-registration. After that, no edits to the rules.

- **Construction choices made BETWEEN pre-reg and execution are allowed and must be documented**: which fundamental discriminants exactly (e.g., if d = 6 has tractable L-function zeros readily available, include it; if not, document why excluded); which LMFDB endpoint is used; how the Kronecker symbol is implemented for composite discriminants. All choices logged in PRE_REGISTRATION.md.

- **No post-hoc rule changes**. If, midway through running, an unexpected pattern emerges (e.g., a different diagonal than predicted), DO NOT redefine the test. Run to completion, report the actual matrix, and let CinC adjudicate.

- **Bonferroni stated up front**: 16 new cells, α_per = 0.05/16. Don't relax it after seeing the data.

- **Anti-circularity check**: all probe minima are computed from the Dirichlet sum; all zero positions from LMFDB; no circular reference where a "match" is built into the data pipeline. Verify by running the (5, 5) cell as a reproduction of Paper 150 v2.0's result; if it reproduces z = −22.91 (or close), the pipeline is sound.

- **Honest reporting**: a clean matrix verdict (refuting the structural prediction) is more valuable than a forced diagonal verdict. Report what comes out.

---

## Deliverable

A markdown report (`findings_paper_203_sub_c.md`) with:

```
## Sub C-1 — The canonical matrix

### Construction choices
- Discriminants used and their conductors
- L-function zero sources and refinement protocol
- Null generation details

### Reproduction of (5, 5) cell from Paper 150 v2.0
- Result at N=1000: z = ? (target ≈ −8.14 from v2.0)
- Pipeline-sound verdict: PASS/FAIL

### Matrix at N=1000
| Norm \ Target | L(χ_3) | L(χ_5) | L(χ_7) | L(χ_13) |
|---|---|---|---|---|
| p² − 3q² | ... | ... | ... | ... |
| p² − 5q² | ... | (5,5) | ... | ... |
| p² − 7q² | ... | ... | ... | ... |
| p² − 13q² | ... | ... | ... | ... |

### Extension to N=10⁴ for diagonal cells
- (3,3), (5,5), (7,7), (13,13) at N=10⁴
- √N scaling check

### Verdict
- Diagonal / Matrix / Mixed / Diagonal-fail
- Bonferroni-corrected interpretation

## Sub C-2 — Effect size vs N at (5,5)

### Construction choices
- N values used; null generation protocol

### Results
| N | mean Δ (observed) | mean Δ (null) | effect size ε(N) |
|---|---|---|---|
| 100 | ... | ... | ... |
| 1000 | ... | ... | ... |
| 10⁴ | ... | ... | ... |
| 10⁵ (if budget allows) | ... | ... | ... |

### Verdict
- Flat / Declining / Stable-then-declining

## Cross-sub flags
Anything interesting that emerged.

## Honest limitations
What wasn't tested; where construction choices may matter.

## Recommended actions for v2.1 and v0.4
Per verdict, what Papers 150 v2.1 and 203 v0.4 should say.
```

---

⌨️ over to you.

The structural prediction is on the record. Paper 196 §5 says the diagonal should appear, with each field's norm hearing its own L-function because each Gauss sum carries field-specific algebraic content (2·cos(π/q) being golden only at q=5). Whether the empirics confirm the derivation is now the question.

Mr A's blade is across two papers and one experiment between us and the answer.

---

## Mr Code's report

*[To be filled in by Mr Code.]*
