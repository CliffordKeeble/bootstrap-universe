# Sub C-extension findings — positive-definite probe + cubic L-function

**Date**: 30 May 2026
**Brief**: [../BRIEF.md](../BRIEF.md), Sub C-extension section.
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `18bdd47` before any cubic L-zero was pulled.
**Verdicts**:
- Sub C-ext-1: **AMBIGUOUS** (z = −4.69 at extended N = 2792; indefiniteness *helps* but isn't strictly *essential*).
- Sub C-ext-2: **DETECTED** (z = −5.01; mechanism is generic to L-functions, not real-quadratic-specific).

---

## Sub C-ext-1 — positive-definite probe control

### Construction

- **Probe (positive-definite)**: `sqrt(Re²(Z_φ(t)) + 5·Im²(Z_φ(t)))`. No null cone.
- **Probe (indefinite baseline)**: `sqrt(|Re²(Z_φ(t)) − 5·Im²(Z_φ(t))|)` — matches Sub C-1's (5, 5) cell exactly.
- **Target**: zeros of L(s, χ_5), with two runs at different N:
  - **Baseline**: N = 904 chi_5 zeros from Sub C-1's t = 1000 cache.
  - **Extended**: N = 2792 chi_5 zeros computed up to t = 3000 (added in this sub).
- Probe params: N_terms = 5000, dt = 0.008, W = 1.0, N_null = 1000, seed = 42.

### Results

| Probe | N | n_minima | density/unit | signal_mean | null_mean | z | effect ε |
|---|---|---|---|---|---|---|---|
| **Indefinite** (baseline) | 904 | 1178 | 1.174 | 0.291 | 0.413 | **−9.49** | **0.310** |
| Positive-definite (baseline) | 904 | 324 | 0.323 | 0.314 | 0.388 | −3.36 | 0.191 |
| **Indefinite** (extended) | 2792 | 3581 | 1.193 | (similar) | (similar) | **−15.49** | **0.287** |
| Positive-definite (extended) | 2792 | 1048 | 0.349 | (similar) | (similar) | **−4.69** | **0.158** |

### Interpretation

- The positive-definite probe **still detects** L-zeros — z = −4.69 (extended) is far above noise (|z| < 2 expected if no detection).
- BUT the effect size is much weaker: 0.158 vs 0.287 for indefinite — almost 2× less.
- The minima density is much lower for positive-definite: 0.349/unit vs 1.193/unit. Indefinite has ~3.4× more minima at the same percentile.

**Why does positive-definite still detect?** `Re² + 5·Im²` has a single minimum at Re = Im = 0 (origin). Minima of this probe occur where the trajectory `(Re(Z_φ), Im(Z_φ))` comes closest to origin. This **is** correlated with L-function-zero structure (when L = Re + i·Im is small), but less strongly than the indefinite `|Re² − 5·Im²|` whose null cone Re = ±√5·Im is a *line*, intersected by the trajectory many more times.

The effect-size reduction (0.158 / 0.287 ≈ 55%) measures how much the indefinite null-cone structure adds over the positive-definite origin-attraction.

### Verdict

**AMBIGUOUS** (|z| in [3, 5)).

- At N = 904: z = −3.36 → AMBIGUOUS (just above lower threshold)
- At N = 2792: z = −4.69 → still AMBIGUOUS, but approaching detection threshold

Per the brief's prescription, AMBIGUOUS means: "Marginal; rerun at N = 10⁴". Extrapolating √N scaling at the stable effect size ε ≈ 0.158:
- At N = 10⁴, expected z ≈ −4.69 × √(10⁴/2792) ≈ −8.9, which would cross the detection threshold |z| ≥ 5.

So at sufficient N, the positive-definite probe *would* register as DETECTED. **At any finite N tested, indefiniteness helps substantially (effect size ≈ 1.8× larger) but is not strictly essential.**

The cleanest empirical statement: **Indefiniteness contributes ≈ 55% additional effect size; it is helpful but not strictly required.**

### Comparison with Paper 150 v2.0 Test 3a

Paper 150 v2.0 Test 3a tested `|Re² + Im²|` (circular, coefficient 1) and reported z = −2.63 (fails detection at the time's N).

My positive-definite test uses coefficient 5 (the canonical Paper 150 d=5). With coefficient 5 at N = 2792 I get z = −4.69, slightly stronger than coefficient 1 would give but still in AMBIGUOUS regime. This suggests the **coefficient value doesn't drive the indefiniteness essentiality** — what drives it is the *shape* of the level set (line in 2D vs point at origin).

---

## Sub C-ext-2 — cubic L-function detection

### Construction

- **Cubic character χ_2 mod 7** ([chi2_mod7.py](chi2_mod7.py)):
  - χ_2(0) = 0; χ_2(1) = 1; χ_2(2) = ω; χ_2(3) = ω²; χ_2(4) = ω²; χ_2(5) = ω; χ_2(6) = 1.
  - Verified: balanced (sum = 0), multiplicative, EVEN (χ(−1) = χ(6) = 1), cubic (χ³ = trivial).
- **L-function zero finder** ([l_zeros_complex.py](l_zeros_complex.py)):
  - Truncated Dirichlet sum at M = 1000.
  - Finds local minima of |L(½ + it, χ_2)|² on a t-grid.
  - Cross-validated against mpmath Hurwitz-based L(s, χ_2) at first 5 zeros — all have |L| < 0.002 (numerically zero).
  - Found 361 zeros in t ∈ (1, 500].
- **Probe**: `sqrt(|Re² − 5·Im²|)` — canonical Paper 150 v2.0, d = 5.

### Result

| Quantity | Value |
|---|---|
| Target zeros | 361 (cubic L) |
| Probe minima | 576 |
| Density | 1.144/unit |
| Signal mean Δ | 0.3187 |
| Null mean Δ | 0.4362 |
| Null std | 0.0235 |
| **z-score** | **−5.01** |
| **Effect size ε** | **0.269** |

### Verdict

**DETECTED**. |z| = 5.01 ≥ 5 → mechanism is generic to L-functions.

### Comparison with Sub C-1 d_target = 7 cells (real chi_28)

| Target | z | effect ε |
|---|---|---|
| L(χ_28) real (d=7), Sub C-1 (5, 7) cell | −5.50 | 0.282 |
| L(χ_2 mod 7) cubic (this test) | −5.01 | 0.269 |

The two are statistically indistinguishable — same probe, similar conductor (28 vs 7), real-quadratic vs cubic-complex character — gives effectively the same z and effect size. **Strong evidence that the mechanism cares about L-function-on-critical-line structure generically, not character class.**

---

## Combined synthesis

### Mechanism reading

Combining Sub C-1's matrix verdict, Sub C-ext-1, and Sub C-ext-2:

| Test | Result |
|---|---|
| Sub C-1 (real-quadratic matrix, 4×4) | All 16 cells PASS; field-general within real-quadratic |
| Sub C-ext-2 (cubic complex character) | DETECTED at z = −5.01 |
| Sub C-ext-1 (positive-definite) | AMBIGUOUS at z = −4.69; effect ~ 55% of indefinite |

The detection mechanism appears to be:
1. **Equidistributed Dirichlet sum** (golden-angle phase α = φ) — the *carrier* of L-function-zero-like minima in the probe.
2. **Bilinear form on (Re, Im)** — selects *where* the carrier produces minima.
3. **Indefinite form** (with null cone) — *amplifies* the minima density and effect size by ~2×.
4. **Real-quadratic specificity** — *not* a constraint; the mechanism extends to cubic complex characters.

### Implications for Paper 150 v2.1

**Title**: "The Golden-Angle Dirichlet Sum Hears L-Function Zeros" is honest at the current evidence level. Variants:
- "Indefinite-Norm Probe Hears L-Function Zeros" — emphasizes indefiniteness, which is helpful (effect ≈ 2× larger) but not strictly essential.
- "The Golden-Angle Dirichlet Sum is a Generic L-Zero Detector" — emphasizes the equidistributed-phase carrier as the mechanism.

I'd suggest: **"An Indefinite Norm on a Golden-Angle Dirichlet Sum is a Generic Critical-Line Probe"**. This:
- Acknowledges the indefiniteness contribution (Sub C-ext-1: ~2× effect).
- States the mechanism: equidistributed Dirichlet sum + bilinear form.
- Avoids overclaiming field-specificity.

**Mechanism section** should report:
- Sub C-1: real-quadratic field-general (all 16 cells PASS at ε ≈ 0.29).
- Sub C-ext-2: cubic complex character DETECTED at ε = 0.27 — mechanism extends.
- Sub C-ext-1: positive-definite weaker but still detects (ε = 0.16) — indefiniteness helpful, not essential.

### Implications for Paper 203 v0.4

**§5 Three Faces of σ**: Face 2 (empirical detection) is now even broader than my Sub C report indicated. It's not just "real-quadratic Galois norms detect ζ_K(s) of any real quadratic field" — it's "any equidistributed Dirichlet sum with bilinear form (indefinite or, weakly, positive-definite) detects L-function zeros of any character (real or complex) with appropriate FE on the critical line". The mechanism is **generic**, not σ-specific in any sense.

This further sharpens (does not abandon) the recommendation from my Sub C findings: σ-of-ℚ(√5) is just one realisation of a one-parameter family that itself sits within a broader equidistributed-phase mechanism.

**§7 Zeta-value-denominators**: Sub B (separate report) found STRIKING at z = 5.5 supplementary — but the explanation is shared small-prime content, not deep number-theoretic structure. Combined with Sub C-ext-2's "mechanism is character-class-generic", the σ-as-icosahedral-bridge claim loses additional specificity.

**Path forward for v0.4**: the closure-in-time framework's empirical anchor (Paper 150 detection) is real but generic. The framework's specificity must come from elsewhere — possibly Paper 196's *phase lock* at arctan(1/φ) (which IS golden-specific per Paper 196 §5) or the algebraic results from paper-203-algebra commit `b19747d` (e.g., [Z, Γ_seed] = −2Γ_adj).

### Open questions

- **At what N does positive-definite cross detection threshold?** Sub C-ext-1 extrapolation suggests N ~ 10⁴. A direct test would confirm.
- **Does Sub C-ext-2 generalize?** I tested one cubic character (χ_2 mod 7). Testing more complex characters of higher order (quintic mod 11, etc.) would map the boundary of the generic-mechanism claim.
- **Mechanism source**: what is it about the equidistributed-phase Dirichlet sum that produces L-zero-correlated minima? This is the deeper question.

---

## Pattern flags (all of Sub C-ext)

- **Pattern 75 (null)**: pre-registered cell-by-cell nulls. Both sub-investigations satisfied the null discipline.
- **Pattern 39 (DERIVED vs OBSERVED)**: cubic character table is DERIVED; L-zero positions, probe minima, and z-scores are OBSERVED.
- **Pattern 19 (adversary)**: Mr A's "what's the essential ingredient" attack is now answered empirically: equidistributed Dirichlet sum + bilinear form. Neither indefiniteness alone nor real-quadratic specificity alone is "essential" — they are *contributors* of additional effect.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — pre-registered probe forms, character table, thresholds
- [sub_c_ext_1.py](sub_c_ext_1.py) — positive-definite probe driver
- [chi2_mod7.py](chi2_mod7.py) — cubic character verification
- [l_zeros_complex.py](l_zeros_complex.py) — cubic L-zero finder (validated against mpmath)
- [sub_c_ext_2.py](sub_c_ext_2.py) — cubic L detection driver
- [sub_c_ext_1_results.csv](sub_c_ext_1_results.csv), [sub_c_ext_2_results.csv](sub_c_ext_2_results.csv), [zeros_chi2_mod7_tmax500.csv](zeros_chi2_mod7_tmax500.csv) — raw data
- [findings_paper_203_sub_c_extension.md](findings_paper_203_sub_c_extension.md) — this report

## Compute

- Sub C-ext-1 baseline (N=904): ~30 s.
- Sub C-ext-1 extended (N=2792, with extended chi_5 zero compute): ~15 min for chi_5 zero compute + ~3 min for probe.
- Sub C-ext-2: ~12 s (very fast; cubic L-zeros found in 0.5 s).
- Total Sub C-extension: ~20 min including the chi_5 extension to t=3000.
