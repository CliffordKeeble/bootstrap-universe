# Paper 157 Rigour Tests — Pre-Lock Scoping

**Status: PRE-LOCK SCOPING.** This note records the lock-free
representation-theory ground truth and the precise residual numbers that
the two pre-registered locks (Test A operator, Test B correction form)
must be chosen against. **No PASS/NULL verdict is declared here** — that
waits for the frozen lock per Pattern 75. Companion: `scoping_rep_theory.py`,
`scoping_output.txt`.

---

## FLAG 0 — arithmetic error in Paper 157 v1.3's headline number

The paper (and the brief, inheriting it) states α⁻¹ = 59√5/(2 log φ) =
**137.06**, residual **0.02%**. High-precision evaluation (40 dp, mpmath):

```
L(1,χ₅) = 2 ln φ / √5      = 0.4304089409640040...   (matches paper)
α⁻¹ = 59 / L = 59√5/(2lnφ) = 137.0789367615257976...   (paper says 137.06)
CODATA 2022                = 137.035999177
residual (abs)             = 0.04293758...
residual (frac)            = 3.1323e-4 = 0.0313%        (paper says 0.02%)
```

**STATUS: DERIVED (the discrepancy is exact arithmetic).** The formula's
true value is 137.079, not 137.06; its true residual is **0.031%**, not
0.02%. The sign is unchanged (overshoot). The 0.02% figure is consistent
with 137.06 having been used as the base (|137.06−137.036|/137.036 ≈
0.018%), so the slip is in evaluating 59√5/(2 log φ) itself, carried into
the residual. The formula is ~1.6× further from CODATA than v1.3 claims.

**Consequence for the brief:** the scope-guard ceiling (§0, "0.02%"),
Test B's pass threshold (§3, "shrink from 0.02% to <0.007%"), and the
"honest 0.02% precision" framing all reference the wrong base. They
should be restated against 0.031% before Test B is locked. To reach
CODATA the formula needs an effective augmentation count of **58.9815**
(not 59) — a downward correction of 0.0185, fractional 3.13e-4.

---

## SANITY (Test A step 2) — augmentation ideal decomposition: CONFIRMED

A₅ irreps {1, 3, 3′, 4, 5}; Σ dim² = 1+9+9+16+25 = 60 = |A₅|. ✓
Character table verified orthonormal (exact, golden values at 5A/5B). ✓

Regular rep contains each irrep with multiplicity = its dimension.
Removing the one trivial copy gives the augmentation ideal:

| irrep | dim | mult | block dim |
|-------|-----|------|-----------|
| 3     | 3   | 3    | 9  |
| 3′    | 3   | 3    | 9  |
| 4     | 4   | 4    | 16 |
| 5     | 5   | 5    | 25 |

Total = 9+9+16+25 = **59** ✓ — **four pairwise-inequivalent isotypic
blocks.** STATUS: DERIVED.

---

## FLAG A — Test A has a structural tension: gap vs equipartition

By Schur's lemma, any A₅-invariant operator acts as a **scalar on each
isotypic block**, and the four blocks {3, 3′, 4, 5} are inequivalent, so
the four block-scalars are **independent**. Class-sum (central) operators
make this concrete — every non-identity class gives a non-uniform
spectrum across the blocks:

```
class |   3   |  3′   |  4  |  5
  2A  |  -5   |  -5   |  0  |  3     NON-uniform
  3A  |   0   |   0   |  5  | -4     NON-uniform
  5A  | 2+2√5 | 2-2√5 | -3  |  0     NON-uniform
```

**The tension:** equipartition (one common eigenvalue across all 59 modes
→ no mode preferred) requires all four block-scalars equal, i.e. the
operator is a multiple of the identity — which carries **no spectral
gap**. Conversely, a non-trivial gap (the brief names λ₁=168) *is* a
distinction between modes — the opposite of equipartition. **Symmetry
does not force equipartition; it leaves a ≥4-parameter family of
invariant weightings, of which equal-weight is one fine-tuned point.**

So "the spectral gap forces equipartition" (H_A) appears
self-undercutting as posed: the gap and equipartition pull opposite ways.
Any operator with a real gap fails the ε=1% uniformity criterion by
order-1 margins (see the table — deviations are hundreds of %, not 1%).
This is a landscape finding for CinC, **not** a verdict: it suggests Test
A as framed lands on H_A0 (null) for any well-defined gapped operator,
and that the operator-lock may need reframing rather than just choosing
among Cayley-Laplacian / Casimir.

---

## Test B scoping — the forced corrections overshoot the gap by 1–3 orders

Target fractional correction (corrected): **3.13e-4**, negative
(overshoot). The natural forced augmentation quantities:

```
1/119 = 8.4e-3    59/119 = 0.496    1/119² = 7.1e-5    1/(59·119) = 1.4e-4
1/239 = 4.2e-3    59/239 = 0.247    1/239² = 1.8e-5    1/(59·239) = 7.1e-5
L/239 = 1.8e-3
```

The single-level ratios (1/119, 1/239, 59/119, 59/239) are 0.4%–50% —
**10×–1600× larger** than the 0.031% gap. The only quantities near the
right scale are *squared/product* terms (1/119², 1/(59·119)), and
choosing which one is exactly the fitted-coefficient freedom H_B0 warns
about. No "best form" is selected here; this only shows that bridging the
gap from {119, 239} without a suppression coefficient is not obviously
available. The single locked form will be run once, pass or fail.

---

## What needs to happen before the locked run

1. **Restate the residual** as 0.031% (FLAG 0) across the brief, or
   confirm the base number — Test B's thresholds depend on it.
2. **Reconsider / lock Test A's operator** in light of the gap-vs-
   equipartition tension (FLAG A) — is the test asking the right question?
3. **Lock Test B's single correction form** against the corrected 3.13e-4
   target.

Items 1–3 are CinC's call (methodology / hypothesis). Mr Code holds the
single locked run until they are frozen.
