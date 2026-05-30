# Sub A findings — signature enumeration scope extension

**Date**: 30 May 2026
**Brief**: [../BRIEF.md](../BRIEF.md), Sub A section.
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `18bdd47` before any matrix enumeration.
**Verdict**: **OBSTRUCTION-THEORETIC**. (1, 3) Minkowski Clifford algebra cannot embed in M₄(ℚ(√5)) at any coefficient bound. Strongest possible negative result.

---

## Method 1 attempt (brute-force search)

I implemented the bounded-coefficient search in `sub_a.py` and `sub_a_v2.py`:

- Computed W = {M ∈ M₄(ℚ(√5)) : {M, Γ⁰} = 0} as the kernel of the linear endomorphism M ↦ {M, Γ⁰} on M₄.
- W has **dimension 8 over ℚ(√5)** (the W-eigenspace of conjugation by Γ⁰).
- Search space at K = 1 (each of 8 coefficients = a + bφ with a, b ∈ {−1, 0, 1}): 9⁸ ≈ 43 million combinations.

Brute-force Sympy enumeration over 43M combinations is computationally infeasible without algebraic pruning. Rather than attempt this, I **switched to Method 2**, which gives a theoretical answer that subsumes any Method-1 result.

## Method 2 — Hilbert symbol obstruction

### Standard Clifford algebra classification

For a field F of characteristic ≠ 2, the Clifford algebra Cl(p, q) over F has structure determined by the parity of (p − q) mod 8 and the Brauer class of a quaternion algebra.

For **Cl(1, 3) over F**:
> Cl(1, 3) over F ≅ M₂(D_F)  where D_F = (−1, −1)_F is the quaternion algebra

For Cl(1, 3) to be isomorphic to M₄(F), we need D_F to be the *split* quaternion algebra, i.e., D_F = M₂(F). This holds iff the Hilbert symbol (−1, −1)_v = +1 at every place v of F (equivalently, the quaternion algebra is split at every place).

### Hilbert symbol computation for F = ℚ(√5)

ℚ(√5) has **two real places** (corresponding to the two embeddings ℚ(√5) → ℝ: √5 ↦ +√5 and √5 ↦ −√5; both give real numbers since √5 ∈ ℝ).

At each real place, F_v = ℝ. The Hilbert symbol (−1, −1)_ℝ asks whether the equation
> −x² − y² = z²
has non-trivial real solutions. Sum of two non-positive quantities equals a non-negative quantity → only solution is x = y = z = 0 → **(−1, −1)_ℝ = −1**.

Both real places of ℚ(√5) yield Hilbert symbol = −1.

Verification (`sub_a_method2.py` output):
> Hilbert symbol at real place '+': −1
> Hilbert symbol at real place '−': −1
> Non-split at both real places: True

### Direct verification: no i ∈ ℚ(√5) with i² = −1

For c = a + b√5 ∈ ℚ(√5):
> c² = (a² + 5b²) + 2ab·√5.

For c² = −1, equate components:
- Rational part: a² + 5b² = −1. Impossible (sum of squares ≥ 0).
- √5 part: 2ab = 0.

Either way, no rational a, b satisfy c² = −1. Therefore ℚ(√5) does not contain a square root of −1, and the standard splitting representation of (−1, −1) as 2×2 matrices does not embed in M₂(ℚ(√5)).

### Conclusion

D = (−1, −1)_{ℚ(√5)} is **non-split** (a non-trivial element of Br(ℚ(√5))). Therefore:

> Cl(1, 3) over ℚ(√5) ≅ M₂(D), a 16-dimensional simple ℚ(√5)-algebra.
> 
> M₄(ℚ(√5)) is also 16-dimensional but corresponds to the **trivial** Brauer class.
> 
> These two algebras are **distinct elements of Br(ℚ(√5))** and hence not isomorphic.
> 
> Therefore Cl(1, 3) cannot be embedded as a subalgebra of M₄(ℚ(√5)).

## Verdict

**OBSTRUCTION-THEORETIC.**

(1, 3) Minkowski Clifford algebra is *not representable* in M₄(ℚ(√5)) at any coefficient bound K — not just at K = 1, 2, 3 but at *every* K. The obstruction is at the level of Brauer classes, not at the level of explicit coefficients.

This is the strongest possible negative result anticipated by the brief.

## Implications for Paper 203 v0.4 §8

Mr A's catch on Paper 203 v0.3 §8 was that the literal-tensor-product enumeration was scope-limited. Sub A confirms — but *strengthens* — Paper 203 v0.3 §8's obstruction claim. Specifically:

- **Paper 203 v0.3 §8 said**: "in the literal tensor-product basis built from the standard 2×2 building blocks, (1, 3) Minkowski is not constructible, only (2, 2) Klein."
- **Sub A upgrades this to**: "in the full ℤ[φ]-linear span of M₄(ℚ(√5)), and indeed in the full ℚ(√5)-linear span, (1, 3) Minkowski is not constructible. The obstruction is at the Brauer-class level."

For v0.4 §8:
1. State the algebraic obstruction in full generality: "Cl(1, 3) over ℚ(√5) is M₂((−1, −1)_{ℚ(√5)}), a non-trivial Brauer class, distinct from M₄(ℚ(√5))."
2. The (2, 2) algebra IS representable in M₄(ℚ(√5)) (corresponds to the trivial Brauer class, splits as M₂(M₂(F))).
3. The "(2, 2) → (1, 3) signature selection" framing must be reformulated:
   - Within ℚ(√5), we are *stuck* in (2, 2). The Brauer class blocks the transition.
   - To "reach" (1, 3), we would need to extend the field — at minimum to a field where (−1, −1) splits. The simplest such extension is ℚ(√5, i) where i² = −1.
4. **The candidate "baryogenesis as signature selection" framing** (Paper 203 v0.3 §8): if interpreted as "matter content creates a field extension where (−1, −1) splits, allowing the (1, 3) algebra to emerge", this becomes a CONCRETE algebraic claim — adjoining a square root of −1 is the signature-selection step.
5. Alternative reading: the (2, 2) → (1, 3) transition requires a field extension ℚ(√5) → ℚ(√5, i). The square root of −1 is what selects the Lorentzian signature from the Klein-space pre-baryonic algebra. This is **a derivation-style argument**, not a metaphor.

## Honest limitations

- **Method 1 (brute-force search) was not completed.** It would have taken ~12 hours at K = 1 in pure Sympy; with algebraic pruning it could be feasible. I did NOT complete it because Method 2 gives a *stronger* answer (no embedding exists at *any* K) and Method 1 cannot contradict Method 2.
- The Hilbert symbol calculation was done at the real places only. Hilbert reciprocity guarantees that the symbol at a finite place would also need to be ±1 with the product over all places = +1. With 2 real places giving (−1)² = +1, the algebra could ALSO be unramified at all finite places — which is consistent with the standard (−1, −1) Hamilton structure. I didn't compute finite-place symbols; the real-place argument alone suffices to establish non-triviality of the Brauer class.
- The result depends on Cl(p, q) structure theory standard in Atiyah–Bott–Shapiro 1964 and Lam (1973) "The Algebraic Theory of Quadratic Forms". I'm citing the result, not proving it from scratch.

## Pattern flags

- **Pattern 39 (DERIVED vs OBSERVED)**: the Hilbert symbol computation (DERIVED, exact), the impossibility of i² = −1 in ℚ(√5) (DERIVED), and the Clifford algebra classification theorem (DERIVED, classical) all combine to give the OBSTRUCTION-THEORETIC verdict.
- **Pattern 19 (adversary)**: Mr A's v0.2 catch on Sub 3 ("literal tensor basis was not exhaustive") is sustained — but the *strengthened* result is that NO basis works. The obstruction is structural at the Brauer level.
- **Pattern 75 (null)**: not applicable in the usual sense; Method 2 is a positive theoretical result rather than a null comparison.

## Files

- [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
- [sub_a.py](sub_a.py) — Method 1 attempt (28-vector "W basis" with redundancies, search infeasible)
- [sub_a_v2.py](sub_a_v2.py) — Method 1 with proper 8-dim W (still infeasible at K=1)
- [sub_a_method2.py](sub_a_method2.py) — Method 2 Hilbert symbol obstruction (verdict source)
- [findings_paper_203_sub_a.md](findings_paper_203_sub_a.md) — this report

## Compute

Method 2: ~5 seconds. Method 1 attempts: ~10 minutes before judging infeasible.
