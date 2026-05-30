# Pre-registration — Sub A (signature enumeration scope extension)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [../BRIEF.md](../BRIEF.md) (CinC, 30 May 2026, Sub A section)
**Binding**: this document is committed to git **before any matrix
enumeration begins**.

## Construction commitments

**Fixed Γ⁰**: I commit to Γ⁰ = Z ⊗ Γ_seed where Z = K·J = diag(1, −1).
This is the unique literal-tensor-product candidate from Sub 3 of
paper-203-algebra with (Γ⁰)² = +√5·I_4 (commit `b19747d`).

**ℤ[φ]-linear span V**: the 36 tensor products A ⊗ B for A, B in the
basis {I, K, J, Z, Γ_seed, Γ_adj}. Each matrix is a 4×4 over ℚ(√5),
with entries in ℤ[φ] = ℤ + ℤφ. The ℤ[φ]-linear span V is a free
ℤ[φ]-module of rank at most 36 (some tensor products may be linearly
dependent over ℤ[φ]).

**Linear independence check (pre-search)**: I'll verify the rank of V
by stacking the 36 matrices as rows of a 36 × 32 matrix over ℤ[φ]
(each 4×4 = 16 entries × 2 components for ℤ + ℤφ → 32 ℤ-components)
and computing the Smith normal form / row-rank. If rank < 36, the
linear-dependence relations are identified and the search space
shrinks accordingly.

## Method 1 — bounded-coefficient search

**Search variables**: 36 ℤ[φ]-valued coefficients (c_{ij} for the 36
basis elements), each c_{ij} = a + bφ with (a, b) ∈ Z(K) where:
- K = 1: a, b ∈ {−1, 0, 1}  (3 × 3 = 9 values per coefficient, but
  excluding (0, 0) when needed)
- K = 2: a, b ∈ {−2, −1, 0, 1, 2}  (25 values per coefficient)

**Total search space**:
- K = 1: 9^36 ≈ 10^34 combinations. Way too large for brute force.

This brute search is infeasible. I therefore commit to a **smarter
search**: enumerate (Γ¹) one at a time within the anticommutes-with-Γ⁰
subspace W ⊂ V, then for each found Γ¹ enumerate (Γ²) ⊂ W ∩
anticommutes-with-Γ¹, etc.

**Step 1**: Compute W (anticommutes with Γ⁰) by linear algebra over
ℤ[φ] (equivalently over ℚ(√5)). The anticommutator condition
{M, Γ⁰} = 0 is a linear condition on M's coefficients. Solve the
linear system and represent W as a free ℤ[φ]-module with explicit basis.

**Step 2**: Within W, enumerate M = Σ c_i · b_i with c_i ∈ K-bounded
ℤ[φ] coefficients on the W basis. For each such M, check:
  - M² = c · I_4 for some scalar c
  - c = +√5 or c = −√5 (the "timelike" / "spacelike" candidates)

The "M² = c·I_4" check is a polynomial in c_i but since c_i are
bounded, it's enumerable. Filter to the elements M with M² = ±√5·I_4.

**Step 3**: For the filtered candidates, build the anticommutation
graph: M ~ N iff {M, N} = 0. Search for 3-cliques (Γ¹, Γ², Γ³)
with all three having M² = −√5·I_4 (signature (1, 3) requires three
spacelike gammas + the timelike Γ⁰).

**Time budget**: K = 1 search at this scale should be feasible.
The actual cost depends on the dimension of W (likely 16-18 over
ℤ[φ], not 36) and on how many candidates pass M² = ±√5·I in the
K-bounded ball.

## Method 2 — algebraic-classification check

The Clifford algebra Cl(1, 3) over a field F has structure:
- Over ℝ: Cl(1, 3) ≅ M_2(ℍ) (algebra of 2×2 matrices over the
  quaternions), a 16-dimensional real algebra.
- Over ℂ: Cl(1, 3) ⊗_ℝ ℂ ≅ M_4(ℂ).

**Question for Sub A**: is Cl(1, 3) representable in M_4(ℚ(√5))
as a sub-algebra with ℤ[φ]-coefficient generators?

This is a representation theory question. The signature (1, 3) over
F is classified by the discriminant of the form and the Witt
decomposition.

Sketch: Cl(1, 3) over F is isomorphic to a matrix algebra over a
quaternion algebra Q_F. For F = ℝ, Q_ℝ = ℍ (since x² + y² + z² = −1
has no real solutions). For F = ℚ(√5), the Hilbert symbol
(−1, −1)_{ℚ(√5)} determines whether the quaternion algebra is split
or division.

I'll commit to computing the Hilbert symbol at the relevant places
and reporting the conclusion: if Q_{ℚ(√5)} = ℍ (division), then
Cl(1, 3) over ℚ(√5) is M_2(ℍ_{ℚ(√5)}), and M_4(ℚ(√5)) cannot contain
it (since ℍ_{ℚ(√5)} ≠ M_2(ℚ(√5)) in this case). If Q_{ℚ(√5)} splits,
M_4(ℚ(√5)) does contain Cl(1, 3).

## Pre-registered thresholds

- **CONSTRUCTIBLE**: Method 1 finds an explicit (1, 3) clique at any
  bound K. Report the clique. Stop the search.
- **NOT CONSTRUCTIBLE at bound K**: Method 1 exhausts at K without
  finding a clique. Report largest K reached.
- **OBSTRUCTION-THEORETIC**: Method 2 establishes that Cl(1, 3) over
  ℚ(√5) is not representable in M_4(ℚ(√5)) (e.g., quaternion algebra
  is non-split at some prime of ℚ(√5)). Strongest negative result.

## What I will NOT do

- Move Γ⁰ to a different choice if the search fails (Γ⁰ is fixed).
- Expand coefficient bound K beyond 2 if infeasible — report
  "NOT CONSTRUCTIBLE at K=2" honestly.
- Re-interpret the verdict thresholds after seeing results.

## Files

- `sub_a.py` — Sympy implementation: ℤ[φ] coefficient arithmetic,
  W subspace, Method 1 search.
- `sub_a_hilbert.py` — Method 2: quaternion algebra over ℚ(√5),
  Hilbert symbol at relevant primes.
- `findings_paper_203_sub_a.md` — long-form report.
