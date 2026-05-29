# Pre-registration — Paper 203 v0.2 algebraic verification

**Date**: 29 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 29 May 2026)
**Binding**: this document is committed to git **before any computation has been run**.
**Notes on context**: Paper 203 v0.1 is not in the repo; I am working from
CinC's summary of the §3–§9 structure in the brief, the Paper 191
companion code at `golden-dirac/golden_dirac.py`, and (with Cliff's
explicit override of the brief's hygiene rule) Mr Adversary's review
itself, which Cliff included in the conversation. Mr A's review tells me
the same three deep questions the brief poses, in his own words; I do
not let it influence the construction choices I commit to below.

---

## Construction commitments

### Sub 1 — i_golden = Γ_adj / 5^(1/4)

**(1.a)** Compute symbolically `(Γ_adj/5^(1/4))²` and check equality to
`-I₂` (and lift to `-I₄` via tensor by `I₂`). Test PASS iff Sympy
simplifies to identical zero matrix.

**(1.b) Galois conjugation σ**. Defined as the unique non-trivial
automorphism of `ℚ(√5)/ℚ`: `σ(√5) = −√5`, equivalently `σ(φ) = ψ =
1 − φ`. Implementation: apply Sympy's substitution `√5 ↦ −√5` entrywise.
Test PASS iff `σ(Γ_adj)` equals `−Γ_adj` as matrices over `ℚ(√5)`.

I flag in advance that this test may fail under literal entry-wise σ.
The entries of Γ_adj are `0`, `−1`, `√5`, `0`. Under σ, only the `√5`
flips sign, giving `σ(Γ_adj) = [[0, −1], [−√5, 0]]`. Compare with
`−Γ_adj = [[0, 1], [−√5, 0]]`. The top-right entry differs (−1 vs +1),
so σ(Γ_adj) ≠ −Γ_adj literally. I'll report what the actual
relationship is (it may equal −Γ_seed, suggesting that σ swaps seed↔adj
rather than negating). If literal (1.b) fails, I'll also test a
*combined* operation σ ∘ (transpose) or σ ∘ (similarity) to see whether
any composite operation reverses Γ_adj as expected.

**(1.c)** Compute `R(θ) = exp(θ · Γ_adj)`. Using `(Γ_adj)² = −√5·I`:

```
R(θ) = cos(5^(1/4) · θ) I + sin(5^(1/4) · θ) · (Γ_adj / 5^(1/4))
```

Verify (i) group composition `R(θ₁) R(θ₂) = R(θ₁ + θ₂)`. For
conserved form: search for a 2×2 matrix `M` (over `ℚ(√5)`, possibly
extending to `ℚ(5^(1/4))` if needed) such that
`R(θ)ᵀ M R(θ) = M` for general θ. Equivalently, `Γ_adjᵀ M + M Γ_adj
= 0` (M is in the "skew-symmetric" subspace under Γ_adj). Solve as a
linear system over the coefficients of M.

**(1.d)** Canonical commutator analogue. Search over the building-block
basis {I, K, J, KJ, Γ_seed, Γ_adj, etc.} for two 2×2 matrices X, P with
entries in ℤ[φ] and X = X^† (Hermitian over ℚ(√5)) such that:
`[X, P] = c · Γ_adj` for some non-zero scalar `c ∈ ℚ(√5)`. Brute-force
the basis combinations. The trace-zero constraint means `[X, P]`
cannot be scalar·I, but can be proportional to Γ_adj (which has trace 0).
PASS iff at least one non-trivial solution found.

**(1.e)** Schrödinger analogue. Take H = K (Hermitian, eigenvalues
±1) as a simple test Hermitian. Evolve `|ψ(t+1)⟩ = U · |ψ(t)⟩` where
`U = exp(− i_golden · H · Δt)` for small `Δt`. Compute U symbolically;
verify `U^† M U = M` (M from 1.c). Equivalently, verify that the
golden norm `⟨ψ | M | ψ⟩` is preserved.

**Null candidates for i**: enumerate matrices X with X² = −I₂ over
ℤ[φ]/ℚ(√5) restricted to the building-block basis:
- `J` itself (J² = −I), the "Pauli" imaginary
- `Γ_adj/5^(1/4)` (the candidate)
- Other combinations like `(Γ_seed + Γ_adj)/something` if any
Apply (1.a)–(1.d) to each. Report which candidates satisfy which roles.

### Sub 2 — fourth gamma Γ⁰

**(2.a)** Enumerate Γ⁰ candidates. I commit to the form Γ⁰ = A ⊗ B
where A, B ∈ {I, K, J, Z=KJ, Γ_seed, Γ_adj} (the building-block basis).
This gives 6 × 6 = 36 candidate tensor products. For each:
- Check anticommutation with Γ¹ = K⊗Γ_seed, Γ² = J⊗Γ_seed, Γ³ = I⊗Γ_adj
- Check (Γ⁰)² = +√5·I₄ (for the brief's timelike requirement)
- Report all valid candidates

I expect (provisionally; falsifiable by computation): `Γ⁰ = Z ⊗ Γ_seed`
works (Z = KJ). I'll enumerate exhaustively rather than assume.

**(2.b)** Iteration tower: compute `Γ_seed^n` for `n = 1, 2, 3, 4, 6`
symbolically. Verify the **operator norm** scales as `5^(n/4)`, where
the operator norm is the largest singular value (= `√(largest eigenvalue
of Γ_seed^n · (Γ_seed^n)^†)`). Use Frobenius norm too as a cross-check.

**(2.c)** Compute `exp(t · Γ⁰)` symbolically using `(Γ⁰)² = √5·I`:

```
exp(t · Γ⁰) = cosh(5^(1/4) · t) · I + sinh(5^(1/4) · t) · (Γ⁰ / 5^(1/4))
```

Verify by Taylor series sanity check (first few terms).

**(2.d)** Identification test. For each `n = 1, 2, 3, 4`, find `t_n`
such that `exp(t_n · Γ⁰)` matches `Γ_seed^n` (in some embedded
4×4 form). Specifically:
- The "embedding" of Γ_seed acting on 2-vectors into the 4D Clifford
  algebra: which Γⁱ corresponds to Γ_seed? (Answer: Γ¹ = K⊗Γ_seed
  includes Γ_seed, and the lift of `Γ_seed` to 4×4 is given by
  some piece.)
- Direct comparison: are `Γ_seed^n` (2×2) and `exp(t · Γ⁰)` (4×4) the
  "same orbit" up to embedding?

Likely outcome: they live in different algebras (2×2 vs 4×4) and have
different character (dilation vs orthogonal boost). Report negative
finding clearly if so.

### Sub 3 — signature enumeration

**(3.a)** Enumerate signatures (p, q) with p + q = 4. For each
signature, attempt to construct four mutually anticommuting 4×4 gammas
with the right signature using the building-block basis (A ⊗ B with
A, B ∈ {I, K, J, Z, Γ_seed, Γ_adj}). This gives a finite search space.

Search algorithm:
1. For each candidate Γ ∈ {A ⊗ B : A, B in basis}, compute Γ² and
   classify as +√5·I (timelike), −√5·I (spacelike), or other.
2. For each (p, q), pick p timelike and q spacelike candidates from
   the classified pool.
3. Check pairwise anticommutation.
4. Record all sets of 4 mutually anticommuting gammas matching (p, q).

**(3.b)** For each (p, q) for which a valid construction exists,
compute (Γ⁰ Γ¹ Γ² Γ³)² symbolically. Should be a scalar multiple of I₄.

**(3.c)** Compare the scalar across signatures. If only (1, 3) gives
−5^(3/2)·I₄, Constraint 3 holds (Mr A's circularity flag withdraws).
If multiple signatures give −5^(3/2)·I₄ (or some signatures give the
same value with sign or factor differences), Constraint 3 demotes
from "uniquely selects" to "consistency check".

## Pre-committed thresholds (re-stated, fixed)

- **Sub 1 PASS** if all of (1.a, 1.b, 1.c, 1.d) hold, including a
  non-trivial canonical commutator.
- **Sub 1 PARTIAL** if (1.a), (1.b), (1.c) hold but (1.d) finds only
  trivial X, P.
- **Sub 1 FAIL** if (1.b) literally fails (σ(Γ_adj) ≠ −Γ_adj
  entrywise) — though I note above this is the *expected* literal
  outcome; I will report the *actual* Galois relationship and let
  CinC adjudicate whether it counts as PASS-with-amendment or FAIL.

- **Sub 2 PASS** if (2.d) finds clean identification between
  iteration tower and exp(t·Γ⁰) orbit.
- **Sub 2 PARTIAL** if partial identification.
- **Sub 2 FAIL** if no clean identification.

- **Sub 3 PASS** if only (1, 3) yields (Γ⁵)² = −5^(3/2)·I₄.
- **Sub 3 PARTIAL** if (1, 3) is one of several but the others can
  be excluded by additional constraints.
- **Sub 3 FAIL** if multiple signatures yield −5^(3/2)·I₄
  indistinguishably.

## What I will NOT do

- Adjust the thresholds after seeing results.
- Switch construction choices after seeing data.
- Hide negative findings under reassuring framing.
- Bring Paper 203 v0.1's specific text into my interpretation (I don't
  have it; I only have CinC's structural summary).

## Files

- `subs.py` — implementation of Subs 1, 2, 3 with Sympy exact arithmetic.
- `findings_paper_203_algebra_v0_1.md` — long-form report.
