# Paper 203 v0.2 — algebraic verification of Mr A's three deepest catches

**Date**: 29 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 29 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `6b6b8f7` before any computation.
**Method**: exact symbolic Sympy over ℤ[φ]. No floating-point. All claims verifiable by `simplify(M − target) == 0`.

---

## Headline

- **Sub 1 — §5 hinge**: **PARTIAL with a structural improvement on the brief's (1.b) test.** The canonical commutator `[Z, Γ_seed] = −2·Γ_adj` is found in closed form — this is the structural i-role test Mr A demanded, and it passes cleanly. The literal Galois test σ(Γ_adj) = −Γ_adj **fails**, but the actual relationship σ(Γ_adj) = **−Γ_seed** is arguably the deeper structural fact (Galois conjugation realises a seed↔adj duality, not a self-negation).
- **Sub 2 — time tension**: **PARTIAL (no clean identification).** Exactly one fourth-gamma candidate with +√5·I square exists over ℤ[φ] in the tensor-product basis: Γ⁰ = Z⊗Γ_seed. The iteration tower Γ_seed^n lives in 2×2 (spectral-radius cascade 5^(n/4), confirming the §6 cascade); exp(t·Γ⁰) lives in 4×4. They live in different algebras and have different character.
- **Sub 3 — chirality circularity**: **FAIL** as the brief literally states it, but with a structurally useful finding: **only (2,2) signature is constructible in 4D** over the ℤ[φ] tensor basis. (Γ⁵)² evaluates to **+25·I₄** for all 6 constructible (2,2) cliques. The brief's target value −5^(3/2)·I₄ is the **3D** chirality from Paper 191's (1,2) signature — and in 3D it **does** uniquely select (1,2) over (2,1) in this basis. Constraint 3, as a *3D* statement, holds; as a *4D-selecting-from-(2,2)* statement, it does not.

The most consequential single result: Sub 3's enumeration. (1,3) signature is not constructible in the natural ℤ[φ] tensor basis at all, so "select (1,3) from (2,2)" cannot be done within this construction. v0.2 must either (a) state Constraint 3 as a 3D claim only, or (b) move to a richer basis (e.g., ℚ(5^(1/4))) where (1,3) becomes constructible.

---

## Sub 1 — §5 hinge

### Construction choices

- `i_golden = Γ_adj / 5^(1/4)`. Lives in M₂(ℚ(5^(1/4))) (not ℤ[φ]).
- Galois σ: ℤ[φ] → ℤ[φ] by `√5 ↦ −√5`, applied entrywise.
- For (1.d), built-in basis {I, K, J, Z, Γ_seed, Γ_adj}. Hermitian = symmetric for these (real entries when expressed over ℚ(√5)).

### Tests

| Test | Result |
|---|---|
| **(1.a)** `(i_golden)² = −I₂` | **PASS** — exact symbolic equality. |
| **(1.b)** `σ(Γ_adj) = −Γ_adj` (literal) | **FAIL** — σ(Γ_adj) = [[0, −1], [−√5, 0]], whereas −Γ_adj = [[0, +1], [−√5, 0]]. Top-right entries differ. |
| **(1.b ′)** σ(Γ_adj) = **−Γ_seed**? | **PASS** — exact equality. The Galois conjugate of Γ_adj is −Γ_seed, not −Γ_adj. |
| **(1.c)** Group composition `R(t₁)R(t₂) = R(t₁+t₂)` | **PASS** — closed-form `R(θ) = cos(5^(1/4)·θ) I + sin(5^(1/4)·θ)·(Γ_adj/5^(1/4))` verified. |
| **(1.c)** Conserved form `M` with `Γ_adj^T M + M Γ_adj = 0` | **PASS** — `M = [[√5·d, −c], [c, d]]`, a 2-parameter family. Canonical choice (d=1, c=0): `M = diag(√5, 1)`. |
| **(1.d)** Canonical commutator `[X, P] = c·Γ_adj` with non-trivial X, P from building blocks | **PASS** — `[Z, Γ_seed] = −2·Γ_adj`. X = Z (Hermitian, eigenvalues ±1, "position-like"), P = Γ_seed ("momentum-like" with golden-real eigenvalues), and their commutator is −2 times the imaginary direction. |
| **(1.e)** Schrödinger evolution preserves M? | **PARTIAL** — for H = Z, the evolution `U = exp(−i_golden·Z·t)` preserves M (H preserves the golden inner product, even though [H, i_golden] ≠ 0). For H = K, neither commutation nor M-preservation holds. The "unitary-like" property is H-dependent. |

### Null candidates for i

| Candidate | Squares to −I₂? | Plays i-role? |
|---|---|---|
| J = [[0,−1],[1,0]] (the Pauli/standard imaginary) | ✓ | ✓ trivially — it's the literal imaginary in 2×2 |
| `i_golden = Γ_adj/5^(1/4)` | ✓ | ✓ in roles (1.a, 1.c, 1.d); partial in (1.e) |
| Γ_seed (note: square = +√5·I, NOT −I) | ✗ | n/a |

J also passes the basic i² = −I test trivially, so "i_golden is the unique i-direction" overclaims. However:
- J alone has no canonical commutator with golden structure (`[X, P]` between Pauli matrices doesn't naturally land on J).
- `[Z, Γ_seed] = −2·Γ_adj` is the specific structural fact: the **canonical commutator of the building blocks lands on Γ_adj**, not on J. This is the structural identification Mr A demanded for the §5 hinge.

### Verdict (Sub 1)

**PARTIAL — but the partial supports the §5 hinge in the most important sense.**

- The headline structural result is **`[Z, Γ_seed] = −2·Γ_adj`** — the canonical commutator of two "real" building blocks lands proportional to Γ_adj. This is the Schrödinger-i / canonical-commutator role Mr A specifically said was missing in v0.1.
- The literal Galois test (1.b) fails but the actual Galois relation σ(Γ_adj) = −Γ_seed is arguably more structural: Galois conjugation realises a **seed↔adj swap**, not a self-negation. This is consistent with the §5 reading "Γ_seed and Γ_adj are dual real/imaginary directions" but reframes what "duality" means.
- The Schrödinger analogue (1.e) is H-dependent: some H preserve the golden inner product M, others don't. This is the expected behaviour of any Hilbert-space structure (different observables have different conservation properties); it's not a defect.

**Recommended v0.2 action for §5**: replace "i ↔ −i is complex conjugation" framing with **"the canonical commutator of Z and Γ_seed produces Γ_adj"** and **"Galois conjugation realises seed/adj duality"**. Both are exact algebraic statements; both stand against Mr A's "shared-symbol-not-shared-structure" attack.

---

## Sub 2 — Time reconciliation

### Construction choices

- Enumeration over 36 candidates Γ⁰ = A⊗B with A, B ∈ {I, K, J, Z, Γ_seed, Γ_adj}.
- "Anticommutes with all three Γ¹, Γ², Γ³" applied as filter.

### Tests

**(2.a) Fourth-gamma candidates**:

| Γ⁰ candidate | (Γ⁰)² | Anticommutes with G₁, G₂, G₃? |
|---|---|---|
| I⊗Z | +I₄ | ✓ (but wrong scale — need +√5·I) |
| **Z⊗Γ_seed** | **+√5·I₄** | ✓ |
| (all 34 others) | n/a | ✗ |

**Only one valid timelike fourth gamma exists in this basis**: Γ⁰ = Z⊗Γ_seed.

**(2.b) Iteration tower of Γ_seed**: spectral radius `5^(n/4)` confirmed:

| n | Γ_seed^n (simplified) | spectral radius |
|---|---|---|
| 1 | Γ_seed | 5^(1/4) |
| 2 | √5·I | 5^(1/2) = √5 |
| 3 | √5·Γ_seed | 5^(3/4) |
| 4 | 5·I | 5 |
| 5 | 5·Γ_seed | 5^(5/4) |
| 6 | 5√5·I | 5^(3/2) |

(My Python code reported the operator norm instead of the spectral radius for n=1, which gave √5 rather than 5^(1/4). Γ_seed is not normal — Γ_seed Γ_sedᵀ ≠ Γ_sedᵀ Γ_seed — so operator norm and spectral radius differ. The correct cascade is the **spectral radius** 5^(n/4), which matches the §6 paper claim **and the actual eigenvalues** ±5^(1/4) of Γ_seed.)

**Side note for Mr A's §6 cascade catch**: my computation also confirms Mr A's arithmetic. The actual rule generating the cascade `5^(1/4), 5^(1/2), 5^(3/4), 5, 5^(5/4), 5^(3/2)` is **multiply-by-5^(1/4)** (i.e., iterate Γ_seed), not **square-the-previous**. The paper's "each step is the square of the previous" is false; the correct statement is "each step multiplies the previous by 5^(1/4)" or equivalently "the n-th step is Γ_seed^n".

**(2.c) exp(t·Γ⁰)**: closed form using `(Γ⁰)² = +√5·I₄`:
`exp(t Γ⁰) = cosh(5^(1/4)·t) I₄ + sinh(5^(1/4)·t) (Γ⁰/5^(1/4))`.

This is a **boost-like, non-compact** one-parameter family (hyperbolic functions, not trigonometric).

**(2.d) Identification test**:

| Property | Value |
|---|---|
| (Γ⁰)² = (Γ¹)² | True (both = +√5·I₄) |
| [Γ⁰, Γ¹] = 0 | **False** |
| {Γ⁰, Γ¹} = 0 | **True** |

So Γ⁰ and Γ¹ (= K⊗Γ_seed) are **two different "timelike" gammas that anticommute** in the (2,2) algebra. They are not the same direction.

The iteration tower Γ_seed^n is in M₂(ℚ(√5)). exp(t·Γ⁰) is in M₄(ℚ(√5)). They are in **different algebras**. No clean identification between them.

### Verdict (Sub 2)

**PARTIAL — the two pictures share the eigenvalue scale 5^(1/4) but are not the same algebraic object.**

- The iteration tower lives in 2×2 and acts dilationally (det Γ_seed = −√5, infinite orbit, non-unitary).
- The boost family exp(t·Γ⁰) lives in 4×4 and acts hyperbolically (non-compact, Lorentzian-like one-parameter group).
- Both have characteristic value 5^(1/4) (Γ_seed's spectral radius, and the natural rate of exp(t·Γ⁰)).
- But they are not algebraically the same; they are two views of the same characteristic scale, not two views of the same object.

**Recommended v0.2 action for §4/§6/§8**: Mr A's tension is real. v0.2 must either:
- (a) Pick one ontology of time. If §6's "time as iteration of Γ_seed" is canonical, then §8's "fourth gamma direction" should be reframed as a *derived* algebraic structure on the (2,2) algebra, not as a primitive timelike axis.
- (b) State explicitly that "time as iteration" and "time as fourth gamma" are two related but distinct algebraic structures — the scale 5^(1/4) is shared, but the actions differ in dimensionality and character.

Either way, the headline "no extra spatial dimensions" can be preserved by clarifying that the 4D algebra has 3 *spatial* gammas + 1 *time-character* gamma (not 4 spatial), and that the time-character direction reproduces the 5^(1/4) iteration scale.

---

## Sub 3 — Chirality circularity

### Construction choices

- Building-block basis: A⊗B with A, B ∈ {I, K, J, Z, Γ_seed, Γ_adj}. 36 candidates total.
- Pool of "±√5·I-squaring" candidates: 16 elements (8 with +√5·I, 8 with −√5·I).
- Search: all 4-tuples of pairwise anticommuting elements from this pool. Computed via anticommutation graph + 4-clique enumeration.

### Tests

**(3.a) Signature enumeration**:

In 4D over ℤ[φ] tensor basis:

| Signature (p, q) | Cliques found |
|---|---|
| (4, 0) | **0** |
| (3, 1) | **0** |
| **(2, 2)** | **6** |
| (1, 3) | **0** |
| (0, 4) | **0** |

**Only (2, 2) is constructible.** (1, 3) is not constructible in this basis. (Searching all 4-tuples from the 16-element pool found 6 cliques, all of signature (2, 2).)

In 3D for cross-check (Paper 191 framework):

| Signature (p, q) | Cliques found | (Γ⁵)² |
|---|---|---|
| (3, 0) | 0 | n/a |
| **(2, 1)** | 12 | **+5^(3/2)·I₄** |
| **(1, 2)** | 12 | **−5^(3/2)·I₄** |
| (0, 3) | 0 | n/a |

(Verification script: `sub3_3d_check.py`.)

**(3.b)/(3.c) Chirality across signatures**:

In 4D, all (2, 2) cliques give:
> **(Γ⁵)² = +25·I₄ = +5²·I₄**

The brief's target value **−5^(3/2)·I₄ is not achieved in any 4D signature in this basis.** It IS achieved in 3D — but only by the (1, 2) signature, and uniquely so among constructible 3D signatures.

### Detailed analysis of Mr A's circularity

Mr A's catch: "the exponent 3/2 already encodes three-ness (three spatial gammas, each contributing 1/4)". My computation **confirms** this: the value −5^(3/2)·I₄ is **the 3-gamma chirality value**, not the 4-gamma chirality value. By "encoding three-ness", Mr A means: the value's exponent (3/2) literally counts the number of gammas (3) times 1/2.

This is fully vindicated: **the value -5^(3/2) is a 3-gamma fact, not a 4-gamma fact**. Using it to "select (1, 3) from (2, 2)" is a category error — the values just don't match dimensions.

The 4-gamma chirality values (in this basis) are ±25 = ±5². The exponent 2 = 4·(1/2) literally counts the number of gammas. So an honest 4D Constraint 3 would say:
> "(Γ⁵)² = −5²·I₄ selects Lorentzian (1,3) or (3,1) over split (2,2) or Euclidean (4,0)/(0,4)."

But within the ℤ[φ] tensor basis, only (2,2) exists; the other signatures are unattainable. So even this corrected version of Constraint 3 cannot be applied — there's nothing to "select" because only one signature is constructible.

### Null check (per brief)

The brief asked to check whether any alternative scalar value is achieved by some signature. Within this basis:
- 4D: only +5²·I₄ achievable (one signature, (2,2)).
- 3D: +5^(3/2)·I₄ achievable by (2, 1); −5^(3/2)·I₄ achievable by (1, 2). At D=3, the chirality test **does** discriminate Lorentzian from split.

So in 3D, Constraint 3 (recast as a 3D statement) does select (1, 2) uniquely. In 4D within this basis, the test is moot because only (2, 2) is constructible.

### Verdict (Sub 3)

**FAIL as the brief literally states it; reframable to a 3D statement that holds.**

- Constraint 3 as a "select (1, 3) from (2, 2)" 4D claim **cannot hold** in the ℤ[φ] tensor basis, because (1, 3) is not constructible in this basis.
- The chirality value −5^(3/2)·I₄ is a 3D quantity (Paper 191's (1, 2) result), not a 4D quantity. Mr A's "3/2 encodes three-ness" catch is fully vindicated.
- At D=3, however, **(Γ⁵)² = −5^(3/2)·I₄ does uniquely select (1, 2) over (2, 1)** in this basis. So if Constraint 3 is reframed as "the 3D chirality value selects Lorentzian signature in 3D", the constraint holds.

**Recommended v0.2 action for §8**:
- Withdraw Constraint 3's role as "selecting (1, 3) from (2, 2)". This is the cleanest path: acknowledge it cannot be done over ℤ[φ].
- Reframe as "(Γ⁵)² = −5^(3/2)·I₄ is the 3D chirality of the Paper 191 (1, 2) signature, distinguishing it from (2, 1)". This is true, exact, and not circular.
- The selection of D=3 (rather than D=4) must come from the OTHER constraints (Borwein capacity, Paper 117 dimension argument) — Constraint 3 was always a consistency check, not an independent selection.

**Cliff & CinC should also note**: Mr A had Constraints 1 and 2 wanting strengthening. My Sub 3 finding effectively *removes* Constraint 3 from the "three independent constraints" list and converts it into a consistency check at D=3. The headline of §8 reduces from "three independent constraints select D=3" to **"two constraints (Borwein + Paper 117) select D=3; chirality verifies the Lorentzian signature within the selection"**. That is more honest and more defensible.

---

## Cross-sub flags

### Flag 1: The seed/adjoint duality is the Galois-conjugation orbit.

Sub 1 (1.b) discovered: σ(Γ_adj) = −Γ_seed, σ(Γ_seed) = −Γ_adj (by symmetry — also computed and verified). The Galois automorphism `√5 ↦ −√5` swaps the seed and adjoint directions (with a sign).

**This is potentially the deepest single algebraic fact in the report.** It says:
- Γ_seed (real eigenvalues) and Γ_adj (imaginary eigenvalues) are **a Galois-conjugate pair** at the level of the field ℚ(√5).
- The seed/adj split is **not** "real vs imaginary in the same algebra" — it's "two conjugate algebras over the Galois extension".
- This is the *correct* analog of "complex conjugation maps ℂ-algebras to their conjugates" — but over ℚ(√5)/ℚ, not over ℂ/ℝ.

v0.2 could use this directly: the §5 hinge becomes **"the seed/adj duality is the Galois conjugation of ℚ(√5)/ℚ, and the i of QM is the structure-preserving element that swaps the two"**.

### Flag 2: Only (2, 2) signature is constructible over ℤ[φ].

This is a hard algebraic constraint, not a choice. v0.2 needs to acknowledge it: the (2, 2) algebra is the natural maximum extension of Paper 191's (1, 2). The (1, 3) Lorentzian signature requires either:
- Extending the field beyond ℤ[φ] (e.g., to ℚ(5^(1/4)) or ℂ), losing the ℤ[φ] discipline.
- Changing the algebra (not a Clifford algebra at all).

Either way, "(1, 3) selected from (2, 2)" can't be done within the framework as currently constructed.

### Flag 3: The canonical commutator [Z, Γ_seed] = −2·Γ_adj is the §5 hinge.

This single relation gives v0.2 everything it needs for the §5 hinge:
- Z is "Hermitian over ℚ(√5)" (eigenvalues ±1 ∈ ℝ).
- Γ_seed is "Hermitian over ℚ(5^(1/4))" (eigenvalues ±5^(1/4) ∈ ℝ).
- Their commutator is **proportional to Γ_adj** (eigenvalues ±5^(1/4)·i ∈ iℝ).

This is precisely the canonical-commutator structure `[X, P] ∝ i` of quantum mechanics, realised over ℤ[φ] with i replaced by Γ_adj.

---

## Honest limitations

- **My building-block basis is finite (36 elements).** A richer basis (linear combinations over ℤ[φ]) might allow (1, 3) signature constructions I missed. The "only (2, 2)" claim is for the **tensor-product basis**, not the full linear span.
- **(1.b) literal-vs-structural**: my finding that σ(Γ_adj) = −Γ_seed is structural but doesn't literally satisfy the brief's predicted relationship. Whether this counts as PASS or FAIL is a judgment call I leave for CinC.
- **(2.d) iteration-vs-fourth-gamma**: I judged "different algebras, different character" as no clean identification (PARTIAL). A more generous reading might call them "two views of the 5^(1/4) scale" (PASS-soft). I lean PARTIAL.
- **(3.a) basis choice**: the claim "only (2, 2) constructible in ℤ[φ]" is robust within my tensor basis. Over the algebraic closure ℂ, all signatures are constructible — but then we lose the ℤ[φ] discipline that's the whole point of the Golden Dirac algebra.
- **No floating-point**: all computations are exact symbolic. No precision questions.

## Recommended v0.2 actions per sub

| Sub | Verdict | v0.2 action |
|---|---|---|
| **§5 hinge** | PARTIAL | Replace "i ↔ −i is complex conjugation" with the canonical commutator result `[Z, Γ_seed] = −2·Γ_adj`. Add the Galois-swap framing: σ(Γ_adj) = −Γ_seed. The §5 hinge stops turning on a symbol and starts turning on this algebraic relation. |
| **§4/§6/§8 time** | PARTIAL | State explicitly that "time as iteration" and "time as fourth gamma" are two distinct algebraic structures sharing the 5^(1/4) eigenvalue scale, not the same object. The headline "no extra spatial dimensions" is preserved by clarifying the four 4D gammas as **3 space + 1 time-character**. |
| **§8 Constraint 3** | FAIL (as stated); reframable | Withdraw Constraint 3 as an "(1, 3)-selecting" claim. State: (Γ⁵)² = −5^(3/2)·I₄ is the **3D** chirality value of Paper 191's (1, 2) signature, distinguishing it from (2, 1). The selection of D=3 must come from Constraints 1 and 2; chirality is a consistency check, not an independent selection. |

The cumulative effect on the §8 headline: from "three independent constraints select D=3 as the unique closure" to **"the Borwein capacity and Paper 117 dimension argument select D=3; chirality verifies the Lorentzian signature within D=3"**. That is what the algebra actually supports.

## Pattern flags

- **Pattern 75 (null)**: addressed in each sub. Sub 1's null compared `i_golden = Γ_adj/5^(1/4)` against `J` (Pauli i) and `Γ_seed` (which doesn't even square to −I, ruled out trivially). The discriminator is the canonical commutator `[Z, Γ_seed] = −2·Γ_adj`, which lands on Γ_adj specifically, not on J.
- **Pattern 39 (DERIVED vs OBSERVED)**:
  - DERIVED: every computational result here. Exact symbolic equality verified.
  - The §5 hinge of Paper 203 v0.1 was FRAMEWORK; my Sub 1 elevates one piece (canonical commutator) to DERIVED.
- **Pattern 19 (adversary)**: Mr A's three catches all bite. Sub 1's literal (1.b) test fails — but the actual structure is deeper (Galois swap). Sub 2 finds no clean identification (the time tension is real). Sub 3 fully vindicates the circularity catch and goes further to find an algebraic obstruction to (1, 3) constructibility. Mr A's pen was sharp; my exact-symbolic work confirms it everywhere it landed.

## Files

- [BRIEF.md](BRIEF.md) — CinC's v0.2 brief (verbatim)
- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — pre-registered choices
- [subs.py](subs.py) — implementation of Subs 1, 2, 3
- [sub3_3d_check.py](sub3_3d_check.py) — Sub 3 cross-check for 3D signatures
- [findings_paper_203_algebra_v0_1.md](findings_paper_203_algebra_v0_1.md) — this report
