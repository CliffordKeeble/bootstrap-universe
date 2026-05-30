# Mr Code triple brief — Subs A, B, and C-extension

**Date**: 30 May 2026
**Brief author**: CinC
**Context**: Three independent investigations to close specific catches before Paper 150 v2.1 and Paper 203 v0.4 draft. Each sub is bounded, runnable independently of the others, and answers a specific pre-registered question with discovery- or null-threshold logic. Total estimated compute: ~one day across all three.

**Reads required for ALL subs** (read once at the start):
- This brief (BRIEF.md)
- The findings from Sub C (`findings_paper_203_sub_c.md`) — for context on what Sub C established and what these sub-investigations extend
- Paper 191 v1.2 (Golden Dirac Algebra) — for Sub A
- Paper 150 v2.0 — for Sub C-extension (the existing machinery to extend)
- Paper 164 v1.2 (Dedekind Bridge, ℚ(ρ)) §2 — for Sub C-extension non-quadratic test

**Reads NOT required** (Pattern 97 — scope protection):
- Mr Adversary's reviews of Papers 150 v2.0 and 203 v0.3 — CinC's summaries below are sufficient
- Paper 203 v0.3 — the framework reading is interpretive; Mr Code stays at the algebraic / computational level

---

## Sub A — Signature enumeration scope extension

### Question

Does the (2,2) → (1,3) signature obstruction found in Sub 3 (paper-203-algebra, commit b19747d) hold in the full ℤ[φ]-linear span in M₄(ℚ(√5)), or is it an artefact of the literal tensor-product basis enumeration?

### Background

Sub 3 enumerated 36 literal tensor products A⊗B with A, B ∈ {I, K, J, Z, Γ_seed, Γ_adj}, found 16 with squares = ±√5·I, and identified 6 anticommuting cliques of 4 elements — ALL with signature (2,2) Klein, NONE with signature (1,3) Minkowski. The (Γ⁵_(4))² = +25·I₄ in all 6 cliques.

Mr Adversary's catch (v0.2 of Paper 203): the literal tensor-product basis is a 36-dimensional ℤ-module inside the 16-dimensional ℚ(√5)-algebra M₄(ℚ(√5)), which has 32 ℤ[φ] dimensions as a ℤ-module. The full span is larger and the (1,3) clique might exist there.

### Test design

**Method 1 — bounded coefficient search within "anticommutes-with-Γ⁰" subspace.**

Step 1. Fix Γ⁰ = Z ⊗ Γ_seed (the unique literal-tensor-product candidate with (Γ⁰)² = +√5·I_4 from Sub 3).

Step 2. Define the ℤ[φ]-linear span of the 36 tensor products:

  V = span_{ℤ[φ]}(A ⊗ B : A, B ∈ {I, K, J, Z, Γ_seed, Γ_adj})

Step 3. Compute the subspace W ⊂ V of matrices M satisfying:
  (a) {M, Γ⁰} = 0 (anticommutes with Γ⁰)
  (b) M² = ±√5 · I_4 (its square is ±√5 times identity)

Step 4. Enumerate all M ∈ W with ℤ[φ] coefficients in {-1, 0, 1}·{1, φ} (4 values per coefficient, bounded search). Then extend to {-2, -1, 0, 1, 2}·{1, φ} (10 values per coefficient) if no clique found.

Step 5. Search for triples (Γ¹, Γ², Γ³) ⊂ W (mutually anticommuting) such that exactly 0 of the three have (Γⁱ)² = +√5·I (the (1,3) signature: Γ⁰ alone positive).

**Method 2 — algebraic-classification check.**

In parallel, check whether Cl(1,3) over ℚ(√5) is representable as a sub-algebra of M₄(ℚ(√5)) with ℤ[φ] coefficients. Standard result: Cl(1,3) ⊗_ℝ ℂ ≅ M₄(ℂ); over ℚ(√5), Cl(1,3) has signature-specific structure that may or may not embed in M₄(ℚ(√5)) with rational-integer × φ coefficients.

If Method 2 gives a clean theoretical answer (yes/no representable), report that; the Method 1 search either confirms (finds the embedding) or is consistent with the negative theoretical answer.

### Pre-registered thresholds

- **CONSTRUCTIBLE**: Method 1 finds a (Γ¹, Γ², Γ³) triple with (Γⁱ)² = −√5·I for all i, anticommuting with each other and Γ⁰. The four matrices form a (1,3) Minkowski clique in the ℤ[φ]-linear span.
- **NOT CONSTRUCTIBLE at bound K**: Method 1 exhausts the search at coefficient bound K (= 1, 2, ...) without finding a clique. Report the largest K tested.
- **OBSTRUCTION-THEORETIC**: Method 2 establishes a representation-theoretic obstruction: Cl(1,3) over ℚ(√5) cannot be realised in the ℤ[φ]-linear span of M₄(ℚ(√5)). Strongest possible negative result.

### Stop-on-find

If Method 1 finds a CONSTRUCTIBLE clique at any bound K, stop the search and report. Mr A's v0.2 catch is then sustained: the literal tensor-product enumeration was not exhaustive; (1,3) Minkowski IS constructible in the ℤ[φ] span. Paper 203 v0.4 §8 can claim the "(2,2)→(1,3) candidate baryogenesis paper" framing with this finding as the scope expansion.

### Stop-on-exhaust

If Method 1 exhausts to K=2 without finding a clique AND Method 2 yields an obstruction theorem, stop and report the strongest result. NOT CONSTRUCTIBLE in the ℤ[φ] span — Paper 203 v0.4 §8's obstruction claim is now scope-honest at the level of the full ℤ[φ]-linear span, with a representation-theoretic backing.

### Implementation notes

- Pure exact symbolic Sympy. No floating-point.
- The ℤ[φ] coefficients can be represented as (a, b) pairs with the multiplication rule (a, b) * (c, d) = (ac + bd, ad + bc + bd) — the standard ℤ[φ] ring arithmetic where φ² = φ + 1.
- The anticommutator {M, N} = MN + NM is a polynomial in the entries; for 4×4 matrices with ℤ[φ] entries, this is exact-symbolic and fast per pair.
- Estimated compute: K=1 search is ~minutes; K=2 search is ~hours; K=3 search may be impractical without algebraic pruning.

### Sub A deliverable

A markdown report (`findings_paper_203_sub_a.md`) with:
- Method, search bound K reached, total clique candidates enumerated
- For the largest K tested: number of valid (Γ¹, Γ², Γ³) triples with each signature class
- Verdict: CONSTRUCTIBLE (with explicit clique) / NOT CONSTRUCTIBLE at K / OBSTRUCTION-THEORETIC
- Implications for Paper 203 v0.4 §8

---

## Sub B — Zeta-value-denominator chance overlap null

### Question

Are the matches between zeta-value denominators and icosahedral subgroup orders (12, 120, 240) striking beyond chance?

### Background

Paper 203 §7 currently notes three matches:
- denom(ζ(−1)) = 12, and the icosahedron has 12 vertices
- denom(ζ(−3)) = 120, and |2I| = 120
- denom(ζ(−7)) = 240 = 2·|2I| = #(E₈ roots)

Mr A's catch: von Staudt–Clausen (1840) gives a complete elementary account of zeta-value denominators (denom(B_{2k}) is the product of primes p with (p−1) | 2k). The matches above are consistent with elementary number theory and require a null comparison before being called striking.

### Test design

**Step 1 — Define the populations.**

Population Z: denominators of ζ(−1), ζ(−3), ζ(−5), ..., ζ(−2k+1) for k = 1, ..., N. Choose N large enough to make the chance-overlap measurement meaningful — N = 50 (so we have ζ(−1) through ζ(−99)) is fine.

Population I: orders of subgroups of the binary icosahedral group 2I (order 120) and the icosahedral group A₅ (order 60). These are determined by Lagrange's theorem to be divisors of 120 (for 2I) and 60 (for A₅), plus doubles (2·|2I| = 240 etc.):
  {1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, 120, 240}

Also include: structurally meaningful counts in the icosahedral / 2I / E₈ realm:
  - 30 = number of edges of the icosahedron
  - 20 = number of faces of the icosahedron
  - 12 = number of vertices of the icosahedron
  - 30 = number of edges = Coxeter number of E₈
  - 240 = number of E₈ roots
  - 120 = number of order-60 vertices of the 600-cell (= |2I|)
  - 168 = spectral gap λ₁ of S³/2I (rendering it a "meaningful icosahedral number")
  - 60 = |A₅|

So Population I (icosahedral, structurally meaningful) = {1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, 120, 168, 240} (15 elements).

**Step 2 — Compute the actual match rate.**

Of the 50 zeta-value denominators in Z, how many match an element of Population I exactly? Compute.

**Step 3 — Compute the expected chance match rate.**

The probability that a uniformly random integer drawn from the natural density of "small integers that are products of primes ≤ p_max" lands in Population I is some computable density. von Staudt-Clausen gives the structure of Z explicitly: denom(B_{2k}) = ∏_{p prime, (p-1) | 2k} p. Compute this for k = 1, ..., 50 and tabulate. Then compute the empirical density of Population I within the same range of integer values.

Specifically: let M = max(Z ∪ I). Within {1, ..., M}, the density of Population I is |I|/M. If Z and I are independent, the expected match rate is approximately |Z| · |I|/M.

**Step 4 — Compute the z-score.**

Compare actual vs expected. Under the null hypothesis of independent populations, the number of matches is approximately Poisson with mean μ = |Z| · |I|/M.

  z = (actual − μ) / √μ

If |z| ≥ 3 (or some pre-registered threshold), matches are striking. If |z| < 2, matches are consistent with chance.

### Pre-registered thresholds

- **STRIKING**: |z| ≥ 3. Matches are 3+ standard deviations above chance. Paper 203 v0.4 §7 can claim "consistent with framework and beyond chance."
- **CONSISTENT WITH FRAMEWORK ONLY**: 2 ≤ |z| < 3. Matches are non-trivial but not striking. v0.4 §7 says "consistent with framework, marginally above chance."
- **CONSISTENT WITH CHANCE**: |z| < 2. Matches are not statistically distinguishable from coincidence. v0.4 §7 must demote the matches to "consistent with chance under von Staudt-Clausen" and the §7 phenomenon is purely elementary.

### Stop-on-fail

If the actual match rate is BELOW the expected rate (i.e., fewer matches than chance), report and stop. This would indicate a structural anti-correlation, which is even more interesting than a positive correlation but requires separate analysis. Flag for CinC.

### Implementation notes

- Pure rational arithmetic; Sympy can compute Bernoulli numbers exactly.
- Estimated compute: minutes.
- Need to be careful about WHICH "icosahedral" integers count. Pre-register the population I list before running the test.

### Sub B deliverable

A markdown report (`findings_paper_203_sub_b.md`) with:
- Populations Z (50 zeta denominators) and I (icosahedral meaningful integers), exactly tabulated
- Actual match count
- Expected match count under Poisson null
- z-score
- Verdict: STRIKING / CONSISTENT WITH FRAMEWORK ONLY / CONSISTENT WITH CHANCE
- Implications for Paper 203 v0.4 §7

---

## Sub C-extension — positive-definite control + non-real-quadratic L-function test

### Question

Two independent mechanism tests:
- **Sub C-ext-1 (positive-definite control)**: Does the indefinite-norm requirement hold for non-unit coefficients? Test |Re² + 5·Im²| as positive-definite probe.
- **Sub C-ext-2 (non-real-quadratic L-function)**: Does the field-general detection extend to L-functions of non-real-quadratic characters? Test |Re² − 5·Im²| against zeros of L(s, χ₂ mod 7) (a cubic character from Paper 164).

### Background

Sub C-1 found |Re² − d·Im²| with d ∈ {3, 5, 7, 13} all detect zeros of L(χ_d') for d' ∈ {3, 5, 7, 13} at uniform effect ε ≈ 0.29 (matrix verdict). Mr Code's mechanism analysis identifies the equidistributed Dirichlet sum as the carrier and the discriminant d as "decoration." Two further tests pin down the actual scope of the mechanism:

**For Sub C-ext-1**: Paper 150 v2.0 Test 3a tested |Re² + Im²| (circular, coefficient 1) and got z = −2.63 (fails). The remaining ambiguity is whether the failure was about coefficient 1 specifically or about positive-definiteness in general. Predicted: the positive-definite |Re² + 5·Im²| also fails. Confirming this lets v2.1 claim "indefiniteness is essential" on direct evidence rather than by inference from Paper 150 Test 3a + Sub C-1.

**For Sub C-ext-2**: Paper 164 establishes ζ_ℚ(ρ)(s) = ζ(s)·L(s, χ₂)·L(s, χ̄₂) where χ₂ is a cubic character mod 7. The L-function L(s, χ₂) is on the critical line and has zeros there. Whether the indefinite-norm Dirichlet-sum probe also detects these cubic-character L-function zeros is open. If yes, the detection mechanism is generic to any L-function with FE on the critical line; if no, the mechanism is specifically real-quadratic in scope.

### Test design for Sub C-ext-1 (positive-definite control)

Use the same probe machinery as Sub C-1 but with the indefinite norm replaced by a positive-definite norm:

- **Probe**: ζ_φ,5+(t) = |Re²(Z_φ(t)) + 5 · Im²(Z_φ(t))|

This is positive-definite: the value is always non-negative, with zero only at Re = Im = 0 (the origin in 2D pair space). No null cone.

- **Target**: zeros of L(s, χ_5) from Sub C-1's cache (904 zeros, t ∈ (1, 1000]).

- **Run**: same machinery as Sub C-1's (5, 5) cell; classification window W = 1.0, N_null = 1000, seed = 42, N_terms = 5000.

- **Compare**: z-score and effect size ε against Sub C-1's (5, 5) cell which had z = −8.90, ε = 0.298.

### Pre-registered thresholds for Sub C-ext-1

- **CONFIRMS-INDEFINITENESS-ESSENTIAL**: |z| < 3 (well below discovery threshold). The positive-definite |Re² + 5·Im²| does NOT detect L(χ_5) zeros; indefiniteness is the essential ingredient.
- **REFUTES-INDEFINITENESS-CLAIM**: |z| ≥ 5 (discovery threshold). The positive-definite norm DOES detect; indefiniteness is not essential. Mechanism is even more generic than Mr Code's analysis suggested.
- **AMBIGUOUS**: 3 ≤ |z| < 5. Marginal; rerun at N = 10⁴ (using cached Riemann zeros if available).

### Test design for Sub C-ext-2 (cubic-character L-function)

- **Probe**: ζ_φ,5(t) = |Re²(Z_φ(t)) − 5 · Im²(Z_φ(t))| — the canonical Paper 150 v2.0 probe.

- **Target**: zeros of L(s, χ₂ mod 7) — the cubic character of Paper 164 with χ₂(a) = ω^k for k determined by a's discrete log mod 7 (where ω = e^(2πi/3)).

Specifically:
  - a = 1: χ₂(1) = 1
  - a = 2: χ₂(2) = ω
  - a = 3: χ₂(3) = ω²  (3 is a primitive root mod 7)
  - a = 4: χ₂(4) = ω²
  - a = 5: χ₂(5) = ω
  - a = 6: χ₂(6) = 1
  - a = 0: χ₂(0) = 0 (trivial extension)

Note: χ₂ is COMPLEX-VALUED, not real. The L-function L(s, χ₂) has zeros on the critical line under GRH. Its functional equation relates L(s, χ₂) to L(1−s, χ̄₂) — so the FE is not self-conjugate.

Implementation note: the L-function evaluation for complex characters requires a slightly different Hardy phase than for real characters. Mr Code can either:
  (a) Use mpmath's L-function library if it handles complex characters directly
  (b) Implement the Hardy phase for complex characters: θ(t, χ) = arg Γ_χ(½ + it) − (t/2) log(N/π), where N is the conductor and Γ_χ depends on whether χ is even or odd. For χ₂ mod 7 with χ₂(−1) = χ₂(6) = 1, the character is EVEN, so Γ_χ(s) = Γ(s/2).

Compute zeros of L(s, χ₂ mod 7) at t ∈ (1, 500] (matching the Sub C-1 d=7 target range). Pre-register the number of zeros found.

- **Run**: classification machinery from Sub C-1 against the new cubic-character zeros.

### Pre-registered thresholds for Sub C-ext-2

- **DETECTED (mechanism is generic to L-functions)**: |z| ≥ 5. The indefinite-norm probe detects cubic-character L-function zeros at discovery threshold. The mechanism is generic to any L-function with FE on critical line, not specific to real-quadratic L-functions. Paper 150 v2.1 title says "L-function zeros" generally.
- **NOT DETECTED (mechanism is real-quadratic-specific)**: |z| < 3. The probe does NOT detect cubic-character zeros. There's a scope boundary at "real-quadratic L-functions." Paper 150 v2.1 title says "real-quadratic L-function zeros" specifically.
- **AMBIGUOUS**: 3 ≤ |z| < 5. Marginal; rerun at higher N or higher t-range; flag for CinC.

### Stop-on-fail

Sub C-ext-1 and Sub C-ext-2 are independent. Either can be run first or both in parallel.

If Sub C-ext-1 returns REFUTES-INDEFINITENESS-CLAIM (positive-definite also detects), this is the more surprising outcome and demands immediate investigation: what IS the essential ingredient if not indefiniteness? Pause Sub C-ext-2 and flag for CinC.

If Sub C-ext-2 returns DETECTED for a cubic L-function, this is mechanistically interesting (the mechanism extends beyond real-quadratic) but doesn't require pausing further work.

### Implementation notes

- Sub C-ext-1 reuses Sub C-1's infrastructure entirely; only the probe norm changes. ~30-60 minutes compute.
- Sub C-ext-2 requires a new L-function zero finder for cubic complex characters. ~2-3 hours including validation against any LMFDB data for L(s, χ₂ mod 7).
- Both can run in parallel.

### Sub C-extension deliverable

A markdown report (`findings_paper_203_sub_c_extension.md`) with two sections (one per sub-investigation):

**Section 1 (Sub C-ext-1, positive-definite control):**
- Construction: probe form, N values, target zeros
- Result: z-score, effect size ε
- Verdict: CONFIRMS-INDEFINITENESS-ESSENTIAL / REFUTES-INDEFINITENESS-CLAIM / AMBIGUOUS

**Section 2 (Sub C-ext-2, cubic L-function):**
- Construction: cubic character verification, L-function zero finder validation, N values
- Result: z-score, effect size ε
- Verdict: DETECTED / NOT DETECTED / AMBIGUOUS

**Combined synthesis section:**
- Implications for Paper 150 v2.1 title and mechanism section
- Implications for Paper 203 v0.4 §5/§7 framing
- Open questions surfaced

---

## Hygiene (all three subs)

- **Pre-registration is binding per sub**. Commit BRIEF.md + per-sub PRE_REGISTRATION_{A,B,C-ext}.md before:
  - Sub A: any matrix enumeration begins
  - Sub B: the populations Z and I are finalised in tabulated form
  - Sub C-ext: any zero data is pulled for cubic character; any probe is run for positive-definite test
- **Pattern 97 (scope protection)**: Mr Code does NOT read Mr A's reviews. The CinC summaries in this brief are sufficient.
- **No post-hoc rule changes** per sub. If a result is surprising, report and let CinC adjudicate before extending.
- **Stop-on-fail respected** per sub.
- **Anti-circularity** verified per sub:
  - Sub A: exact symbolic; no data dependency
  - Sub B: zeta denominators are von Staudt-Clausen-computed exactly; icosahedral integers are group-theoretic
  - Sub C-ext: pipeline reproduction implicit via reuse of Sub C-1 machinery
- **Honest reporting**: a clean refutation of any pre-registered threshold is more valuable than a forced confirmation. Report what comes out.

---

## Deliverable summary

Three markdown reports:
- `findings_paper_203_sub_a.md` (signature enumeration scope extension)
- `findings_paper_203_sub_b.md` (zeta-denominator chance overlap null)
- `findings_paper_203_sub_c_extension.md` (positive-definite + cubic L-function)

Plus the per-sub pre-registration documents and any supporting computation files (analogous to Sub C's chars.py, l_zeros.py, etc.).

---

## Compute budget

- Sub A: half a day if K=2 search needed; faster if Method 2 obstruction theorem found quickly
- Sub B: under an hour
- Sub C-ext: half a day combined (Sub C-ext-1 quick; Sub C-ext-2 needs cubic-character zero finder)

**Total: ~one day of compute.**

---

## Order

Independent — run in any order or parallel. Recommended priority by impact on paper drafting:

1. **Sub B first** (quickest; result determines whether Paper 203 v0.4 §7 retains or demotes the matches).
2. **Sub C-ext in parallel with Sub A** (both substantive but independent).
3. **Sub A last** if K=2 search is long; Paper 203 v0.4 §8 can hedge until Sub A reports.

---

⌨️ over to you.

Three bounded investigations, each with discovery- or null-threshold logic. Each answers a specific pre-registered question. None depends on the others. Total compute ~one day. Outputs feed directly into Paper 150 v2.1 and Paper 203 v0.4 drafting.

🐕☕⬡

---

## Mr Code's reports

*[To be filled in per sub.]*
