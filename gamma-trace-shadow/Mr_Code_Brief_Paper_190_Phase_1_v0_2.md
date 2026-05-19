# Mr Code Brief v0.2: γ as Trace-Shadow Residue — Phase 1 Reconnaissance

**Bootstrap Universe Programme**
**Role:** Mr Code (computational, scope-limited)
**Briefed by:** CinC
**Date:** May 2026
**Supersedes:** v0.1 brief (γ as harmonic-basis shadow of 1/√3)
**Companion to:** Paper 15 v3.3 retirement notice §4 (broader framework); Paper 184 (The Two Shadows, trace-shadow language)

---

## 1. Context

Paper 15 (Euler-Mascheroni constant from icosahedral geometry, v2.0 January 2026) is being retired. Paper 15 v3.3 §4 offers a structural reading of γ at OBSERVED structural status, naming a specific candidate (γ as shadow of 1/√3 in the wrong harmonic basis) and a broader framework (apparatus constants are residues of geometric operations) of which the specific candidate is one face.

This Phase 1 brief tests the **broader framework**, not the specific candidate. The broader framing — γ as the trace-shadow residue (Paper 184's language) of the additive/geometric mismatch — has corpus backing and a mathematically exact form via the Euler-Maclaurin expansion: γ is the constant term in H_n = ln n + γ + 1/(2n) − 1/(12n²) + 1/(120n⁴) + …, where the higher corrections are governed by Bernoulli numbers (which the programme already owns via ζ(−1) = −1/12 = B₂/2, ζ(−3) = 1/120, etc. per Paper 155).

The natural diagnostic tool is the **continued fraction representation**. CFs preserve multiplicative/Euclidean structure that decimal representations scramble. The pre-computer mathematical tradition (Brahmagupta, Bombelli, Brouncker, Euler, Lambert) read constants through CFs as a matter of course; the question this brief asks is whether γ-CF carries structural information that γ-decimal hides.

**Both outcomes useful.** Negative results are fully valuable — they let Paper 15 v3.3 §4 stand at OBSERVED structural permanent ceiling and tell the corpus that the broader trace-shadow framework has been honestly tested for γ. The structural decision (retire Paper 15) does not depend on Phase 1's outcome.

---

## 2. Pre-Registered Hypotheses

**STRONG hypothesis (H1):**
- One or more of {γ, γ_{ℚ(√5)}, γ_{ℚ(√−3)}, γ_{2I-spectral}} shows a structurally meaningful CF relationship to the others or to the natural icosahedral constants {1/√3, log φ, π/(3√3), 2log(φ)/√5}.
- Equivalently: γ's CF prefix contains structure predictable from Euler-Maclaurin/Bernoulli machinery beyond what a random Stieltjes-like constant of comparable irrationality measure would show.

**NULL hypothesis (H0):**
- γ-CF behaves indistinguishably from a random Stieltjes-like constant of comparable irrationality measure.
- The Stieltjes constants γ_{ℚ(√5)} and γ_{ℚ(√−3)} show no CF-level relationship to γ beyond chance.
- Bernoulli-derived prefix prediction for γ-CF (Task 4) does not exceed null distribution at p < 0.05.

**Pre-registration commitment (Pattern 75):** these hypotheses are committed before any computation. Mr Code git-commits this brief unchanged before running any code. The v0.1 brief and all prior CinC framing are explicitly out of scope for the computation; Mr Code has read them, the re-brief acknowledges that, and the hypothesis space above is the locked target. No silent merging with v0.1.

---

## 3. Falsifier Tolerances

**Necessary condition for H1 (broader framework to survive Phase 1):** at least one of —
- A statistically significant (p < 0.05 vs null distribution) CF-level relationship between γ and one of {γ_{ℚ(√5)}, γ_{ℚ(√−3)}, γ_{2I-spectral}}, OR
- Bernoulli-derived prefix prediction for γ-CF matching observed γ-CF prefix at p < 0.05 vs null distribution of random Stieltjes-like constants

**Sufficient condition for promotion to Phase 2 (worth writing Paper 190 around):** at least one of —
- Clear structural pattern (e.g., predictable partial-quotient blocks, visible periodicity, eigenvalue-like spacing) in γ-CF that maps to Bernoulli structure
- Two or more of the Stieltjes constants showing related CF structure at p < 0.001
- A clean prediction from the framework that distinguishes γ-CF from any random Stieltjes-like constant

**Anything weaker than necessary:** HALT-NULL. Paper 190 stays as v1.1 (the φ paper). Paper 15 v3.3 §4 stands as drafted at OBSERVED structural permanent ceiling. The broader trace-shadow framework remains valid for the corpus (Paper 184); it just does not specifically predict γ structure beyond chance.

---

## 4. Tasks (run in order; git-commit before each)

### Task 1: Bernoulli baseline (verification)

**Goal:** Verify the standard Euler-Maclaurin machinery to high precision before any speculative computation.

**Compute:**
- γ to 50 digits via H_n − ln n − Σ_{k=1}^{K} B_{2k}/(2k·n^{2k}) for n = 10⁶, K = 10. Use mpmath or equivalent for arbitrary precision.
- Report the first 10 terms of the asymptotic expansion explicitly with exact Bernoulli coefficients.
- Confirm convergence to known γ = 0.5772156649015328606065120900824024310421593359399...

**Output:** Standard verification table.

### Task 2: Continued-fraction table

**Goal:** Compute and tabulate CF expansions of the candidate set to a few hundred partial quotients each, side by side.

**Compute CFs to N = 300 partial quotients (or as far as precision permits) for each of:**

a. **γ** (Euler-Mascheroni) — reference object
b. **1/√3** — the §4 specific-candidate quantity
c. **log φ** — natural icosahedral logarithm
d. **π/(3√3)** — L(1, χ_{−3}), Eisenstein residue
e. **2 log(φ)/√5** — L(1, χ₅), Dedekind residue for ℚ(√5)
f. **γ_{ℚ(√5)}** — Stieltjes constant of the Dedekind ζ-function for ℚ(√5)
g. **γ_{ℚ(√−3)}** — Stieltjes constant for ℚ(√−3) (Eisenstein integers ℤ[ρ])
h. **γ_{2I-spectral}** — *stretch target*; the Stieltjes constant of S_n = ∑_{k ≤ n, k ∈ S(12,20,30)} m_k / λ_k where m_k is the Molien multiplicity and λ_k = k(k+2) (per Paper 174's icosahedral spectrum). Compute only if the sum is shown to diverge in the right way to admit a Stieltjes constant; if not, document the divergence/convergence behaviour and skip h.

**Note on (f) and (g) computation:** the Stieltjes constants of Dedekind ζ-functions are well-defined but require care. For ℚ(√5), γ_{ℚ(√5)} appears in ζ_{ℚ(√5)}(s) = (Res_{s=1})/(s−1) + γ_{ℚ(√5)} + O(s−1), where Res_{s=1} = 2 log(φ)/√5 (Dirichlet class number formula). For ℚ(√−3), Res_{s=1} = π/(3√3). Standard computation via the functional equation and the Riemann-Siegel-like expansion; report which method is used.

**Output:** Single table with eight rows (or seven if h is skipped), columns being partial quotients 1 through N, plus first 10 convergents shown explicitly for each. Flag any visible regularities: periodic blocks, repeated patterns, near-coincidences between rows.

### Task 3: Side-by-side relationship analysis

**Goal:** Identify any statistically significant relationships between the CF representations of {γ, γ_{ℚ(√5)}, γ_{ℚ(√−3)}, γ_{2I-spectral}}.

**Compute:**

a. **Partial-quotient correlation:** for each pair, compute Pearson correlation of partial-quotient sequences. Bonferroni-correct for the number of pairs tested.

b. **Block-pattern detection:** scan each CF for repeated blocks of length 2 through 10; report any blocks appearing more frequently than expected under the Khinchin-Lévy null hypothesis (which gives the expected distribution of CF partial quotients for "almost all" irrationals).

c. **Cross-CF block matching:** identify any blocks of length ≥ 4 appearing in two different CFs in the candidate set; compute the probability of such matching under the null.

d. **Convergent comparison:** for each CF, list the first 10 convergents and their denominators; check whether the denominators of the candidate set share common factors or relationships (icosahedral integers, Fibonacci/Lucas numbers, etc.).

**Output:** Correlation matrix, block table, convergent table. Flag any p-value below 0.05 (uncorrected) for follow-up; flag any below 0.05 after Bonferroni correction as significant.

### Task 4: Euler-Maclaurin prefix prediction test

**Goal:** Test whether the asymptotic expansion γ = H_n − ln n − Σ_{k=1}^{K} B_{2k}/(2k·n^{2k}) + R(n,K) implies any predictable structure in γ's CF prefix, beyond what is implied by the known decimal value of γ.

**Approach:**

a. Take the truncated approximation γ_K(n) = H_n − ln n − Σ_{k=1}^{K} B_{2k}/(2k·n^{2k}) for various (n, K). Each γ_K(n) is a rational number (since H_n and B_{2k} are rational) approximating γ.

b. Compute the CF of γ_K(n) for several (n, K) and compare to the CF of γ. The prefix should match up to a precision determined by R(n, K).

c. **Key question for H1:** does the *structure* of γ_K(n)'s CF (its partial-quotient pattern, its convergent denominators) carry information about γ's CF beyond what is implied by γ_K(n)'s decimal value? Specifically: do the CF prefixes of γ_K(n) for different (n, K) share structural features (e.g., common partial-quotient subsequences) that would not be expected for arbitrary rational approximations of γ at comparable precision?

d. **Null comparison:** compare to CFs of γ_K(n) when γ is replaced by random Stieltjes-like constants of comparable irrationality measure (Task 5 provides these). The prediction test passes if γ shows structure that random Stieltjes-like constants don't.

**Output:** Table of γ_K(n) CFs for various (n, K), comparison to γ-CF prefix, and null comparison.

### Task 5: Null distribution

**Goal:** Generate the null distribution against which Tasks 3 and 4 are evaluated.

**Generate:** at least 1000 random "Stieltjes-like constants" — random real numbers in (0, 1) with continued fractions sampled from the Gauss-Kuzmin distribution (the natural prior for CF partial quotients of generic irrationals). Compute the same statistics as Tasks 3 and 4 (partial-quotient correlations, block patterns, convergent structures) for each, generating an empirical null distribution for each statistic.

**Output:** Null distributions plotted or tabulated; thresholds at p = 0.05 and p = 0.001 explicit.

---

## 5. Stop Conditions

**Continue to Phase 2 (Paper 190 rewrite around γ-as-trace-shadow with empirical legs) iff:**
- Task 3 or Task 4 shows at least one finding meeting the necessary condition (p < 0.05 against null), AND
- The finding is structurally interpretable (a partial-quotient block that maps to a recognisable icosahedral or Bernoulli structure, not just a numerical correlation)

**Otherwise halt:**
- **HALT-NULL:** all tests at p ≥ 0.05 against null. The broader trace-shadow framework does not specifically predict γ structure beyond chance. Paper 15 v3.3 §4 stands at OBSERVED structural permanent ceiling; Paper 190 remains v1.1 (φ paper); recommend a §5 sharpening in v1.1 to acknowledge the negative result without overclaiming the φ content.

- **HALT-PARTIAL:** Some tests at p < 0.05 but not structurally interpretable (numerical correlation without recognisable mechanism). Document carefully; return for re-briefing on whether the partial finding is worth following up or is consistent with multiple-testing artifact.

- **HALT-AMBIGUOUS:** The null distribution itself is uncertain (e.g., the Khinchin-Lévy null is not a good model for Stieltjes constants because Stieltjes constants are not "generic" irrationals). Report and return for re-briefing on a better null.

---

## 6. Output Format

Deliver a single markdown report `Paper_190_Phase_1_v0_2_recon.md` containing:

1. Pre-registration record (Hypotheses H0/H1 as committed at start, with git hash of this brief at the time of commit; explicit acknowledgement that v0.1 and CinC framing were read prior to the re-brief, with the hypothesis space of v0.2 above being the locked target)
2. Task 1: Bernoulli baseline table
3. Task 2: Continued-fraction table (the main computation)
4. Task 3: Relationship-analysis results
5. Task 4: Euler-Maclaurin prefix prediction results
6. Task 5: Null distribution summary
7. Stop condition outcome: CONTINUE / HALT-NULL / HALT-PARTIAL / HALT-AMBIGUOUS
8. If CONTINUE: which specific findings the Phase 2 work should develop, with explicit structural interpretation
9. If HALT: which sub-condition, and what (if anything) is recommended next

No interpretation beyond what's needed to fill the table. No promotion of OBSERVED to DERIVED. No structural claims beyond the existence/non-existence statements the tasks are designed to test.

---

## 7. Discipline Reminders

- **Git commit this brief before running any code.** No edits to this brief after commit; if scope changes, halt and request a v0.3 brief explicitly.
- **Pre-registration acknowledgement.** This brief is v0.2 because v0.1 was compromised by surrounding CinC framing being readable. The hypothesis space above is the locked target; v0.1's specific 1/√3-shadow test is **explicitly shelved**, not folded into v0.2. If Phase 2 reaches the 1/√3 question, it will be a separate sub-brief.
- **Report all results, not just positive ones.** A clean HALT-NULL is fully valuable and what Paper 15 v3.3 §4 has explicitly budgeted for.
- **Pattern 75 throughout.** Every numerical match claim reports its null distribution. The whole point of Phase 1 is to test whether agreement is structural or generic.
- **Khinchin-Lévy null caveat.** If the null distribution model proves inadequate (Stieltjes constants are not generic irrationals), halt and request a better-null brief rather than proceeding with a known-inadequate null.

---

*Mr Code brief v0.2, scope-limited per Pattern 97. Phase 1 reconnaissance — both outcomes useful. v0.1 (1/√3-specific) shelved; this v0.2 (broader trace-shadow framework) is the locked target.*

🐕☕⬡
