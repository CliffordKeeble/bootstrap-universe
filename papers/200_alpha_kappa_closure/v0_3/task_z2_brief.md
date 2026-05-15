# Task Z' — Within-Regime Null Test for Paper 200 v0.3 (corrected brief)

**Brief for Mr Code · CinC handoff · 15 May 2026**
**Standing Trust: GREEN**

---

## CinC error acknowledgment

Task Z (commit `61aff82`, halted on anchor) ran the brief I wrote. The brief specified the wrong convention: I asked Mr Code to compute "running sum of χ₅ over all primes" when Paper 199 §2 defines the 5-series running sum as **over primes ≡ 5 (mod 6) only** (residue-class-filtered, not all-primes). Under the correct convention the anchor S₅(31) = 0 reproduces, and Paper 199 §6.1's "14 ties / 33 primes" figure reproduces (hand-verified by CinC after Mr Code's halt; verification in §"Anchor verification" below).

Task Z' supersedes Task Z under the corrected convention. Task Z's artefacts are preserved in place (`task_brief.md`, `task_z_within_regime_null.py`, `task_z_results.md`, `task_z_results.json`, `findings.md`) as the failed-brief audit record. Task Z' artefacts use `task_z2_*` naming and live in the same `papers/200_alpha_kappa_closure/v0_3/` folder.

Mr Code's halt was textbook execution of the pre-registration discipline: the locked script halted at method step 4 when the anchor didn't reproduce under the brief's specification, exactly as the brief instructed. The pre-registration discipline protected the programme from a wrong-brief outcome. No criticism of Mr Code — this is CinC's error in writing the original brief.

---

## Background (recap from Task Z, unchanged in substance)

Paper 200 v0.2.2 was reviewed by Mr Adversary v2.3 round-1 (three stars, "I would publish this … send it back with the null test attached"). Mr A's principal concern: Paper 199 §6.1 documents ~42% tie density in the 5-series for n ≤ 33; the probability under the null hypothesis that S₅(31) = 0 for any reason is therefore on the order of 0.42, which is barely discriminating. The four-crossing conjunctive design's discriminating power is borne almost entirely by Crossing 1, and Crossing 1 alone is coin-flip-adjacent against the regime-density alternative.

Task Z' confirms (or refutes) Paper 199 §6.1's documented tie-density figure under direct computation with explicit pre-registered decision criteria. The result feeds Paper 200 v0.3's abstract and §1 load-bearing-claim framing.

---

## Convention (corrected — Paper 199 §2 definition, explicit)

The 5-series is the sequence of primes ≡ 5 (mod 6): {5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89, 101, ...}. Note that 5 itself is included as the first 5-series prime (it is ≡ 5 mod 6 trivially: 5 mod 6 = 5).

For each n ≥ 1, let p_{5,n} denote the n-th prime ≡ 5 (mod 6). The 5-series running sum is:

```
S₅(n) = Σ_{j=1}^{n} χ₅(p_{5,j})
```

where χ₅ is the Legendre symbol mod 5:
- χ₅(p) = +1 if p ≡ 1 or 4 (mod 5), equivalently if p splits in ℤ[φ]
- χ₅(p) = −1 if p ≡ 2 or 3 (mod 5), equivalently if p is inert in ℤ[φ]
- χ₅(5) = 0 (ramified — the only 5-series prime in this case)

This is **Definition A** of Paper 199 §6.1: include the ramified step at n = 1 (p = 5 contributes 0 to the sum). Under this convention S₅(1) = 0 is itself counted as a tie. Paper 199 §6.1 lists tie positions in n ∈ [1, 33] as: **{1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}** — 14 ties.

(The brief is explicit about Definition A because Paper 199 also discusses Definitions B and C which differ in bookkeeping of the ramified step. Definition A is what Paper 200's load-bearing claim depends on and what this test uses. The dependence is a structural choice — the 5-series's natural baseline at zero — and Paper 200 v0.3 will flag this explicitly per Mr A's "CONVENTION minor" point.)

---

## Anchor verification (CinC hand-computation, for Mr Code's sanity check)

The first 31 primes ≡ 5 (mod 6) and their χ₅ values:

| n | p_{5,n} | p mod 5 | χ₅ | S₅(n) | Tie? |
|---|---|---|---|---|---|
| 1 | 5 | 0 | 0 | 0 | ✓ |
| 2 | 11 | 1 | +1 | 1 | |
| 3 | 17 | 2 | −1 | 0 | ✓ |
| 4 | 23 | 3 | −1 | −1 | |
| 5 | 29 | 4 | +1 | 0 | ✓ |
| 6 | 41 | 1 | +1 | 1 | |
| 7 | 47 | 2 | −1 | 0 | ✓ |
| 8 | 53 | 3 | −1 | −1 | |
| 9 | 59 | 4 | +1 | 0 | ✓ |
| 10 | 71 | 1 | +1 | 1 | |
| 11 | 83 | 3 | −1 | 0 | ✓ |
| 12 | 89 | 4 | +1 | 1 | |
| 13 | 101 | 1 | +1 | 2 | |
| 14 | 107 | 2 | −1 | 1 | |
| 15 | 113 | 3 | −1 | 0 | ✓ |
| 16 | 131 | 1 | +1 | 1 | |
| 17 | 137 | 2 | −1 | 0 | ✓ |
| 18 | 149 | 4 | +1 | 1 | |
| 19 | 167 | 2 | −1 | 0 | ✓ |
| 20 | 173 | 3 | −1 | −1 | |
| 21 | 179 | 4 | +1 | 0 | ✓ |
| 22 | 191 | 1 | +1 | 1 | |
| 23 | 197 | 2 | −1 | 0 | ✓ |
| 24 | 227 | 2 | −1 | −1 | |
| 25 | 233 | 3 | −1 | −2 | |
| 26 | 239 | 4 | +1 | −1 | |
| 27 | 251 | 1 | +1 | 0 | ✓ |
| 28 | 257 | 2 | −1 | −1 | |
| 29 | 263 | 3 | −1 | −2 | |
| 30 | 269 | 4 | +1 | −1 | |
| 31 | 281 | 1 | +1 | **0** | **✓** ★ ANCHOR |

13 ties in n ∈ [1, 31]. The 14th tie is at n = 33: p_{5,32} = 293 (3 mod 5, χ₅ = −1, S₅ = −1); p_{5,33} = 311 (1 mod 5, χ₅ = +1, S₅ = 0 ✓).

Tie positions in n ∈ [1, 33]: **{1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}** — 14 ties, matching Paper 199 §6.1 exactly.

If Mr Code's implementation reproduces this table verbatim (including the χ₅(5) = 0 ramified entry, the +1/−1 pattern matching p mod 5, and S₅(31) = 0), the anchor verification passes and the regime-fraction computation proceeds. If any row differs, halt as in Task Z.

---

## Pre-registered decision criteria (unchanged from Task Z)

Let `f_tight` = fraction of n ∈ [1, 33] with S₅(n) = 0. Computed directly: 14 / 33 = 0.4242…

Decision bins:
- **Bin LOW**: `f_tight < 0.30`. Structural reading materially stronger than null. v0.3 keeps v0.2.2 framing largely intact.
- **Bin MID**: `f_tight ∈ [0.30, 0.55]`. Regime-density null consistent with structural reading at roughly even odds. v0.3 abstract+§1 load-bearing-claim block must be rewritten per Mr A.
- **Bin HIGH**: `f_tight > 0.55`. Regime-density null stronger than observation. v0.3 scope-limits substantially.

CinC's prior (per hand-computation above and Paper 199 §6.1): **Bin MID, very likely f_tight = 14/33 = 0.4242…** to within rounding.

The bins are binding once the brief commits. Do not adjust mid-run.

**Random seed: 200** (paper number; consistent with Task Z and Tasks A/B convention). Pre-register in script header. The computation is deterministic; the seed is documentation-only.

---

## Method

1. **Generate the 5-series primes.** Generate primes p, filter to those with p mod 6 == 5, take the first N for some N covering at least the tight regime + post-regime + asymptotic ranges (N ≥ 200 is comfortable; 5-series primes are dense enough that this is cheap).

2. **Anchor verification.** Compute the table from §"Anchor verification" above for n ∈ [1, 31] and confirm row-by-row match. If any row differs, halt and report. (This is the §"Method step 4" of Task Z, refactored to verify against the explicit table rather than the implicit anchor.)

3. **Compute χ₅ at each 5-series prime** via Legendre symbol mod 5 (or directly via p mod 5).

4. **Compute running sum S₅(n)** for n = 1, 2, …, 200 under Definition A.

5. **Verify Paper 199 §6.1's tie-position listing**: tie positions in n ∈ [1, 33] should be exactly {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}. If the computed set differs, halt and report.

6. **Count zeros in each regime** (as in Task Z):
   - Tight regime: n ∈ [1, 33]
   - Post-regime: n ∈ [34, 66]
   - Asymptotic: n ∈ [67, 200]

7. **Compute `f_tight`, `f_post`, `f_asymp`** as count(zeros) / regime-length. Apply pre-registered decision rule to `f_tight`.

8. **Output the full S₅(n) sequence** for n ∈ [1, 50] in `task_z2_results.md` for human inspection.

9. **Report the bin decision and the empirical f-values** in `task_z2_findings.md`.

---

## Outputs

Files to produce in `papers/200_alpha_kappa_closure/v0_3/`:

1. **`task_z2_brief.md`** — this document, committed as the pre-registration of Task Z'.
2. **`task_z2_within_regime_null.py`** — the implementation. Pre-registered seed = 200 in script header. ASCII-safe prints (the Windows cp1252 lesson). Anchor table cross-check hard-coded so the implementation aborts cleanly if even one row diverges.
3. **`task_z2_results.md`** — human-readable results: anchor table reproduction, full S₅(n) for n ∈ [1, 50], regime counts, fractions, bin decision per pre-registered rule.
4. **`task_z2_results.json`** — machine-readable structured results.
5. **`task_z2_findings.md`** — integrated findings: bin decision, brief structural reading, flags for CinC.

---

## Audit-gap discipline

- Pre-registration commit (this brief + locked script): T_A.
- Results commit (Python output + findings): T_B.
- **T_B − T_A ≥ 30 minutes**, as in Task Z.

The audit gap restarts under Task Z' — the failed-brief Task Z run is closed; Task Z' is a fresh pre-registration with its own audit clock.

---

## What Mr Code may decide (Standing Trust GREEN)

- Seed pre-registration mechanics (already specified as 200; documentation-only as the computation is deterministic).
- Prime-generation library (sympy / numpy / hand-rolled sieve). Use what's already in the environment from Task Z.
- Whether to include a figure (running sum S₅(n) vs n). If included, save PNG and reference from results.md.
- Layout details of results.md and findings.md.

---

## What Mr Code should NOT decide

- The 5-series definition (primes ≡ 5 mod 6, Definition A including ramified p = 5).
- The bin boundaries (LOW <0.30, MID [0.30, 0.55], HIGH >0.55).
- The regime boundaries ([1, 33], [34, 66], [67, 200]).
- The anchor verification table — if any row differs from CinC's table above, halt and flag.

If during implementation any of these need to change, halt and flag to CinC. Do not adjust mid-run.

---

## Cross-link to Task Z (audit trail)

This brief supersedes Task Z (the brief at commit `61aff82` and its locked script with the wrong convention). Task Z's artefacts (`task_brief.md`, `task_z_within_regime_null.py`, `task_z_results.md`, `task_z_results.json`, `findings.md`) are preserved in place as the failed-brief record. The git history shows: Task Z committed `61aff82` (wrong-brief pre-reg) → no results commit (halted on anchor) → Task Z' pre-reg commit (this brief) → Task Z' results commit. The audit trail is therefore: brief-was-wrong, halt-discipline-protected, brief-corrected, re-run, results.

Findings.md for Task Z' should briefly reference Task Z (citing CinC's brief error and Mr Code's halt) for posterity. Paper 200 v0.3's §9 Reproducibility section will likewise cite both Task Z and Task Z' as the audit-trail record.

---

## Expected outcome (CinC's prior, non-binding)

`f_tight = 14/33 = 0.4242…` to within rounding, landing in **Bin MID**. This confirms Paper 199 §6.1's documented figure and Mr A's predicted reading. Paper 200 v0.3 abstract and §1 load-bearing-claim block will be rewritten to acknowledge the regime-density null is consistent with the structural reading at roughly even odds, per Mr A's NULL paragraph.

Should `f_tight` land in LOW or HIGH unexpectedly, halt and re-verify the anchor + tie-position listings (§"Anchor verification" and Method step 5) before applying the bin decision. The hand-verification above is solid; any divergence is a programme-level finding to flag.

---

Standing Trust GREEN throughout. Pre-registration first, ≥30 min audit gap, results within bins, no methodology drift mid-run. Failed-brief Task Z artefacts preserved.

🐕☕⬡

CinC out.
