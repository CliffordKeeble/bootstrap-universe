# Task Z — Within-Regime Null Test for Paper 200 v0.3

**Brief for Mr Code · CinC handoff · 12 May 2026**
**Standing Trust: GREEN** (per CLAUDE.md and prior Mr Code work this session)

---

## Background

Paper 200 v0.2.2 (in `/mnt/user-data/outputs/paper_200_v0_2_2.md`) was reviewed by Mr Adversary v2.3 round-1 in fresh context. Mr A returned three stars with "I would publish this … send it back with the null test attached and we shall see where it stands." The principal concern is in his NULL paragraph:

> Paper 199 §6.1 records a ~42% tie density in the 5-series for n ≤ 33. The probability under the null hypothesis that S₅(31) = 0 for any reason whatever is therefore on the order of 0.42. That is not a discriminating observation; it is consistent with the structural reading and consistent with chance at roughly even odds. […] Before publication, run the simple within-regime null I have already half-described: across all 30-step trailing windows ending in the tight-oscillation regime (n ≤ 33), what fraction terminate at S₅ = 0? If the answer is anywhere near 0.42, the load-bearing claim of Paper 200 needs to either acquire a second independent anchor or be scope-limited correspondingly.

The Paper 200 load-bearing claim — first staggered gap at 10⁵⁹–10⁶¹ in mod-6 prime products = α → κ closure event, with χ₅ reset at S₅(31) = 0 as conjunctive empirical signature — is currently anchored on Crossing 1 (the S₅(31) = 0 observation at the boundary-crossing prime n = 31). Crossings 2-4 are predicted-non-zero observations which the 7-series satisfies by default (Paper 199 §4: two ties in five thousand primes). The discriminating power of the four-crossing design is borne almost entirely by Crossing 1.

If the null fraction is ~0.42 as Paper 199 §6.1 documents, Crossing 1 alone is barely better than coin-flip discrimination against the regime-density alternative. Paper 200 v0.3 must reflect this in the abstract and §1 load-bearing-claim block, regardless of which way the test lands. The test fixes the framing; the framing direction is what's at stake.

This task runs the explicit computation Mr A asked for, with pre-registered decision criteria. Result feeds Paper 200 v0.3.

---

## Question

**What fraction of values of S₅(n) for n in the tight-oscillation regime (n ≤ 33) are exactly zero?**

Under Definition A (Paper 199's convention: include the ramified prime p = 5 with χ₅(5) = 0; running sum indexed by prime ordinal n where p_n is the n-th prime, starting p_1 = 2).

Auxiliary: compute the same fraction for two further regimes for comparison:
- Post-regime n ∈ [34, 66] (still small-n but past the tight oscillations Paper 199 §6.1 documents)
- Asymptotic n ∈ [67, 200] (range where Paper 199 §6.1 reports tie density drops substantially)

---

## Pre-registered decision criteria

**The decision rule is binding once committed to git. Do not move the bin boundaries after results.**

Let `f_tight` = fraction of n ∈ [1, 33] with S₅(n) = 0.

- **Bin LOW**: `f_tight < 0.30`. The structural reading is materially stronger than the regime-density null. Paper 200 v0.3 may keep the v0.2.2 abstract/§1 framing largely intact, with a tightening sentence acknowledging that the test was run and the structural reading survives.
- **Bin MID**: `f_tight ∈ [0.30, 0.55]`. The regime-density null is consistent with the observation at roughly even odds. Paper 200 v0.3 abstract and §1 load-bearing-claim block must be rewritten to acknowledge that Crossing 1 is barely discriminating against the regime-density alternative, and Paper 200 either finds a second independent anchor (out of scope for v0.3) or scope-limits the claim correspondingly.
- **Bin HIGH**: `f_tight > 0.55`. The regime-density null is stronger than the observation. Paper 200 v0.3 must substantially scope-limit, reframing the structural reading as "consistent with the regime statistics but not selected by them" and reducing the strength of the load-bearing-claim block accordingly.

Mr A's reading (per Paper 199 §6.1 figure of ~14/33 ≈ 0.424) predicts Bin MID. The test confirms or refutes this prediction.

**Random seed: 200** (paper number convention from prior tasks). Pre-register in script header.

---

## Method

1. **Generate primes.** Sieve of Eratosthenes (or sympy.primerange) for the first 200 primes. Verify p_1 = 2, p_5 = 11, p_31 = 127 (cross-check against any standard prime table).

2. **Compute χ₅(p) for each prime.** The Dirichlet character mod 5:
   - χ₅(p) = +1 if p ≡ 1, 4 (mod 5)
   - χ₅(p) = −1 if p ≡ 2, 3 (mod 5)
   - χ₅(p) = 0 if p ≡ 0 (mod 5), i.e. p = 5 (the ramified prime)

3. **Compute running sum S₅(n) = Σ_{k=1}^{n} χ₅(p_k)** for n = 1, 2, …, 200 under Definition A.

4. **Verify the structural anchor: S₅(31) = 0.** This is the observed datum from Paper 199 / Paper 200 §5. If it doesn't reproduce, halt and report — something is wrong with the convention or the primes.

5. **Count zeros in each regime:**
   - Tight regime: n ∈ [1, 33]
   - Post-regime: n ∈ [34, 66]
   - Asymptotic: n ∈ [67, 200]

6. **Compute `f_tight`, `f_post`, `f_asymp`** as count(zeros) / regime-length.

7. **Apply pre-registered decision rule** to `f_tight` and report the bin.

8. **Output the full S₅(n) values** for n ∈ [1, 50] in `results.md` for human inspection — Mr A may want to see the actual sequence.

9. **Compute Paper 199 §6.1 cross-check**: extract the count of n ∈ [1, 33] with S₅(n) = 0 and verify it matches the "14 ties in 33 primes" figure Paper 199 §6.1 reported (if my reading of Mr A's interpretation is correct). If the actual count differs from 14, flag this — Paper 199's framing may need its own audit.

---

## Outputs

Files to produce in `papers/200_alpha_kappa_closure/v0_3/`:

1. **`task_brief.md`** — this document, committed as the pre-registration.
2. **`task_z_within_regime_null.py`** — the implementation. Pre-registered seed = 200 in script header. ASCII-safe prints (lesson from Task A/B Windows encoding gremlin).
3. **`task_z_results.md`** — human-readable results: full S₅(n) table for n ∈ [1, 50], counts per regime, fractions, bin decision per pre-registered rule.
4. **`task_z_results.json`** — machine-readable structured results.
5. **`findings.md`** — integrated findings: pre-registered bin decision, brief structural reading of the result, flags for CinC.

---

## Audit-gap discipline

- Pre-registration commit (this brief + bin definitions in script header): commit hash `<HASH_A>`, timestamp T_A.
- Results commit (Python output + findings.md): commit hash `<HASH_B>`, timestamp T_B.
- **Audit gap must satisfy T_B − T_A ≥ 30 minutes.** Per programme convention. Use this time to sanity-check the JSON, write the findings narrative, and verify the implementation against the brief.

If the audit gap is short and you're confident in the result, hold the results commit until the timer fires. Don't shortcut the discipline.

---

## What Mr Code may decide (within Standing Trust GREEN)

- Seed pre-registered in script header (already specified as 200, but if a different seed is more natural for this computation, flag and choose it).
- Choice of prime-generation library (sympy / numpy / hand-rolled sieve). Use what's already in the environment.
- Whether to include a figure (running sum as a function of n). If included, save as PNG in the v0_3 folder and reference from results.md.
- Layout of results.md and findings.md. Use Task A/B as the template; deviate where it makes the report cleaner.

---

## What Mr Code should NOT decide

- The bin boundaries (LOW <0.30, MID [0.30, 0.55], HIGH >0.55) are pre-registered and binding.
- Definition A is binding (include ramified p = 5, contributes 0 to the sum). Definition B would shift indices and is out of scope for this task.
- The regime boundaries (tight n ∈ [1, 33], post n ∈ [34, 66], asymptotic n ∈ [67, 200]) are pre-registered.

If during implementation any of these need to change, halt and flag to CinC. Do not adjust mid-run.

---

## Expected outcome (CinC's prior, for Mr Code's reference only — not binding on decision)

Mr A's quoted figure of ~0.42 (14 ties / 33 primes from Paper 199 §6.1) suggests `f_tight` will land in Bin MID. This is my prior but the test is what fixes the framing — if I'm wrong, the framing moves accordingly. Run the test, apply the rule, report the bin.

A LOW result would substantially strengthen Paper 200 v0.3.
A MID result confirms Mr A's reading and forces v0.3 framing tightening per the decision rule.
A HIGH result would force substantial scope-limiting on Paper 200 v0.3.

All three are publishable outcomes for Paper 200; only the framing differs.

---

## Flag for Mr Code

If the implementation reveals that Paper 199 §6.1's "14 ties" figure does not match the count of S₅(n) = 0 events in [1, 33] (i.e., "ties" in Paper 199 means something other than what I've assumed here — e.g., consecutive equal values, or some other tie-definition), flag this to CinC. Paper 199 §6.1 may need its own audit if the definitions don't line up. This is not blocking for Task Z's decision rule (which uses the directly-computed `f_tight` regardless of what Paper 199 §6.1 meant by "ties"), but it's a programme-level finding worth recording in findings.md.

---

Standing Trust GREEN throughout. Pre-registration first, ≥30 min audit gap, results within bins, no methodology drift mid-run.

🐕☕⬡

CinC out.
