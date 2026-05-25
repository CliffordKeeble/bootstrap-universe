# Task AA — Definition A/B/C Verification for Paper 200 v0.5 §5.1

**Brief for Mr Code · CinC handoff · 16 May 2026**
**Standing Trust: GREEN**

---

## Background

Paper 200 v0.4 was reviewed by Mr Adversary v2.3 round-3 (four stars, fifth-star path itemised across seven points). The most substantive of the seven is round-3 item 2 (STATUS/NULL — §5.1 Definition B parenthetical), where Mr A observes:

> The convention note on line 222 ("S_5(30) under Definition B = S_5(31) under Definition A on the next non-ramified prime, which is +1, after p_{5,32} = 293 contributes χ_5(293) = −1") was the sentence I had to re-read three times. On my reading, if Def B is "exclude the ramified prime from indexing", then S_5(30) under Def B and S_5(31) under Def A should be arithmetically equal (since χ_5(5) = 0 contributes nothing in either), and the value should be 0, not +1. … I would ask Mr Code to verify the Def-A vs Def-B values at the relevant indices and rewrite the parenthetical as a small table or as a fully-worked one-liner.

On CinC re-check, Mr A's diagnosis appears correct: the v0.4 parenthetical value (+1) is wrong, and the substantive claim it supports ("Definition A is load-bearing because Def B would give a different value") does not survive. Paper 199 §3 explicitly states that Definitions A, B, C agree on the substantive observation and differ only in bookkeeping of the ramification step. The v0.4 §5.1 framing overstates.

Task AA verifies the Def-A / Def-B / Def-C values at the relevant indices under pre-registration discipline, in order to (i) confirm or refute CinC's correction of the v0.4 figure, (ii) provide a small verified table suitable for direct inclusion in v0.5 §5.1, and (iii) close round-3 item 2 with the audit-trail discipline that all Paper 200 computational claims now carry.

This is a verification task, not a new null test. The pre-registration is of CinC's predicted values; Mr Code confirms or denies. If denied, v0.5 §5.1's framing needs to find its support elsewhere — Mr A's exact phrasing on the consequence.

---

## Convention specifications (Paper 199 §3, formal)

Mr Code should read Paper 199 §3 at `/mnt/user-data/uploads/Paper_199_RunningSumOf5_Series_Primes_v1_2.md` for the formal Definition A / B / C statements before computing. CinC's working understanding from Paper 199 §3 (subject to Mr Code's authoritative reading):

- **Definition A**: running sum S_5(n) = Σ_{j=1}^{n} χ_5(p_{5,j}) over the first n primes ≡ 5 (mod 6), starting with p_{5,1} = 5 (ramified, contributes χ_5(5) = 0). Indexing includes the ramification step.
- **Definition B**: as Def A but with the ramification step excluded from indexing — index by non-ramified primes only. The j-th non-ramified prime is p_{5,j+1} under Def A's indexing. S_5(j) under Def B = Σ_{k=1}^{j} χ_5(p_{B,k}) where p_{B,k} = p_{5,k+1} under Def A.
- **Definition C**: Paper 199 §3's principal definition. CinC's working reading: exclude the ramification step from both indexing and the window — sum over the first N non-ramified primes for some N consistent with Paper 199's "33-prime window" framing. Mr Code should confirm Def C's exact statement against Paper 199 §3.

If Mr Code's reading of Paper 199 §3 differs from the above, the brief's pre-registered predictions need to be re-anchored to the authoritative definitions. Halt and flag if so.

---

## Question

For each definition D ∈ {A, B, C}, compute S_5(n) under D for n covering the tight-oscillation regime [1, 33] and the boundary-crossing prime p = 281. Specifically:

1. Confirm S_5(31) = 0 under Def A (the anchor, already verified in Task Z').
2. Compute S_5(j) under Def B at the index j corresponding to the boundary-crossing prime p_{5,31} = 281 — i.e., S_5(30) under Def B if indexing shifts by 1.
3. Compute S_5(k) under Def C at the corresponding index k (per Paper 199 §3's Def C convention).
4. Verify Paper 199 §3's robustness claim: that A, B, C agree on the substantive observation, with tie positions corresponding up to indexing shifts.
5. Produce a small verified table suitable for direct inclusion in Paper 200 v0.5 §5.1.

---

## Pre-registered predictions

**Random seed: 200** (paper number; documentation-only as the computation is deterministic).

CinC's predictions, binding for the decision rule:

| Observation | Predicted value | Rationale |
|---|---|---|
| S_5(31) under Def A at p_{5,31} = 281 | **0** | Anchor, verified in Task Z' (c2c5eeb) |
| S_5(30) under Def B at p_{B,30} = 281 | **0** | Indices shift by 1 (skip ramified j=1); S_5(j) under Def B = S_5(j+1) under Def A; at j=30, this equals S_5(31) under Def A − χ_5(5) = 0 − 0 = 0 |
| S_5(?) under Def C at index for p = 281 | **0** | Per Paper 199 §3's robustness claim; index depends on Def C convention |
| Tie positions {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33} under Def A | 14 ties | Reproduced from Task Z' (c2c5eeb) |
| Tie positions under Def B | 13 ties, shifted by 1 | Same substantive observations, ramification step (j=1 under Def A) excluded from Def B indexing |
| Tie positions under Def C | per Paper 199 §3 | Mr Code's authoritative reading |

---

## Pre-registered decision rule

**The decision rule is binding once committed to git.**

- **Bin CONFIRM**: all three definitions give the same substantive observation at correspondingly-shifted indices, matching CinC's predictions to the digit. Paper 200 v0.5 §5.1 rewrites the parenthetical with the verified table; the claim "Definition A is load-bearing because Def B would give a different value" is retracted; the honest reframing is "Definition A is the natural convention for the reset reading because it makes the 5-series's natural zero baseline at j=1 (from p=5's ramification) explicit, with the substantive observation robust across A, B, C per Paper 199 §3."

- **Bin MISMATCH-MINOR**: predictions hold but Paper 199 §3's Def C convention turns out to be different from CinC's working reading. v0.5 §5.1 rewrites the parenthetical with the correct Def C accounting; no impact on the substantive claim direction (Def A's load-bearing status remains "natural convention for reset reading" not "substantive dependence").

- **Bin MISMATCH-MAJOR**: any of the predicted values is wrong, or Paper 199 §3's robustness claim does not hold under direct computation. Halt and report — this would indicate either an error in Paper 199 or a deeper misunderstanding of the conventions. Paper 200 v0.5 cannot be drafted until this is resolved.

CinC's prior: Bin CONFIRM, very likely. The Bin MISMATCH-MAJOR scenario would be a programme-level finding.

---

## Method

1. **Read Paper 199 §3.** Confirm authoritative statements of Definitions A, B, C. If divergent from CinC's working understanding above, flag and proceed under Paper 199's authoritative definitions.

2. **Generate 5-series primes.** First 50 primes ≡ 5 (mod 6); cross-check that p_{5,1} = 5, p_{5,5} = 29, p_{5,31} = 281, p_{5,33} = 311 (CinC anchor table from Task Z' / Paper 119 v2.0 verification).

3. **Compute χ_5(p) for each prime.** Legendre symbol mod 5; χ_5(5) = 0.

4. **Compute S_5(n) under each definition** for n covering [1, 50]:
   - Def A: include the ramification step
   - Def B: exclude the ramification step from indexing (shift by 1)
   - Def C: per Paper 199 §3's authoritative statement

5. **Cross-check Task Z'**: S_5(31) under Def A = 0; 14 ties at {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33} under Def A. If divergent, halt — something is wrong with the script.

6. **Produce verification table** showing for each definition: (i) the index at which the substantive observation "running sum = 0 at the boundary-crossing prime p = 281" sits, (ii) the value of S_5 at that index, (iii) cross-confirmation that the same primes are summed (just at different indexing). The table should be directly insertable into Paper 200 v0.5 §5.1.

7. **Apply pre-registered decision rule** and report the bin.

---

## Outputs

Files to produce in `papers/200_alpha_kappa_closure/v0_5/` (new directory for v0.5 work):

1. **`task_aa_brief.md`** — this document, pre-registered.
2. **`task_aa_definitions_verification.py`** — implementation. Pre-registered seed = 200 in script header. ASCII-safe prints (Windows cp1252 lesson). Anchor cross-check against Task Z' values hard-coded.
3. **`task_aa_results.md`** — verification table, regime tie listings under each definition, full S_5(n) sequences for n ∈ [1, 50] under each definition.
4. **`task_aa_results.json`** — machine-readable structured results.
5. **`task_aa_findings.md`** — integrated findings: bin decision, brief structural reading, recommended v0.5 §5.1 parenthetical text (Mr Code may propose; CinC will adopt or revise).

---

## Audit-gap discipline

- Pre-registration commit (this brief + locked script): T_A.
- Results commit (Python output + findings): T_B.
- **T_B − T_A ≥ 30 minutes**, standard programme rule.

If a wrong-brief halt fires mid-task (per the protocol revision in `Cliff_Patterns_Working` after the Task Z / Task Z' sequence), let the original audit clock fire before committing the halt record. Discipline is symmetric.

---

## What Mr Code may decide

- Choice of prime-generation library (sympy / numpy / hand-rolled). Use what's already in the environment from Task Z'.
- Whether to produce additional auxiliary tables (e.g., side-by-side S_5(n) under all three definitions). Helpful for v0.5 §5.1 if compact.
- Layout of results.md and findings.md.

## What Mr Code should NOT decide

- The pre-registered predictions (above table) are binding once committed.
- The bin boundaries are pre-registered.
- Paper 199 §3's Def A / Def B / Def C statements are authoritative — if Mr Code's reading of Paper 199 §3 differs from CinC's working reading, halt and flag rather than proceed under a different convention.

If anything in the implementation needs to drift, halt and flag.

---

## Expected outcome (CinC's prior, non-binding)

**Bin CONFIRM.** Paper 199 §3's robustness claim is verified under direct computation; Mr A's diagnosis of the v0.4 +1 figure as wrong is confirmed; v0.5 §5.1 rewrites with the verified table; the load-bearing-claim direction shifts from "Def A is substantively load-bearing" to "Def A is the natural convention for the reset reading; the observation is robust across definitions."

Bin MISMATCH-MAJOR would be programme-level news worth flagging immediately (and unblocking Paper 200 v0.5 from drafting until resolved).

---

Standing Trust GREEN throughout. Pre-registration first, ≥30 min audit gap, results within bins, no methodology drift mid-run. Discipline symmetry per the Cliff_Patterns_Working entry from the Task Z / Task Z' sequence.

🐕☕⬡

CinC out.
