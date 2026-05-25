# Paper 200 v0.5 Task AA — findings

**Bin MISMATCH-MINOR.** All substantive predictions hold; Paper 199 §3's Def A/B/C robustness claim is verified under direct computation; Mr A's diagnosis that the v0.4 §5.1 `+1` figure is wrong is confirmed; the only mismatch is a labelling swap between CinC's working reading of "Def B" and Paper 199 §3's actual Def B vs Def C.

- Pre-registration commit: `84b14f4`
- Locked script: `task_aa_definitions_verification.py` (unmodified after run)
- Authoritative source for definitions: Paper 199 §3, lines 122–128 and 54 of `Paper_199_RunningSumOf5_Series Primes_v1_2.md`

---

## What Paper 199 §3 actually says

Verbatim (Paper 199 lines 124–126):

- **Def A** (include n = 1): count positions k ∈ {1, ..., 33} with `S_5(k) = 0`. Natural-order count = 14.
- **Def B** (exclude n = 1): count positions k ∈ {2, ..., 33} with `S_5(k) = 0`. Natural-order count = 13.
- **Def C** (omit p = 5 entirely; 32-step sum over non-ramified primes): natural-order count = 13.

And line 128: *"the three differ only in bookkeeping of the ramification step and are not three independent tests."*

The key structural point — and the one that resolves Mr A's confusion in v0.4 §5.1 — is that **Def A and Def B use the same running sum**. They differ only in whether the tie at index `k = 1` (which is `S_5(1) = 0` from `χ_5(5) = 0`) is counted. **Def C uses a different running sum** (over non-ramified primes only), with index `j` corresponding to Def A index `j + 1`.

CinC's working reading had Def B as "shift indices by 1". That description matches Paper 199's Def C, not Def B. The substantive prediction (running sum = 0 at the boundary-crossing prime p = 281 under all three definitions) holds in both readings; only the label is swapped.

---

## Verification cascade — all substantive checks pass

| Check | Expected | Computed | Status |
|---|---|---|---|
| `p_{5,1}` | 5 | 5 | **OK** |
| `p_{5,5}` | 29 | 29 | **OK** |
| `p_{5,31}` | 281 | 281 | **OK** |
| `p_{5,33}` | 311 | 311 | **OK** |
| `S_A(31)` | 0 | 0 | **OK** |
| Def A tie set in [1, 33] | 14 positions | identical | **OK** |
| Def A tie count (Paper 199 §3 window [1, 33]) | 14 | 14 | **OK** |
| Def B tie count (Paper 199 §3 window [2, 33]) | 13 | 13 | **OK** |
| Def C tie count (Paper 199 §3 window [1, 32]) | 13 | 13 | **OK** |
| `S_5 = 0` at p = 281 under all three defs | True | True | **OK** |

Paper 199 §3's robustness claim is verified to the digit.

---

## Running sum at the boundary-crossing prime — the verification table for v0.5 §5.1

| Definition | Index where p = 281 sits | `S_5` at that index | Tie? |
|---|---:|---:|:---:|
| **Def A** | k = 31 | **0** | ✓ |
| **Def B** | k = 31 (same indexing as A) | **0** | ✓ |
| **Def C** | j = 30 (after omitting ramified p = 5) | **0** | ✓ |

The same underlying observation — "the running sum is zero at the boundary-crossing prime p = 281" — appears under all three definitions, at the index appropriate to each definition's indexing convention.

---

## Tie-position structure across definitions

- **Def A** [1, 33], 14 ties: `{1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}`
- **Def B** [2, 33], 13 ties: `{3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}` (= Def A tie set minus `{1}`)
- **Def C** [1, 32], 13 ties: `{2, 4, 6, 8, 10, 14, 16, 18, 20, 22, 26, 30, 32}` (= Def A tie set in [2, 33], each index minus 1)

The Def B and Def C tie sets are *the same set of underlying primes*, indexed differently. The Def A tie set is the same plus the ramified-step tie at `k = 1`.

---

## CinC literal-prediction scoring

| Prediction | CinC value | Computed | Match |
|---|---:|---:|:---:|
| `S_A(31)` | 0 | 0 | ✓ |
| `S_B(30)` under Paper 199 §3 Def B | 0 | **−1** | ✗ |
| `S_C` at boundary p = 281 (j = 30 in Def C) | 0 | 0 | ✓ |
| Def A 14-tie set in [1, 33] | 14 | 14 | ✓ |
| Def B 13 ties under "shift by 1" reading | 13 | 13 | ✓ (rationale differs, count happens to coincide) |

The single literal failure (`S_B(30)`) is the labelling slip: CinC's prediction was the value under Paper 199 §3's Def C at index 30 (which is 0), labelled as the value under Def B at index 30. Under Paper 199 §3's actual Def B, `S_B(30) = S_A(30) = −1` (per the Task Z' table).

The Def B 13-tie count happens to match because, by coincidence, the count of ties under Paper 199 §3's actual Def B (counts in window [2, 33]) equals the count CinC predicted under the "shift by 1" reading (which would also give 13 ties under that re-indexed sequence). The number is right; the route to it differs.

---

## Mr A's round-3 item 2 — verdict

Mr A's diagnosis is **confirmed**.

> "the value should be 0, not +1"

Correct under both Paper 199's actual Def B (which gives `S_B(30) = −1`) and CinC's intended Def C (which gives `S_C(30) = 0`). Either way, the v0.4 figure of `+1` is wrong. The substantive claim it was supporting — "Definition A is load-bearing because Def B would give a different value" — fails for a deeper reason than the arithmetic: under Paper 199 §3, Def A and Def B *use the same running sum*, so claiming Def A is substantively load-bearing relative to Def B was never going to survive. The honest reframing is that Def A is the *natural* convention for the reset reading (since the j=1 zero baseline from the ramified prime makes the "reset" framing explicit at the start), with the substantive observation robust across A, B, C per Paper 199 §3.

---

## Proposed v0.5 §5.1 parenthetical text (Mr Code's draft; CinC to adopt or revise)

To replace the v0.4 parenthetical that triggered round-3 item 2, the following short table is directly insertable. It (a) reflects Paper 199 §3 verbatim, (b) carries the verified table from Task AA, and (c) re-frames Def A's role as "natural convention" rather than "substantively load-bearing."

> The convention chosen here is Definition A of Paper 199 §3, which includes the ramified prime `p = 5` (contributing `χ_5(5) = 0`) as the first step of the running sum. Definition A makes the 5-series's natural zero baseline at `n = 1` explicit, which is the reset reading the Crossing 1 anchor `S_5(31) = 0` is most naturally stated against. The substantive observation — that the running sum is zero at the boundary-crossing prime `p_{5,31} = 281` — is robust across all three Paper 199 §3 definitions, which differ only in the bookkeeping of the ramification step:
>
> | Definition | Same prime p = 281 at index | Running sum at that index |
> |---|---:|---:|
> | Def A (include `n = 1`) | `k = 31` | **0** |
> | Def B (include `n = 1` in the sum but not in the tie window) | `k = 31` | **0** |
> | Def C (omit `p = 5` entirely; sum over non-ramified primes) | `j = 30` | **0** |
>
> Per Paper 199 §3 the three definitions "differ only in bookkeeping of the ramification step and are not three independent tests"; Definition A is chosen here as the natural convention for the reset reading, not as a substantive load-bearing choice.

The earlier v0.4 framing — claim that Def A was load-bearing because Def B would give a different value — is retracted. The arithmetic that motivated it (the `+1` figure) was a Def B/Def C labelling slip on CinC's part, caught by Mr A.

---

## Audit-gap record

- Task AA pre-registration commit `84b14f4`: T_A = 2026-05-16T14:40:27+01:00.
- Task AA results commit: T_B at commit time, ≥30 min after T_A.

---

## Flags

1. **PRIMARY (decision):** **Bin MISMATCH-MINOR**. Substantive predictions all hold; Paper 199 §3 robustness verified; Mr A's diagnosis confirmed. The labelling correction (CinC's "Def B" reading is Paper 199's Def C) is documented; no impact on the substantive claim direction. v0.5 §5.1 should adopt the verified table above and reframe Def A as "natural convention" rather than "substantively load-bearing."

2. **PROCESS (recognition-log candidate):** Task AA's pre-registration scored CinC's predictions both substantively *and* literally. The substantive predictions passed cleanly; the literal mismatch caught the Def B/Def C labelling slip in CinC's working reading without needing the kind of HALT-and-rebrief cycle the Task Z → Task Z' sequence required. Same brief, two scoring layers — and the pre-registration discipline still caught the slip. The mechanism worked, just at a finer grain.

3. **MINOR (worth noting for v0.5 §5.1 prose):** Paper 199 §3 itself uses Def C as primary ("Definition C is the principal one as it measures what the paper claims to measure"). Paper 200 v0.5 §5.1 — which states the Crossing 1 anchor — naturally uses Def A (as the reset reading). The asymmetry is consistent: Paper 199 measures the split/inert balance directly (Def C), Paper 200 anchors the reset reading at the natural zero baseline (Def A). No conflict; just a difference in what each paper foregrounds.

🐕☕⬡

Mr Code, in band.
