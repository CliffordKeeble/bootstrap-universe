# Paper 200 v0.3 Task Z — findings

**Outcome: HALT on structural anchor.** No bin decision was reached.
The pre-registered locked script halted at brief method step 4 (anchor
verification) and did not proceed to compute regime fractions.

- Pre-registration commit: `61aff82`
- Locked script: `task_z_within_regime_null.py` — unmodified after the
  run. The script is the artefact of record.
- Brief's halt instruction (method step 4): *"If it doesn't reproduce,
  halt and report — something is wrong with the convention or the
  primes."*
- Brief's adjust-mid-run instruction (NOT decide section): *"If during
  implementation any of these need to change, halt and flag to CinC.
  Do not adjust mid-run."*

Both instructions point the same way. The script halted; this report is
the flag.

---

## What was checked, and what failed

**Pre-flight checks (passed):**

| Check | Expected | Computed | Status |
|---|---|---|---|
| `p_1` | 2 | 2 | OK |
| `p_5` | 11 | 11 | OK |
| `p_31` | 127 | 127 | OK |

The prime generator is correct.

**Anchor check (failed):**

| Quantity | Expected | Computed (Def A) | Computed (Def B diag.) |
|---|---|---|---|
| `S_5(31)` | 0 | **−6** | **−5** |

Definition A is the brief's binding convention (ramified prime `p = 5`
included with `χ_5(5) = 0`, 1-indexed by prime ordinal). Definition B
is the diagnostic alternative (skip the ramified prime; re-index). The
anchor reproduces under neither.

**Cross-check (failed):**

| Quantity | Expected (per brief) | Computed (Def A) |
|---|---|---|
| Count of `n ∈ [1, 33]` with `S_5(n) = 0` | 14 | **0** |

Not 14 ties, not 1 tie — zero ties. `S_5(n)` is strictly negative
across the entire tight-oscillation regime under Definition A
(values in {−6, −5, −4, −3, −2, −1}; never 0).

This is the secondary flag the brief explicitly anticipated:

> If the implementation reveals that Paper 199 §6.1's "14 ties" figure
> does not match the count of S₅(n) = 0 events in [1, 33] (i.e., "ties"
> in Paper 199 means something other than what I've assumed here ...),
> flag this to CinC.

But the secondary flag is now subsumed by the primary one: the anchor
itself fails, not just the §6.1 cross-check.

---

## Convention audit — character is consistent with the repo

To rule out a chi_5 implementation error, I cross-checked against
`golden-chebyshev/character_projections.py`, which encodes the mod-5
character table for prior Paper 196 work:

```
chi_1:  r=1: +1,  r=2: −1,  r=3: −1,  r=4: +1     (Legendre (r/5))
```

This is identical to my Definition A values. The character is the
standard real Dirichlet character mod 5, the Legendre symbol (r/5),
with the ramified prime contributing 0. Implementation is consistent
with the rest of the programme.

---

## What the data actually shows for `n ∈ [1, 33]`

Under Definition A, `S_5(n)` over the tight-oscillation regime:

```
n:     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
S_5:  -1  -2  -2  -3  -2  -3  -4  -3  -4  -3  -2  -3  -2  -3  -4  -5  -4

n:    18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33
S_5:  -3  -4  -3  -4  -3  -4  -3  -4  -3  -4  -5  -4  -5  -6  -5  -6
```

Range: `[−6, −1]`. **`S_5` never crosses zero in `[1, 33]`.** The
sequence is mildly oscillating but with a strong negative trend driven
by an early imbalance: primes ≡ 2 or 3 (mod 5) outpace primes ≡ 1 or 4
(mod 5) in this regime. By `n = 31` there are 18 negative contributions
vs 12 positive, plus one zero at `p = 5`, giving `S_5(31) = 12 − 18 = −6`.

The "tight oscillation" framing of Paper 199 §6.1 may describe a real
phenomenon — `S_5` does oscillate visibly between −2 and −6 here — but
the oscillation is around a negative offset, not around zero. The
"~42% tie density" figure Mr A quoted from Paper 199 §6.1 does not
correspond to `S_5(n) = 0` events under Definition A.

---

## Candidate causes (for CinC to consider)

Listed without preference. I am NOT picking one and running with it
— that's exactly the "adjust mid-run" move the brief forbids. The
purpose here is to give CinC a starting set.

1. **`S_5(31) = 0` may be a Paper 200 §5 framing slip.** If Paper 200
   v0.2.2 §5 reports `S_5(31) = 0` as the observed datum, the datum
   itself may not be reproducible under the convention the same paper
   defines. Worth re-reading Paper 200 §5 to check whether the value
   is asserted or derived, and which prior result it cites.

2. **"Ties" in Paper 199 §6.1 may mean something other than `S_5 = 0`.**
   Candidates: (a) `S_5(n) = S_5(n − k)` for some lag k, i.e. revisits
   to a previously-attained sum value; (b) sign changes of `S_5(n)`
   (zero crossings, which still wouldn't happen here under Def A since
   `S_5` stays negative); (c) local maxima/minima of `S_5(n)`;
   (d) a different statistic entirely. Under Def A, `S_5` visits −3
   eight times in `[1, 33]`, −4 nine times, −2 four times — there are
   many "revisits" but they're not zeros. If Paper 199 §6.1 actually
   reports revisits-to-some-recurring-value, the "14" figure becomes
   plausible (e.g., −4 alone is hit 9 times, and combining −3 and −4
   gives 17). This is speculation; the actual count depends on the
   exact definition.

3. **The "5-series" in Paper 199 may not be `S_5` as defined here.**
   Could be a different running sum: e.g., a *signed* count weighted
   differently, or `T_5(n) = Σ χ_5(p_k) · log(p_k)` (the Chebyshev-style
   sum used in `golden-chebyshev/`), or one of the order-4 mod-5
   characters (chi_2 or chi_3 in `character_projections.py`). Worth
   checking which sum Paper 199 §6.1 actually plots.

4. **Indexing may differ.** Brief specifies "prime ordinal n, p_1 = 2".
   Alternatives: indexing by integer n (so `S_5(31)` sums χ_5 over
   *integers* 1..31, not primes), or by odd-prime ordinal (so `n = 1`
   means `p = 3`, shifting everything by one). Under each, the anchor
   index 31 picks up a different `S_5` value.

5. **Convention may use a different ramification rule.** Some authors
   define `χ_5` to map all multiples of 5 to 0 in the running sum over
   integers; if we're summing over primes specifically, the only
   multiple of 5 that's prime is 5 itself, which is what Def A already
   handles. But if the running sum were over integers, the anchor
   computation changes completely. Worth confirming "primes" is the
   intended domain.

My instinct from the data is that cause (2) or (3) is most likely:
Paper 199's "14 ties" framing isn't naming `S_5 = 0` events, and Paper
200 §5's `S_5(31) = 0` assertion may have inherited the framing without
the data check.

---

## Implications for Paper 200 v0.3

This is bigger than Task Z. The brief framed Task Z as a within-regime
null on the load-bearing claim. The HALT outcome says the load-bearing
*observation* (`S_5(31) = 0`) itself doesn't reproduce under the
convention the brief defines. That's not a null test verdict — that's
a question about whether the Crossing 1 anchor exists under the claimed
convention at all.

**Three orthogonal possibilities, all worth Mr Adversary v2.3 attention
in round-2:**

- The convention named in Paper 200 §5 / Paper 199 §6.1 is correct, and
  the brief's Definition A is a slightly off restatement. → CinC
  re-states Definition A, Task Z re-runs cleanly.
- The convention is correct as stated, and Paper 200 §5's `S_5(31) = 0`
  assertion is wrong. → Paper 200 v0.3 needs the §5 datum audited
  before any null-test work attaches to it.
- The convention is correct as stated, the datum is correct, and the
  brief's identification of `n = 31` with `p_31 = 127` is the slip
  (e.g., the intended anchor is at a different `n`). → Re-identify
  the anchor index, re-run Task Z against it.

I cannot tell which from inside the implementer role.

---

## What I did not do

- Did not modify the locked script (`task_z_within_regime_null.py`).
- Did not adjust the convention to make the anchor pass.
- Did not pick one of the candidate causes and produce a "best guess"
  bin decision under it.
- Did not compute regime fractions (`f_tight`, `f_post`, `f_asymp`)
  under Definition A as a fallback, because the brief's decision rule
  is contingent on anchor verification passing.
- Did not skip the audit-gap discipline — the results commit waits until
  ≥30 min after the pre-registration commit (T_A = 2026-05-15T18:41:30+01:00).

---

## Audit-gap record

- T_A (pre-registration commit `61aff82`): 2026-05-15T18:41:30+01:00
- T_B (results commit): to be recorded at commit time; ≥30 min after T_A.

---

## Flags

1. **PRIMARY:** Structural anchor `S_5(31) = 0` does not reproduce under
   Definition A as stated in the brief. Task Z cannot proceed to a bin
   decision until the anchor convention is clarified.

2. **SECONDARY (subsumed by primary):** Paper 199 §6.1's "14 ties /
   33 primes" figure does not match the count of `S_5(n) = 0` events
   in `[1, 33]` under Definition A (computed count: 0). Paper 199 §6.1
   may need a definition audit independent of Task Z.

3. **PROGRAMME-LEVEL:** Paper 200 §5's `S_5(31) = 0` assertion is the
   anchor of Crossing 1 — currently the load-bearing crossing of the
   four-crossing design. If §5 cites a Paper 199 result that doesn't
   compute the way the brief reads it, this is a chain-of-citation
   audit issue that goes back at least to Paper 199 §6.1.

🐕☕⬡

Mr Code, halted at anchor.
