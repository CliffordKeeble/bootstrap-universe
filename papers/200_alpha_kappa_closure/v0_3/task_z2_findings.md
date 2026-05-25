# Paper 200 v0.3 Task Z' — findings

**Bin decision: MID.** `f_tight = 14/33 = 0.4242…` — CinC's prior confirmed to within rounding. Paper 199 §6.1's "14 ties / 33 primes" figure reproduces exactly; Mr A's NULL reading is the right reading; Paper 200 v0.3's abstract and §1 load-bearing-claim block must be rewritten per the brief's Bin MID branch.

- Task Z' pre-registration commit: `4ae5b36`
- Locked script: `task_z2_within_regime_null.py` — unmodified after run.
- Failed-brief Task Z record (cross-link): commit `61aff82` (wrong-convention pre-reg), commit `85d4ab8` (halt artefacts).

---

## Verification cascade — all checks passed

| Check | Expected | Computed | Status |
|---|---|---|---|
| Anchor table rows 1–31 | CinC hand-table verbatim | 31/31 rows match | **OK** |
| `S_5(31)` | 0 | 0 | **OK** |
| Tie-position set in [1, 33] | {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33} | identical | **OK** |
| Tie count in [1, 33] | 14 (Paper 199 §6.1) | 14 | **OK** |

The implementation reproduces CinC's hand-computation row for row. The convention is settled.

---

## Regime fractions

| Regime | n range | zeros | length | fraction |
|---|---|---:|---:|---:|
| **tight** | [1, 33] | 14 | 33 | **0.4242** |
| post-regime | [34, 66] | 5 | 33 | 0.1515 |
| asymptotic | [67, 200] | 25 | 134 | 0.1866 |

`f_tight = 0.4242` falls cleanly inside the pre-registered MID band `[0.30, 0.55]`.

**Decision against pre-registered rule: Bin MID.**

---

## Structural reading

The three regime fractions tell a fairly clean story:

- The tight regime is roughly **2.5× to 2.8× denser in `S_5 = 0` events** than the post-regime and asymptotic windows. The "tight oscillation" framing in Paper 199 §6.1 is empirically real — there is something distinctive about the [1, 33] window.
- But the *strength* of that distinction is exactly what Mr A's NULL paragraph identifies as the load-bearing weakness: within the tight regime, the unconditional probability of `S_5(n) = 0` is already ~0.42. The Crossing 1 datum `S_5(31) = 0` is one observation drawn from a regime where ~42% of indices have that property. As a discriminator against the regime-density alternative, it is roughly coin-flip.
- A small but worth-noting auxiliary: `f_post = 0.1515` is slightly *below* `f_asymp = 0.1866`. The post-regime is not just "a continuation of the tight regime with lower density" — it is, by this sample, marginally below even the asymptotic rate. This is small-N (n=33), so I would not lean on it; but it weakens any framing that posits a smooth monotonic decay from tight → post → asymptotic. The tight regime appears closer to a regime-specific phenomenon than to the head of a slow decay.

None of this is reframing the bin decision (which is set by `f_tight` alone, per the pre-registered rule). It is the kind of side-observation that the brief invited.

---

## What the decision rule implies for Paper 200 v0.3 (per Task Z brief, unchanged in Z')

Bin MID branch:

> Paper 200 v0.3 abstract and §1 load-bearing-claim block must be rewritten to acknowledge that Crossing 1 is barely discriminating against the regime-density alternative, and Paper 200 either finds a second independent anchor (out of scope for v0.3) or scope-limits the claim correspondingly.

This is now binding. The choice between *second independent anchor* and *scope-limit* is CinC's, not Mr Code's.

---

## Cross-link to Task Z (audit-trail summary)

Task Z (commit `61aff82`) ran a brief in which the 5-series running sum was defined over all primes rather than primes ≡ 5 (mod 6). Under that wrong convention `S_5(31) = −6`, not 0. The locked script halted at the anchor check exactly as the brief instructed, and the halt was committed as the failed-brief audit record (commit `85d4ab8`). Task Z' (this task) was issued with the corrected convention; the audit-gap clock reset; the implementation reproduces both CinC's hand-computation and Paper 199 §6.1's figure.

Reading the audit trail from outside, the sequence is:

1. Brief was wrong (CinC error).
2. Pre-registration discipline forced the locked script to halt at the anchor.
3. Halt was recorded faithfully; no methodology drift to make the wrong brief "work".
4. Brief was corrected.
5. Re-run under the corrected convention; result confirms the prior.

The discipline did exactly the job it exists to do.

---

## Audit-gap record

- Task Z' pre-registration (this brief + locked script): `4ae5b36`, T_A = 2026-05-15T19:09:38+01:00.
- Task Z' results commit (this findings + results files): T_B at commit time, ≥30 min after T_A.

---

## Flags

1. **PRIMARY (decision):** `f_tight = 14/33 = 0.4242` → **Bin MID**. Paper 200 v0.3 abstract and §1 must be rewritten per brief; either second-anchor or scope-limit choice is CinC's.

2. **SUPPORTING:** Paper 199 §6.1's tie-density figure reproduces exactly under direct computation (14 ties at positions {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}). No §6.1 audit needed; the citation Paper 200 §5 leans on is sound.

3. **MINOR (worth mentioning to CinC, not load-bearing):** `f_post = 0.1515 < f_asymp = 0.1866`. The post-regime is, by this sample, slightly less zero-dense than the asymptotic regime — i.e. the tight regime's elevated density does not appear to taper monotonically into the asymptotic regime. Small-N caveat (post and tight are both length-33 samples; the asymptotic regime is length-134). Not a contradiction of anything in Paper 199; just a texture note.

4. **PROCESS (not for the paper, for the recognition log):** Task Z's wrong-brief halt and Task Z'`s clean re-run together demonstrate that the pre-registration + anchor-verification + audit-gap discipline catches CinC-side brief errors as cleanly as it catches implementer-side bugs. The discipline is symmetric.

🐕☕⬡

Mr Code, in band.
