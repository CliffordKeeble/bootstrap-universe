# Paper 157 Rigour Test v0.2 — Gate Verdict

**Pre-registration:** `paper_157_rigour_prereg_v0_2.md` (frozen, this commit).
**Script:** `gate_evaluation.py` · **Output:** `gate_output.txt`.
**v0.1 scoping input:** commit 0c79bc6 (`SCOPING.md`).

---

## VERDICT — NULL BY CONSTRUCTION

**STATUS: NULL (gate halt).** No candidate weighting rule (W1/W2/W3)
yields a **forced, peek-free** block weighting at the scale the residual
requires. Per brief v0.2 §2/§6, the §2 forcing gate halts the test before
§3 is reached. This is the brief's anticipated *safe* outcome (§6, row 3):
the 0.031% residual is the honest precision ceiling of a convergence
result, not a forced block-weight correction.

**157 consequence:** stays **OBSERVED** at the corrected 0.031%. The
convergence argument (Pattern 9) is untouched; no numerology trap is
acquired. Equipartition is confirmed the *wrong* mechanism (v0.1), and
block-weighting is confirmed *not forced* (v0.2) — so the residual has no
forced structural cause from the A₅ augmentation ideal alone.

---

## Why the gate is decisive — the scale obstruction

The CODATA-matching effective count is **N = 58.9815** (deficit 0.0185
from the uniform 59), a fractional perturbation of **3.13×10⁻⁴**. A
passing forced weighting must therefore land the four block weights at
**1 ± ~0.001**.

Every weighting built from A₅ representation theory — characters,
dimensions, or operator spectra — produces deviations of **O(1)**, three
orders of magnitude too large. Bridging that gap needs a fitted small
parameter, which is exactly the **four-knobs trap** the brief forbids
(§1, §4; the Paper 90 lesson: "constructible from icosahedral quantities"
≠ forced). This is a generic structural fact, not the result of comparing
candidate fits — all candidates fail for the same reason.

## Per-candidate gate evaluation (forcedness + scale)

| Candidate | Forced & peek-free? | Scale of weights | Gate |
|---|---|---|---|
| **W1** χ₅ character weighting | No — the χ₅→block contraction (Σ\|C\|χ_block(C)χ₅(ord)) gives raw weights **{18, 18, −16, 10}**; turning these into multiplicative weights ~1 needs a free shift+scale (two knobs). | O(1), one negative | **FAIL** |
| **W2** dimension / 2nd augmentation | No — integer mode removal gives N=58 (≈1.7%, ~55× the gap, wrong side). A forced fractional deficit = 0.0185 has no clean A₅/augmentation form (1/59=0.0170 is 9% off; 1/54 fits but 54 is not a tower quantity → reverse-engineered = peeking). | integer steps / none | **FAIL** |
| **W3** spectral-gap eigenvalue ratios | No — class-sum eigenvalues per block are O(1) integers/algebraics incl. 0 and negatives (e.g. 2A: {−5,−5,0,3}); forming weights needs free shift+scale, and the generating-set/operator choice is itself unforced (2A vs 3A vs 5A give different spectra). | O(1) | **FAIL** |

Note (W1 detail): χ₅(ord) annihilates the 5-fold classes (χ₅(5)=0), so
the two 12-element vertex classes drop out of the contraction entirely —
a structurally interesting fact, but it does not rescue the scale.

## Discipline notes

- The gate was evaluated on **forcedness and scale only**; **N_eff was
  not computed against the ±0.003 window** for any candidate (that is §3,
  unreached). Showing all three candidates fail the gate documents the
  null — it is not a multiple-comparison search for a winner (brief §3.6).
- Single pass, frozen spec, no re-spec to force a result.
- The §5 arithmetic rein-back of Paper 157 is reported separately and
  stands regardless of this verdict — see `REINBACK_157.md`.
