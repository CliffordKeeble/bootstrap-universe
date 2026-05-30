# Findings — Skeleton-Fixed Null for α⁻¹

**Investigation:** broader-family null over the skeleton `F = eᵃ − b√c − d + 1/n`
with icosahedral integers, holding the *functional form* fixed and varying the
five integer parameters.
**Serves:** Paper 2 v8.1 §4.1 (broader-family null discharge) and §10 (numerology
answer); Paper 127 §2 (retroactive skeleton-fixed support).
**Brief:** [`Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md`](Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md),
pre-registered 30 May 2026, seed 42, locked thresholds.
**Date run:** 30 May 2026.

---

## Question

Within the family `eᵃ − b√c − d + 1/n` over uniformly random integer parameters
in bounded ranges (a∈1–10, b,c,d∈1–30, n∈2–200; 53,730,000 tuples), what
fraction land within 50 ppm of 137.036? If decisively below the pre-registered
1% threshold, the canonical (5, 6, 3, 1, 66) choice is *structurally
distinguished within the skeleton*, not generically achievable.

## Method

Two scripts, run in sequence (the order is the point — pre-registered anchor
first, stronger layer second):

1. **`paper_002_s4_1_skeleton_null.py`** — the pre-registered Monte-Carlo run.
   N = 10⁶, seed 42, locked thresholds. Ran exactly as committed, **unedited**.
2. **`paper_002_s4_1_skeleton_null_exhaustive.py`** — confirmatory secondary,
   added *after* the committed run landed (so it is visibly post-result, not
   post-hoc tuning). Enumerates the **entire** 53.73M-tuple space exactly,
   removing seed dependence and making the tight 50 ppb window meaningful.

## Result

| Window | Accounting | Count | Fraction | Pre-reg threshold |
|---|---|---|---|---|
| 50 ppm (MC, N=10⁶, seed 42) | sampled | 8 | 8.0×10⁻⁶ | 1×10⁻² |
| 50 ppm (exhaustive) | raw tuples | 498 / 53.73M | 9.27×10⁻⁶ | 1×10⁻² |
| 50 ppm (exhaustive) | **distinct values** | **341 / 53.73M** | **6.35×10⁻⁶** | 1×10⁻² |
| 50 ppb (exhaustive) | raw tuples | 4 / 53.73M | 7.44×10⁻⁸ | — |
| 50 ppb (exhaustive) | **distinct values** | **2 / 53.73M** | **3.72×10⁻⁸** | — |

**Pre-registered decision — p < 1% → claim structurally distinguished — holds
decisively, now exact rather than estimated.** The MC headline (8.0×10⁻⁶) and
the exact fraction (9.27×10⁻⁶) agree. **STATUS: the broader-family null is
DISCHARGED at the skeleton-fixed level.**

---

## Three findings beyond the headline

### 1. Effective dimensionality — only a = 5 contributes (exact)

All 498 matches sit at **a = 5**; a ∈ {1,2,3,4,6,7,8,9,10} contribute **exactly
zero**. Structurally: e⁴ ≈ 54.6 cannot reach 137 even with the smallest
subtraction, and e⁶ ≈ 403 cannot descend to 137 even with the largest
(max b√c + d ≈ 194). So a is *effectively one-dimensional* in this null.

The honest figure to quote is therefore the **a = 5 conditional**:
**498 / 5,373,000 = 9.27×10⁻⁵** — still ~100× below the 1% threshold, i.e. two
orders of margin under the most demanding honest accounting. The unconditional
9.27×10⁻⁶ is the conditional diluted by the 1/10 prior on a, and should not be
quoted as if the extra order of magnitude were earned.

**STATUS: DERIVED** (exhaustive, exact).

### 2. The b√c term is degenerate — report distinct values, not raw tuples

`b√c` is many-to-one onto surd values: `6√3 = 3√12 = 2√27` are one number
written three ways. The search space ranges (b,c) over *all* integer pairs, so a
single surd is counted up to its representation multiplicity, inflating raw
tuple counts — worst at the tight end, where it matters most:

- 50 ppm: 498 raw tuples → **341 distinct values** (×1.46 inflation).
- 50 ppb: 4 raw tuples → **2 distinct values** — three of the four are the
  canonical `6√3` surd in disguise (the (5,6,3,1,66), (5,3,12,1,66),
  (5,2,27,1,66) rows all give F = 137.036005772 identically).

**Any count over a skeleton with a b√c term must report distinct values (or
restrict to square-free c).** The §4.1 figures use the distinct-value accounting.

**STATUS: STRUCTURAL** (representational degeneracy of the surd term).

### 3. The canonical formula is *not unique* at 50 ppb — one genuine competitor

After collapsing the degeneracy, **two distinct values** survive at 50 ppb:

| Tuple | Formula | F | vs TARGET | vs CODATA 2018 |
|---|---|---|---|---|
| (5, 6, 3, 1, 66) | e⁵ − 6√3 − 1 + 1/66 | 137.036005772 | +4.21×10⁻⁸ | **+4.9×10⁻⁸ (≈+49 ppb)** |
| (5, 1, 29, 6, 125) | e⁵ − √29 − 6 + 1/125 | 137.035994295 | +4.16×10⁻⁸ | **−3.5×10⁻⁸ (≈−35 ppb)** |

Against the *actual* CODATA value 137.035999084 (not the brief's rounded
TARGET = 137.036), the two occupants **straddle CODATA on opposite sides** —
canonical ~49 ppb above, competitor ~35 ppb below — both inside the 50 ppb
window. The competitor is real and lands fractionally *tighter* on the brief
target.

Its icosahedral reading is **strained**: 29 = E−1 is mild, 6 = D! is fine, but
125 = 5³ is an unusual icosahedral structural quantity and there is no clean
(D!, D, 1, D!·(V−1)) story as the canonical has. It exists *numerically* without
the canonical's *structural narrative*.

**Implication for the paper:** the skeleton is structurally distinguished (2
distinct values in 53.73M = 3.7×10⁻⁸), but the canonical tuple is **not uniquely
so at this precision**. §10 should read "structurally distinguished; one
numerical competitor at the same precision, lacking the canonical's icosahedral
structural story" — not "canonical uniquely matched." Naming the competitor
costs less than pretending it isn't there. §4.1's narrower four-member 1/n-family
uniqueness claim (with (a,b,c,d) = (5,6,3,1) fixed) is **unaffected** — a
different, narrower audit.

**STATUS: OBSERVED** (a numerical near-coincidence without derived structure).

---

## Decisions made beyond the brief

- **Exhaustive enumeration added** as a confirmatory secondary after the
  committed MC run (CinC-approved; the run-as-committed-then-improve order keeps
  pre-registration intact). The committed MC script was **not edited**.
- **Distinct-value collapse threshold:** round F to 1e-7 (50 ppm) / 1e-8
  (50 ppb). Surd-identical tuples collapse cleanly; no borderline cases.
- **CODATA column** added to the CSV and the 50 ppb table (CinC's refinement —
  the straddle is only visible against the true CODATA value, not the rounded
  target).

## Files

- [`Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md`](Mr_Code_Brief_Paper_002_section_4_1_skeleton_null_v0_1.md)
  — pre-registered brief (unedited since commitment).
- [`paper_002_s4_1_skeleton_null.py`](paper_002_s4_1_skeleton_null.py)
  — committed MC script (8/10⁶ at 50 ppm, seed 42).
- [`paper_002_s4_1_skeleton_null_exhaustive.py`](paper_002_s4_1_skeleton_null_exhaustive.py)
  — confirmatory exact enumeration; writes `matches_50ppm.csv`.
- [`matches_50ppm.csv`](matches_50ppm.csv) — all 498 matches with rel-to-target,
  rel-to-CODATA, 50 ppb flag, and distinct-value id.

## Lessons captured (for future skeleton-fixed briefs)

1. **Surd canonicalisation** — for any null with a b√c term, report distinct
   values or restrict to square-free c; representational degeneracy inflates
   counts at the tight end where it matters most.
2. **Effective-dimensionality conditional reporting** — when one parameter is
   effectively-dimensional (only some values can contribute), quote the
   conditional fraction, not the unconditional flattered by the dilution.
3. **Exhaustive when the space allows it** — 53.7M is feasible in seconds;
   MC + rule-of-three is for when it isn't. Pre-registered MC is the discipline
   anchor; exhaustive is the stronger claim layered on top.

🐕☕⬡
