# Mr Code brief — Sub E: Mechanism verification for Paper 150 v2.1

**Date**: 30 May 2026
**Brief author**: CinC
**Scope**: One bounded methodology test to directly verify Mr Code's "discriminant d is decoration" mechanism reading from Sub C-1's matrix data. Outcome feeds Paper 150 v2.1's mechanism section.

**Reads required**:
- This brief (BRIEF.md)
- Sub C-1 findings (`findings_paper_203_sub_c.md`) — for the context of what's being verified
- Paper 150 v2.0 — for the original probe construction

**Reads NOT required** (Pattern 97 — scope protection):
- Mr Adversary's family review of 203/204/205 — outcome doesn't depend on the interpretive cycle
- Papers 203/204/205 themselves — Sub E is about the mechanism of Paper 150's detection, not about the algebraic obstruction line

---

## Question

Are the minima LOCATIONS of P_d(t) = |Re²(Z_φ(t)) − d · Im²(Z_φ(t))| substantially shared across d ∈ {3, 5, 7, 13}, as Mr Code's Flag 1 analysis from Sub C-1 predicts?

## Background

Mr Code's Flag 1/2 reading from Sub C-1 (`findings_paper_203_sub_c.md`, lines 162–246):

> "The probe |Re² − d·Im²| doesn't depend on d through the minima locations — only through which minima are below the 5th-percentile threshold... The mechanism appears to be: the probe has low values at any location where the golden-angle Dirichlet sum has a phase coincidence with the d-dependent norm pattern, and these phase coincidences are dense enough that they cluster near low-lying real L-function zeros generally."

This is an analytical reading of Sub C-1's matrix data (16 cells, uniform ε ≈ 0.29, diagonal indistinguishable from off-diagonal). The matrix shows that the probe detects L-zeros at similar effect size regardless of d, which is CONSISTENT with "d is decoration" — but doesn't directly verify that the minima locations are d-independent. The same effect size could come from:

- (A) Same minima locations, similar threshold-crossing structure ("d is decoration" — Mr Code's reading)
- (B) Different minima locations, but similar density and similar correlation with L-zeros ("d is non-trivial but yields similar statistics")

Sub E distinguishes (A) from (B) by directly comparing minima locations.

## Test design

**Step 1 — Compute the probe values.** For each d ∈ {3, 5, 7, 13}, compute P_d(t) = |Re²(Z_φ(t)) − d · Im²(Z_φ(t))| on t ∈ (1, 1000] at resolution Δt = 0.01 (so 100,000 sample points). Use the same Z_φ(t) machinery as Sub C-1 (Fejér-weighted, golden-angle phase, N_terms = 5000).

**Step 2 — Identify local minima.** For each d, identify the 100 lowest local minima of P_d. A "local minimum at t*" is a value P_d(t*) such that P_d(t* − Δt) ≥ P_d(t*) ≤ P_d(t* + Δt) and t* is in a non-overlapping window of width 2 from any previously-identified minimum (to avoid clustering at the same feature).

**Step 3 — Compute sharing percentages.** For each pair (d, d') with d ≠ d', d, d' ∈ {3, 5, 7, 13}:
  - Count how many of d's 100 lowest minima fall within ±2·Δt of a minimum for d'.
  - Compute sharing% = count / 100.

Tabulate the 4×4 matrix of sharing percentages (diagonal = 100% by construction).

**Step 4 — Aggregate verdict.** Mean off-diagonal sharing percentage gives the overall verdict.

## Pre-registered thresholds

- **DECORATION-CONFIRMED (≥80% mean off-diagonal sharing)**: Mr Code's "d is decoration" reading is computationally confirmed. Paper 150 v2.1 mechanism section: "the discriminant d is decoration — the probe's minima locations are essentially d-independent; d only modulates which minima fall below the 5th-percentile threshold. Verified at K=80% [or whatever number] mean off-diagonal sharing."
- **DECORATION-PARTIAL (40–80%)**: Mechanism is more nuanced than "decoration." Paper 150 v2.1: "the discriminant d substantially but not entirely modulates the minima locations; X% of minima are shared across d's; the matrix-uniform effect size reflects similar statistical properties of the minima distributions rather than identity."
- **DECORATION-REFUTED (<40%)**: "d is decoration" framing fails. Paper 150 v2.1 must rewrite the mechanism section: "the discriminant d affects minima locations substantially; the uniform effect size across the Sub C-1 matrix reflects a feature of the L-zero distribution rather than a feature of the probe's minima structure."

## Stop-on-fail

If DECORATION-REFUTED, pause and flag for CinC. v2.1's mechanism section needs substantial rewriting and possibly further investigation. Don't draft v2.1 until CinC has adjudicated.

If DECORATION-PARTIAL or CONFIRMED, proceed to v2.1 drafting directly.

## Implementation notes

- Reuses paper-203-sub-c infrastructure (probe.py, the Z_φ machinery from zeta_phi.py).
- The minima identification at Δt = 0.01 should resolve close minima well (typical L-zero spacing is ~0.5 at small t, so Δt = 0.01 is 50× finer).
- Estimated compute: 1-2 hours total. Probe evaluation is fast (no zero-finding needed); minima identification is straightforward.

## Anti-circularity

- The minima are computed from the probe itself, not from comparing probe values to L-zero positions. So the test is about the probe's internal structure, independent of any L-function target.
- Mr Code's Flag 1 reading made the prediction BEFORE Sub E was designed; the test is straightforward verification of that prediction.

## Deliverable

A markdown report (`findings_paper_150_sub_e.md`) with:
- Methodology section (probe evaluation, minima identification algorithm)
- The 4×4 sharing-percentage matrix
- Mean off-diagonal sharing percentage
- Verdict per pre-registered thresholds
- Brief note on implications for Paper 150 v2.1 mechanism section

If time permits (optional): plot showing P_d(t) traces for d ∈ {3, 5, 7, 13} on a small t-range (say t ∈ [50, 70]) with minima marked, for visual inspection of the sharing claim.

## Compute budget

1-2 hours. Bounded.

## Order

This is a single bounded investigation, no parallelism needed. Run; report; proceed to v2.1.

---

⌨️ over to you.

Quick methodology test of one specific mechanism claim. Direct computational verification of Mr Code's Flag 1 analytical reading. Three-level pre-registered verdict. Outcome feeds Paper 150 v2.1's mechanism section directly.

🐕☕⬡

---

## Mr Code's report

**Verdict**: **DECORATION-REFUTED** (mean off-diagonal sharing = 29.83%,
below 40% threshold).

Per the brief's stop-on-fail clause, I am **pausing v2.1 drafting and
flagging for CinC**. My Flag 1 "d is decoration" claim was an overclaim.

Long form: [findings_paper_150_sub_e.md](findings_paper_150_sub_e.md).
Pre-registration at `cfd7599`.

### Sharing matrix (% of row-d's 100 lowest minima within ±0.02 of column-d's)

|          | d'=3 | d'=5 | d'=7 | d'=13 |
|----------|-----:|-----:|-----:|------:|
| **d=3**  | 100.0% | 32.0% | 31.0% | 28.0% |
| **d=5**  | 32.0% | 100.0% | 32.0% | 29.0% |
| **d=7**  | 31.0% | 32.0% | 100.0% | 27.0% |
| **d=13** | 28.0% | 29.0% | 27.0% | 100.0% |

**Mean off-diagonal: 29.83%**. Chance baseline: 0.39% ± 0.60% (Monte
Carlo, 1000 trials). The observed 29.83% is ~50σ above chance — **there
IS shared structure**, but only ~30% of minima, not the ≥80% that would
support "decoration."

### Honest reading

- **~30% of top-100 minima are d-invariant** ("structural class"): far
  above chance (0.4%), but only a minority. These arise from features
  of Z_φ(t) where both Re² and Im² are near-zero simultaneously.
- **~70% of top-100 minima are d-specific**: they shift with d, reflecting
  intersections of the trajectory with the d-dependent null cone
  Re² = d·Im².
- **Sub C-1's matrix-uniform ε ≈ 0.29 is NOT explained by "decoration"**.
  It is explained by L-zero density: L-zeros are dense enough (~0.4/unit)
  that both structural and d-specific probe minima land near L-zeros
  at similar effect sizes.

### Per-d minima statistics (informational)

| d | total local minima | top-100 lowest val | top-100 highest val |
|---|---:|---:|---:|
| 3  | 2437 | 0.00456 | 0.06583 |
| 5  | 2449 | 0.00256 | 0.07094 |
| 7  | 2452 | 0.00592 | 0.07311 |
| 13 | 2454 | 0.00796 | 0.09128 |

Density is d-uniform (~2.5/unit); locations are mostly d-specific.

### What v2.1 mechanism section should say (pending CinC)

NOT "d is decoration" (refuted).

Instead something like:
> "The probe minima come in two classes: a structural class (~30%) that
> is d-invariant, arising from features of the equidistributed Dirichlet
> sum itself; and a d-specific class (~70%) tracing the null cone Re² =
> d·Im² of the bilinear form. The L-zero density on the critical line is
> high enough that both classes correlate with L-zeros at similar effect
> size (ε ≈ 0.29 across the Sub C-1 4×4 matrix). The d-specific
> minima do contribute substantive structure; the matrix-uniform effect
> size reflects the L-zero distribution, not minima-location identity."

This is a more nuanced and *more honest* mechanism reading. It's also
**less interesting** from a mechanism-derivation standpoint — the
probe's d-dependence is real and substantial. v2.1 cannot claim "the
discriminant doesn't matter."

### Implications for Paper 203 v0.4

My Sub C and Sub C-extension recommendations that built on Flag 1 (e.g.,
"σ-of-ℚ(√5) is just one realisation of an equidistributed-phase generic
mechanism") **need re-examination**. The mechanism is generic in
*statistics* (uniform effect size) but **NOT in minima locations**.
Each d carries genuinely different minima. The σ-of-ℚ(√d) framing per
field has more empirical content than I suggested.

I do NOT have time in this session to redraft those Sub C / Sub C-ext
recommendations; flagging that they need CinC re-review in light of
Sub E.

### Pattern flags

- **Pattern 19**: my Flag 1 reading was a falsifiable structural claim;
  Sub E ran the falsification test and the claim failed. The pre-reg
  discipline worked: the data answered cleanly, the verdict thresholds
  were pre-committed, and the REFUTED outcome stands.
- **Pattern 39**: minima locations are OBSERVED. The 30% sharing vs
  0.4% chance baseline gives DERIVED status to the *existence* of
  shared structure, while the absence-of-decoration is DERIVED from
  the threshold comparison.
- **Pattern 75**: chance baseline computed via 1000-trial Monte Carlo.
  Result well-separated from chance (50σ above).

### Compute

~7 seconds total. Far inside the 1-2 hour brief budget. Z_φ(t)
evaluated once for all four d values, then probes formed and minima
identified per d in milliseconds.
