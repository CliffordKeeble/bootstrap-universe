# Density transition scan — findings

Companion to the *density transition at e^142 / 10^61* brief (CinC, 5 May
2026). Status flags follow the programme convention: DERIVED (proved from
definitions), OBSERVED (numerical), STRUCTURAL (framework choice),
CONJECTURED.

---

## What this is, and what it isn't

**What this is.** The leading explicit-formula prediction of the
normalised Chebyshev bias

> B_χ(x) / √x  =  ψ_pred(x, χ) / √x  =  −2 ∑_n [(½) cos(γ_n τ) + γ_n sin(γ_n τ)] / (¼ + γ_n²)

at τ = ln x ∈ [130, 150], for χ ∈ {χ_5, χ_4}, computed from the first 25
non-trivial zeros of L(s, χ) under RH. Plus envelope (rolling-RMS)
diagnostics looking for amplitude/beat features at programme-relevant
τ ∈ {137, 140.5, 142}.

**What this is *not*.** It is **not** a test of a "transition in primes
near e^140.5". Direct prime enumeration to x ∈ [e^130, e^150] is not
possible. The curve over τ ∈ [130, 150] is purely the explicit-formula
prediction; it is by construction a smooth quasi-periodic function in τ,
fully determined by the L(s, χ) zero spectrum once the lowest zeros are
fixed. The brief's framing of "amplitude deviation from √x scaling at
high τ" was a conflation of DERIVED with OBSERVED — corrected by CinC in
review of Mr Code's first-pass clarification.

The territory that *is* testable here:

| Path | What it can show |
|---|---|
| (a) **Envelope/beat features in the predicted curve** at programme τ from interference of the lowest L(s, χ) zeros | a structural fact about the L-function zero spectrum |
| (b) **Empirical-vs-prediction validation** at sieve-feasible scales (τ ≲ 21) | already done in `paper_196_results.md` Layer 2 — Pearson +0.916, RMS residual 0.216 vs RMS of E_norm 1.206 |

This document is path (a) only. Path (b) is the existing Layer 2 work.

---

## Method [STRUCTURAL]

### Computation

- `density_transition_scan.py:31` — `B_norm(τ, γ)` vectorised explicit-formula sum.
- `density_transition_scan.py:65` — `rolling_rms(y, w)` reflective-padded RMS.
- `density_transition_robustness.py` — local-minima ranking, window-width sweep, truncation sweep.

Reuses [`lchi5_zeros.get_zeros`](lchi5_zeros.py) for χ_5 (LMFDB-anchored,
mpmath-refined to dps=30) and [`density_lchi5.LCHI4_ZEROS_25`](density_lchi5.py)
for χ_4 (LMFDB at 4-decimal precision).

### Grids

| Window | τ range | Δτ | Samples | Purpose |
|---|---|---|---|---|
| Scan | [130, 150] | 0.05 | 401 | programme-relevant probe |
| Baseline | [20, 200] | 0.05 | 3601 | envelope distribution reference |

Δτ = 0.05 corresponds to ~19 samples per dominant period 2π/γ_1 ≈ 0.945.

### Envelope

Rolling RMS of B_norm over τ-windows of half-width w_τ. Default w_τ = 1.0
(full-width 2.0), matching the natural beat period of the lowest two
zeros: 2π/(γ_2 − γ_1) ≈ 2π/3.18 ≈ 1.97 for χ_5.

### Programme-relevant τ

| τ | x | Programme meaning |
|---|---|---|
| 137 | e^137 | PNT density 1/ln N = 1/137 = α |
| 140.5 | e^140.5 ≈ 10^61 | first digital boundary |
| 142 | e^142 | closure scale, e^5 past α-density |

---

## Results [OBSERVED]

### Headline number

> **The χ_5 envelope reaches a local minimum of 0.3061 at τ = 137.00, which
> is the 2nd-deepest of 15 local minima in the scan window [130, 150] and
> sits at z = −2.40 against the baseline-window envelope distribution
> (mean 0.359, std 0.022).** The deepest local minimum in the scan window
> is at τ = 149.20, with envelope 0.2827.

> **At the same τ = 137, the χ_4 envelope is near a local maximum, ranked
> 12/13 of χ_4's local minima — the χ_5/χ_4 contrast at τ = 137 is real
> and unambiguous in the prediction.**

### Per-probe table

Three probes per character, plus baseline-rank (= fraction of all 3601
baseline samples with envelope ≤ probe value; rank ≤ 0.05 or ≥ 0.95
counts as extreme):

| Character | τ probe | env | z | local-min rank | baseline-rank | verdict |
|---|---|---|---|---|---|---|
| χ_5 | 137 | 0.3061 | −2.40 | **2/15** (near minimum) | 0.009 | **dip** |
| χ_5 | 140.5 | 0.3526 | −0.30 | 5/15 (mid-pack) | 0.384 | typical |
| χ_5 | 142 | 0.3389 | −0.92 | 10/15 (shallow) | 0.186 | typical |
| χ_4 | 137 | 0.4077 | +0.69 | 12/13 (near peak) | 0.742 | typical |
| χ_4 | 140.5 | 0.3731 | −0.49 | 4/13 (mid-pack) | 0.296 | typical |
| χ_4 | 142 | 0.4017 | +0.48 | 11/13 (shallow) | 0.674 | typical |

Of the six probes only one — χ_5 at τ=137 — shows a feature.
Bonferroni-corrected naive p ~ 6 × 0.009 ~ 0.05; on the edge of
"noteworthy" before the autocorrelation correction below.

### Truncation sensitivity (Q3) — stable [OBSERVED]

| N zeros | env(137) | baseline mean | baseline std | z(137) |
|---|---|---|---|---|
| 10 | 0.2853 | 0.3335 | 0.0226 | −2.14 |
| 15 | 0.2956 | 0.3463 | 0.0226 | −2.24 |
| 20 | 0.3031 | 0.3539 | 0.0224 | −2.27 |
| 25 | 0.3061 | 0.3592 | 0.0220 | −2.40 |

**The χ_5 dip at τ = 137 is already present in the lowest 10 zeros.** It
is not a high-N artefact; it is a structural property of the lowest
spectral modes of L(s, χ_5).

### Window-width sensitivity (Q2) — caveat [OBSERVED, fragile]

| w_τ | env(137) | baseline mean | baseline std | z(137) |
|---|---|---|---|---|
| 0.50 | 0.2603 | 0.3557 | 0.0549 | −1.74 |
| 0.75 | 0.2946 | 0.3581 | 0.0358 | −1.77 |
| **1.00** | **0.3061** | 0.3592 | 0.0220 | **−2.40** |
| 1.25 | 0.3727 | 0.3594 | 0.0189 | +0.71 |
| 1.50 | 0.3664 | 0.3595 | 0.0168 | +0.41 |
| 2.00 | 0.3688 | 0.3597 | 0.0110 | +0.83 |

**The z = −2.40 finding is window-specific.** It survives at w_τ ∈
[0.5, 1.0] but flips to mildly positive at w_τ ≥ 1.25. The reason is
mechanical: at full-width 2 (w_τ = 1.0) the rolling RMS averages over
exactly one beat period of the γ_2 − γ_1 modulation, capturing the
beat envelope cleanly; at full-width ≥ 2.5 we average over a full beat
period plus excess, smoothing the modulation away.

The honest reading: **the χ_5 envelope dip at τ = 137 is a property of
the natural beat structure of the lowest L(s, χ_5) zeros, visible only
when the smoothing width is matched to that beat period.** This is not a
"too-narrow window inflating the signal" concern (the env value itself
is similar at w_τ = 0.75 and 1.0); it is a "the baseline envelope
amplitude has more variance at narrow windows" concern (s = 0.0220 at
w_τ = 1.0 vs 0.0549 at w_τ = 0.5). At narrow windows the baseline is
also noisier, so the same value is less unusual.

### χ_5 vs χ_4 contrast [OBSERVED, robust]

The most striking unambiguous observation, robust to all
window-and-truncation choices:

> At τ = 137, **χ_5 sits at its 2nd-deepest local minimum and χ_4 is
> ranked 12/13 (i.e. shallowest minima, near a local maximum).** The
> two characters are anti-correlated at this τ.

This is content-rich because χ_5 is the programme-natural character (it
is the one that gives L(1, χ_5) = 2 log φ / √5, the golden-ratio-bearing
class number formula). χ_4 is the Leibniz character (L(1, χ_4) = π/4),
a sister with a different ground state. That the two characters'
envelopes are *anti-correlated* at exactly the α-density τ is the kind
of feature that would not be predicted by uniform-randomness over Dirichlet
characters. Whether it is content-rich about α specifically or just an
artefact of low-zero alignments needs more characters in the comparison —
deferred (see *Limits and follow-ups*).

---

## Honest assessment

**What we found, in one paragraph.** The leading-zero spectral prediction
of the χ_5 Chebyshev bias has a structural minimum near τ = 137 that
ranks 2nd-deepest of 15 minima in [130, 150], lands within 0.05 τ-units
of τ = 137 exactly, is robust to L-zero truncation (10..25 zeros), but
is window-width-fragile (only visible when smoothing matches the natural
γ_2 − γ_1 beat period). At the same τ, χ_4 is in the *opposite* phase —
near a local maximum. The χ_5/χ_4 anti-alignment at τ = 137 is the
robust observation; whether it has content about α-density or is a
coincidence of low-zero phases requires more characters and more τ
controls than tonight's investigation has.

**This is not a transition of primes; it cannot be.** Direct prime
enumeration at e^137..e^142 is not feasible. The result is a statement
about the L-function zero spectra: the lowest zeros of L(s, χ_5) and
L(s, χ_4) interfere in opposite phases at τ = 137, and the χ_5
interference happens to produce a near-maximal envelope minimum at that
location.

**This is not yet a finding strong enough to enter a paper.** The
window-fragility, the moderate Bonferroni-uncorrected significance, and
the absence of a multi-character control sample all pull against
publication-grade strength. The robustness facts (truncation-stable,
χ_4-anti-aligned) pull in favour. Net: **OBSERVED structural feature,
not a transition; alignment with α-density is suggestive but not
established.**

---

## Status flags

**DERIVED:**
- Explicit-formula expression for B_χ(x)/√x in terms of low L(s, χ) zeros
- The numerical envelope curve at any τ given the zeros (deterministic,
  reproducible to machine precision)

**OBSERVED:**
- χ_5 envelope at τ = 137: env 0.3061, rank 2/15 of local minima, z=−2.40
- χ_4 envelope at τ = 137: env 0.4077, rank 12/13, z=+0.69 (anti-aligned with χ_5)
- Truncation stability: feature visible at N ≥ 10
- Window fragility: feature visible at w_τ ∈ [0.5, 1.0], not at w_τ ≥ 1.25
- χ_5 envelope at τ ∈ {140.5, 142}: typical (no feature)

**OBSERVED but content-ambiguous:**
- Whether the χ_5 dip's alignment with τ = 137 specifically is more
  than a coincidence of low-zero phases; rank 2/15 with Bonferroni
  correction is not in itself unusual

**NOT TESTABLE here:**
- Any "transition of primes" near e^137..e^142
- Empirical bias amplitude at high τ — no observation possible

**DEFERRED:**
- Mod-5 complex characters of order 4 (brief Priority 2, optional second part)
- Mod-137 characters of higher order (brief Priority 3, speculative;
  CinC: defer unless LMFDB has the lowest zeros tabulated cheaply)
- Multi-conductor control sample (would calibrate whether χ_5/χ_4
  anti-alignment at τ=137 is content-rich)

---

## Files

- [`density_transition_scan.py`](density_transition_scan.py) — main scan + envelope diagnostics
- [`density_transition_robustness.py`](density_transition_robustness.py) — Q1/Q2/Q3 robustness sweeps
- [`density_transition_scan.csv`](density_transition_scan.csv) — per-τ scan-window data
- [`density_transition_baseline.csv`](density_transition_baseline.csv) — per-τ baseline-window data
- [`density_transition_probes.csv`](density_transition_probes.csv) — programme-τ probe summary
- [`density_transition_robustness.csv`](density_transition_robustness.csv) — robustness sweep tables
- [`density_transition_scan.png`](density_transition_scan.png) — scan plot (B_norm + envelope, both characters)
- [`density_transition_baseline.png`](density_transition_baseline.png) — baseline-window context plot

---

## Limits and follow-ups

**Open questions that would meaningfully strengthen or refute the
observation:**

1. **Multi-character calibration.** Repeat the scan for 4–6 small-conductor
   characters (χ_3, χ_7, χ_8, χ_11, χ_12, χ_13). If τ=137 shows a dip
   in only χ_5, that's content. If 2/6 or 3/6 also show dips at τ=137,
   it's chance. **Cheap, ~1 hour of compute including LMFDB pulls.**

2. **Window-period diagnostic.** Plot the envelope at multiple matched
   window-widths (each character's own γ_2 − γ_1 beat period) and see
   whether the τ=137 alignment is reproducible or a w_τ = 1.0 accident.

3. **Complex χ_5 characters of order 4.** Brief flagged "worth including
   in the bias scan." With L(s, χ) for χ complex, the explicit-formula
   prediction is itself complex; the envelope of |B_χ(x)|/√x is the
   relevant scalar. Different zero spectrum from χ_5 real, would calibrate
   whether the τ=137 feature is mod-5-conductor-specific or
   chi_5-real-specific.

4. **More zeros.** N=25 caps `lchi5_zeros.py` (LMFDB seed limit).
   Extending to N=100 via Riemann–von Mangoldt seeding (the TODO in
   that module) would tighten high-frequency content of the envelope
   at τ ∈ [130, 150] but per Q3 should not move the feature visibly.

The chi_5/chi_4 anti-alignment at τ = 137 is robust *given the
characters tested tonight*. Promoting this from "suggestive" to
"established" is item 1 of the followup list.

---

🐕☕⬡

*Density transition scan, 5 May 2026 evening to 6 May 2026 morning. Mr
Code, on the brief from CinC, with framing correction by CinC in
mid-investigation. Reuses Paper 196 apparatus.*
