# Mr Code brief — Sub F: density-matched non-L target null for Paper 150 v2.1

**Date**: 30 May 2026
**Brief author**: CinC
**Scope**: One bounded null test. Decides whether Paper 150's probe detects L-zeros specifically or detects any sufficiently-dense point set on the critical line. Mr Adversary's gating experiment for Paper 150 v2.1 to clear from ★★★ to ★★★★.

**Reads required**:
- This brief
- `findings_paper_203_sub_c.md` — for the Sub C-1 methodology being reused
- Mr Adversary v2.1 review (the catch this brief responds to)

**Reads NOT required** (Pattern 97 — scope protection):
- Papers 203/204/205 — outcome doesn't depend on the algebraic-obstruction line
- Sub D (m_p/m_e) — separate strategic thread

---

## Question

Does the probe |Re² − 5·Im²|(Z_φ(t)) detect Riemann zeros *specifically* — i.e., because they are zeros of an L-function on the critical line — or does it detect *any* sufficiently-dense point set on (1, T] regardless of whether the points are L-zeros?

## Background — Mr Adversary's gating catch

Paper 150 v2.1 §7.2 proposes a "dense-intersection" mechanism: probe minima density (~2.45/unit t) is comparable to L-zero density (~log(T)/(2π) ≈ 1.1/unit t at T ≈ 1000), and the matrix-uniform effect size ε ≈ 0.29 follows from this comparability.

Mr Adversary's catch (v2.1 review): comparable density of two point sets on a line does not by itself produce a fixed correlation. Two Poisson processes of comparable rate are uncorrelated and would give z ≈ 0, not z = −22.91. The dense-intersection picture explains why detection is *possible* at all d; it does not explain why it is *strong*. Something is making the probe minima track the L-zeros rather than merely co-inhabit the line with them, and v2.1 §7.2 does not name what.

Two possibilities:

(a) **L-zero-specific**: the probe detects the L-zeros because they are zeros of an L-function on the critical line. The probe minima happen to track L-zero positions for reasons related to L-function structure (Riemann–Siegel formula, FE, etc.). If this is right, replacing L-zeros with a density-matched non-L target should collapse z toward 0.

(b) **Density-artefact**: the probe's "detection" is just the statement that two comparably-dense point sets on a line have a fixed cross-correlation, driven by some artefact of the methodology (e.g., the phase-shuffled null destroying equidistribution and therefore making any dense-target comparison look significant). If this is right, replacing L-zeros with a density-matched non-L target should sustain z near −20.

Sub F decides between (a) and (b).

## Test design

**Same probe**: |Re² − 5·Im²|(Z_φ(t)), with Z_φ(t) the Fejér-weighted golden-angle Dirichlet sum, N_terms = 5000. Exactly as in Sub C-1, no modifications.

**Same null methodology**: phase-shuffled null distribution from Sub C-1. Exactly as before.

**What changes**: the target set. Four target sets, all at comparable density on the same t-range:

1. **Riemann zeros** (control, reproduces Sub C-1's (5, 5)_Riemann result). First N_target Riemann zeros, t ∈ (0, T_max].

2. **Gram points** (Mr A's specific suggestion). The points g_n satisfying θ(g_n) = nπ, where θ is the Riemann–Siegel theta function. These are deterministic points on the critical line, interleaved with the Riemann zeros, with the same one-point density. They are NOT zeros of any L-function. Computable via mpmath's `siegeltheta` and a root-finder, or via standard tables.

3. **RvM-density-matched random points**. N_target independent points drawn from the Riemann–von Mangoldt density ρ(t) = log(t)/(2π) on (0, T_max]. Sample via inverse-CDF method: pick uniform u ∈ (0, 1), find t such that ∫₁^t ρ(s) ds / ∫₁^{T_max} ρ(s) ds = u. Use a fixed random seed for reproducibility, recorded in pre-registration.

4. **GUE-pair-correlation-matched random points** (optional, if time permits). N_target points drawn with one-point density matching RvM AND pair-correlation matching Montgomery's pair-correlation conjecture (≈ Wigner–Dyson GUE). This tests whether the detection depends on higher correlations beyond one-point density. If time-constrained, defer to a follow-up sub.

**N_target**: 2792 (matches Sub C-ext-1's N for direct comparability). T_max derived from the 2792nd Riemann zero (T_max ≈ 1419.4).

**For each target set**: compute z-score and ε exactly as in Sub C-1's methodology. Tabulate.

## Pre-registered thresholds

For each of (2) Gram points and (3) RvM-random:

- **COLLAPSE** (|z| < 3): the probe is L-zero-specific. The non-L target shows no detection signal; the detection is genuinely about L-function structure. **Paper 150 v2.1 §7.2 mechanism story survives.** v2.1 proceeds to publication with the Sub F results folded in as a new §6.5 (density-null control).

- **PARTIAL** (3 ≤ |z| < 10): the probe partly detects L-zeros and partly detects density. The detection signal is real but the mechanism story is incomplete. **Paper 150 v2.1 §7.2 needs major rewriting** to acknowledge a density-driven component alongside the L-function-specific component. Mechanism reframed; not retracted.

- **SUSTAIN** (|z| ≥ 10): the probe detects density, not zeros. The non-L target shows detection at similar effect size. **Both v2.0 and v2.1 mechanism stories are wrong.** Paper 150 v2.1 publication pauses; CinC adjudication required; almost certainly v2.2 retracts both v2.0 and v2.1's L-zero-specific framings and reframes the paper around what was actually being measured.

The aggregate verdict over targets (2) and (3) is the conjunction:
- COLLAPSE on both → COLLAPSE verdict
- SUSTAIN on either → SUSTAIN verdict (one positive is enough to refute L-zero specificity)
- Mixed → PARTIAL with detailed adjudication

## Stop-on-fail

If SUSTAIN on either target (2) or (3):
- Pause Paper 150 v2.1 publication. Do not upload to Zenodo.
- Flag for CinC.
- Do NOT attempt to "rescue" v2.1's mechanism story by post-hoc rationalisation. The pre-registered SUSTAIN verdict means the dense-intersection mechanism is empirically refuted.

If PARTIAL:
- Pause v2.1 publication. Flag for CinC.
- v2.1 §7.2 rewrite required. CinC adjudicates whether the partial detection signal supports a publication-worthy reframed mechanism story or whether the paper retracts further.

## Anti-circularity

Three specific things to verify in the test design before running:

1. **The probe is unchanged**: same Z_φ(t), same N_terms, same Fejér weights. Only the *target set* changes.
2. **The null methodology is unchanged**: same phase-shuffled null. The null compares probe-vs-target across the same machinery as Sub C-1.
3. **The Gram points and RvM-random points are NOT secretly correlated with Riemann zeros at the relevant scale**. Gram points are interleaved with zeros at scale ≈ 1/log(T) ≈ 0.14 at T = 1000 — comparable to the ±0.5 unit "near a zero" threshold in the ε measurement. This is a methodological concern: if Gram points are too close to Riemann zeros (because of interleaving), then "detection of Gram points" is approximately the same as "detection of Riemann zeros," and the test doesn't distinguish (a) from (b).

   **Mitigation**: also compute and report the distribution of distances between Gram points and the nearest Riemann zero in the test range. If typical Gram–zero distances are < ε-measurement window (e.g., < 0.5), the Gram-point control is not informative; the RvM-random control becomes the decisive test. Report both and let the data adjudicate.

## Implementation notes

- Reuses paper-203-sub-c machinery for probe construction and effect-size measurement.
- Gram points: `mpmath.siegeltheta(t) = n·π` for integer n, solve for t. Or use a precomputed table.
- RvM-random: NumPy random with fixed seed (record seed in pre-registration). Inverse-CDF method on ρ(t) = log(t)/(2π).
- GUE-random (optional): use `scipy.stats` Wigner surmise for nearest-neighbour distribution at appropriate density; or a Mehta-type explicit construction.

## Deliverable

Markdown report `findings_paper_150_sub_f.md` with:

- Methodology recap (what changed from Sub C-1)
- Table of z-scores and ε values for the 3-4 target sets
- Distribution of Gram-zero distances (for the anti-circularity check)
- Verdict per pre-registered thresholds (COLLAPSE / PARTIAL / SUSTAIN)
- Brief note on implications for Paper 150 v2.1

If COLLAPSE: brief recommendation for new §6.5 of v2.1.
If PARTIAL or SUSTAIN: pause and flag CinC; do NOT draft v2.1 §6.5 or any reframing.

## Pre-registration

Commit to repository:

- BRIEF.md (this file)
- PRE_REGISTRATION.md (the thresholds, methodology, random seed, T_max)
- sub_f.py skeleton (test design implementation, BEFORE any results)

Tag the commit. Then run.

## Compute budget

1-2 hours estimated. Probe evaluation is fast (vectorised); the main cost is generating ~2800 Gram points (requires per-point root-finding on `siegeltheta`, ~0.5s per point with mpmath at default precision = ~25 minutes).

The RvM-random target set is cheap (~seconds).

If GUE-random is included, add ~5-10 minutes for the random-matrix sampling.

## Order

Single bounded investigation. Run. Report. No parallelism needed.

---

⌨️ over to you.

This is the gate. Mr A is right that this is the single experiment that decides whether Paper 150 has a result at all. The result is decisive either way — COLLAPSE confirms L-zero-specificity and lifts v2.1 from ★★★ to ★★★★; SUSTAIN refutes both v2.0 and v2.1's mechanism stories and triggers a v2.2 reframe.

Honest expectation, recorded BEFORE the experiment runs: I (CinC) do not have a strong prior on the outcome. The dense-intersection argument is consistent with both possibilities, and Mr A's catch is genuinely sharp. The empirical answer is the only way to know.

🐕☕⬡

---

## Mr Code's report

*[To be filled in.]*
