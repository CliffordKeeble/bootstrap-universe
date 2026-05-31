# Pre-registration — Sub F (density-matched non-L target null)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Binding**: this document is committed to git **before** Gram points are
computed, RvM-random points are sampled, or the probe is run against any
non-Riemann target.

---

## Honest expectation, recorded BEFORE the experiment runs

CinC reports "no strong prior". I (Mr Code) record my prior here:

After Sub E's refutation of the Flag 1 "d is decoration" claim
(`paper-150-sub-e/findings_paper_150_sub_e.md`, commit `3056092`),
I take Sub C-1's matrix at less face value. The dense-intersection
argument (probe minima density ~ L-zero density ⇒ correlation possible)
is **necessary but not sufficient** for the observed z = −22.91. Mr A
is right.

My weak prior: I lean toward COLLAPSE on the RvM-random target (because
random points are genuinely uncorrelated with the probe minima), and
toward PARTIAL on the Gram-point target (because Gram points are
deterministic and interleaved with Riemann zeros at small scale,
which may make the test less discriminating). But this is genuinely
weak — I've been wrong recently (Flag 1), and I'm not banking on it.

The empirical answer is what matters.

## Construction commitments

**Probe**: `P(t) = sqrt(|Re²(Z_φ(t)) − 5·Im²(Z_φ(t))|)` with
N_terms = 5000, α = φ, Fejér σ_n = 1 − n/N. Identical to Sub C-1's
(5, 5) probe.

**Probe grid**: t ∈ (1, T_max + 5] at Δt = 0.008 (matches Sub C-1
`run_N10000.py` resolution).

**Minima detection**: 3-point local minima, 5th percentile threshold,
deduplication window 0.3. Identical to Sub C-1 `probe.find_minima`.

**Match window**: W = 1.0. Identical to Sub C-1.

**MC null**: N_null = 1000 trials, uniform random fake "minima" at
the same density as observed. Seed = 42 (matches Sub C-1).

**N_target = 2792**. T_max = t-coordinate of 2792nd Riemann zero
≈ 1419.4. Sub C-ext-1's N for direct comparability.

## Target sets (4 commitments)

### Target 1 — Riemann zeros (control)

First 2792 non-trivial Riemann zeta zeros (imaginary parts), pulled
via `mpmath.zetazero`. The Sub C-1 / Sub C-ext-1 (5, 5)_Riemann control.

I already have 1900 Riemann zeros cached at
`paper-203-sub-c/riemann_zeros_10000.csv`. Need to extend by ~890.

### Target 2 — Gram points

Gram points g_n are defined by:
> θ(g_n) = n·π for n = 0, 1, 2, ...

where θ is the Riemann–Siegel theta function. I solve via
`mpmath.findroot` with bracketing.

mpmath has `mp.siegeltheta(t)`. Compute g_n for n = 1, 2, …, ~2792.
First Gram point g_0 ≈ 17.85.

### Target 3 — RvM-density-matched random points

Riemann–von Mangoldt one-point density:
> ρ(t) = (1/2π) · log(t / (2π))

(more accurate than `log(t)/(2π)` for finite t; the constant matters
for matching precise density). For consistency with the brief I'll use
the brief's stated `ρ(t) = log(t)/(2π)`.

Sampling procedure (inverse-CDF):
1. Compute CDF F(t) = ∫₁^t (log(s)/(2π)) ds = (1/(2π)) · (t·log(t) − t − (1·log(1) − 1)) = (1/(2π))(t·log(t) − t + 1).
2. Normalise: F̃(t) = F(t) / F(T_max).
3. For 2792 uniform u ∈ (0, 1), solve F̃(t) = u for t numerically via bisection.
4. Sort the resulting t-values.

**Seed = 20260530** (pre-registered date).

### Target 4 — GUE-pair-correlation-matched (optional)

If time permits after Targets 1-3. Otherwise deferred.

Generate eigenvalues of a Gaussian Unitary Ensemble (GUE) matrix of
appropriate size such that the sample is scaled to match the
one-point density log(t)/(2π) on (1, T_max] AND has the GUE pair
correlation Σ₂(L) ≈ (1/π²)(log(L) + γ + 1) at small L.

I commit to attempting this only if Targets 1-3 take less than
1 hour total.

## Anti-circularity checks (mandatory)

Before running the probe match against each target, verify:

1. **Probe unchanged**: assert that `compute_probe(t, d=5, N=5000)` from
   `paper-203-sub-c/probe.py` is the function used; no reimplementation,
   no modifications.

2. **Null unchanged**: assert that `mc_null(...)` from
   `paper-203-sub-c/probe.py` is used; same seed, same procedure.

3. **Gram-Riemann interleaving distribution**: compute
   `min_distance(g_n, nearest_zeta_zero)` for all 2792 Gram points,
   report distribution. If median is below W=1.0, then Gram detection
   ≈ Riemann detection by interleaving, and the Gram-point test is
   not informative about (a) vs (b) — the RvM-random control becomes
   the decisive test.

   Pre-register: I will report the median and the quantile
   distribution (10th, 50th, 90th percentile). The verdict on the
   Gram-point target is contingent on this check.

## Pre-registered thresholds (re-stated from brief)

Per non-L target (Gram and RvM-random):

- COLLAPSE: |z| < 3
- PARTIAL: 3 ≤ |z| < 10
- SUSTAIN: |z| ≥ 10

Aggregate verdict over the two non-L targets:
- COLLAPSE on both → COLLAPSE verdict
- SUSTAIN on either → SUSTAIN verdict
- Mixed → PARTIAL with adjudication

**Stop-on-fail**:
- SUSTAIN on either: pause v2.1 publication, flag CinC, no rescue.
- PARTIAL: pause v2.1 publication, flag CinC for adjudication.
- COLLAPSE on both: recommend new §6.5 of v2.1 incorporating Sub F.

## What I will NOT do

- Adjust target set definitions after seeing data.
- Change the seed.
- Run additional non-pre-registered targets.
- Attempt to rescue v2.1's mechanism story by post-hoc rationalisation
  if SUSTAIN.
- Read Mr Adversary's v2.1 review (per Pattern 97 scope protection).

## Files

- [BRIEF.md](BRIEF.md), [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
- `targets.py` — compute Riemann (extension), Gram, RvM-random sets
- `sub_f.py` — main driver: probe + matching + null per target
- `findings_paper_150_sub_f.md` — long-form report
