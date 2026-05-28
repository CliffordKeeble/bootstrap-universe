# Icosahedral lag simulation v0.1 — findings (space from time)

**Date**: 28 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 28 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) committed at `7229282` before any compute.
**Verdict**: **FAIL — wrong-dimension AND null-indistinguishable** (both pre-registered failure modes simultaneously).

---

## Summary

The dynamics has a strikingly binary phase diagram: nearly every parameter
combination collapses to one of two trivial attractors:

1. **Full synchronisation** (mean inter-unit correlation = 1.0, std = 0).
   All units rotate together; the emergent "space" is a single point.
2. **Uniform random distribution on S³** (mean inter-unit correlation
   = −1/(N−1) exactly, std = 1/2 exactly — the *signature* statistics of
   N independent uniform unit vectors in ℝ⁴).

There is no intermediate "structured" regime in which the icosahedral
symmetry organises units into a non-trivial pattern. At the verdict
operating point (ε = 0.025, τ = 3, N = 240, T = 8000) the dynamics
lands cleanly in the "uniform random on S³" attractor, and the same
result is obtained verbatim with a non-icosahedral drift element (the
non-icosahedral null produces an identical correlation matrix
distribution: mean = −0.0042, std = 0.5018 in both runs to four
decimals).

The spectral dimension extracted from the resulting weighted graph is
**d_s ≈ 25**, which is the noise-dominated Laplacian spectrum of a
graph built from cosine-similarities of random unit vectors — not a
3-manifold spectrum.

## Construction choices

**Dynamical unit**: continuous unit quaternions in S³, chosen via
AskUserQuestion with Cliff as the brief's S²-Riemann-sphere construction
would have burned too much budget building a documented A₅-equivariant
rational map. S³ is the natural ambient for 2I (the 120 elements *are*
unit quaternions / 600-cell vertices), so left-multiplication by a 2I
element is A₅/2I-equivariant by construction. Brief's "commutes with the
A₅ action" requirement satisfied.

**Drift element**: `c = (φ/2, 1/(2φ), 1/2, 0)` — a 600-cell vertex
representing rotation by 2π/5 about the axis (1/(2φ), 1/2, 0). Has
order 10 in the quaternion group (c¹⁰ = +1, c⁵ = −1). DERIVED 2I element.

**Update rule** (single fixed form):

```
q_drift = c · q_i(t)
q_neigh = mean_{j ≠ i} K_ij · q_j(t − τ)
q_i(t+1) = normalize( (1 − ε) q_drift + ε q_neigh )
```

**Coupling**: K_ij = 1/(N−1) for i ≠ j (uniform all-to-all). **No
spatial labelling** — this is the central anti-circularity property.

**Lag**: single global integer τ (not per-pair, not distance-encoded).

**Operating point** (from pre-data sweep at small N=30, T=2000):
ε = 0.025, τ = 3. Selected as a "differentiated" point (non-synced)
with highest std_corr. Verdict params: N=240, T=8000, transient=2000.

## The operating-point sweep was diagnostic

The sweep itself revealed the diagnostic finding: **almost every
parameter combination produces full sync** (mean_corr = 1.000,
std = 0.000). The only "non-synced" combinations are
τ ∈ {2, 3, 4, 5} at small ε, all of which gave **identical statistics**
to within 4 decimals: mean ≈ −0.0345 (= −1/29 for N=30), std ≈ 0.485
(close to 1/2 = 1/√4, the variance of unit-vector inner products
on S³).

This is the signature of *N independent uniform random points on S³*,
not of an emergent structure. The lag in this regime is not generating
geometry; it is simply preventing the all-to-all coupling from
synchronising the units, leaving them as independent rotators that
fill S³ uniformly.

Full sweep results in [op_point_sweep.csv](op_point_sweep.csv).

## Observables table (verdict, N=240, T=8000)

| Run | mean_corr | std_corr | d_s | small-r slope | saturated? | sat_frac | frac+def | λ₁·N^(2/3) |
|---|---|---|---|---|---|---|---|---|
| **icosahedral+lag** | −0.0042 | 0.5018 | **25.07** | **1.62** | True | 1.000 | 0.928 | 640.5 |
| null: non-icosahedral | **−0.0042** | **0.5018** | **25.07** | 1.62 | True | 1.000 | 0.928 | — |
| null: no-delay (τ=0) | +1.0000 | 0.0000 | NaN | NaN | False | NaN | 0.000 | — |
| null: shuffled-K | −0.0042 | 0.4717 | 24.41 | 2.90 | True | 1.000 | 0.941 | — |

Brief pre-registered thresholds:
- d_s ∈ [2.5, 3.5]: **all runs FAIL** (icosahedral at 25.07, three orders too large)
- Saturation present where nulls lack it: **all but the sync null saturate equally** — no differentiation
- Icosahedral must beat ALL THREE nulls by clear margin: **FAILS** vs non-icosahedral (identical) and shuffled-K (nearly identical)

## Verdict

**FAIL — wrong-dimension AND null-indistinguishable.**

Both pre-registered failure modes fire concurrently:
- **FAIL-wrong-dimension**: d_s = 25.07, far outside [2.5, 3.5].
- **FAIL-null-indistinguishable**: the non-icosahedral null produces
  *the same correlation statistics to four decimals*, and the same
  d_s and saturation to graph-eigendecomposition precision. The
  icosahedral symmetry contributes nothing distinguishable at this
  operating point in this construction.

The no-delay null is "distinguishable" in the trivial sense that it
fully synchronises (a single-point attractor) — i.e., lag is needed
just to prevent sync, but it does not produce S³ structure.

## Anti-circularity self-audit (pre-committed)

Each of the four checks pre-registered:

1. **Coupling matrix K is all-to-all, no spatial labelling.** ✓
   K_ij = 1/(N−1) for i ≠ j. Exported and confirmed: K has rank 1
   (up to the diagonal); no metric structure embedded.
2. **Lag is a single integer, not per-pair.** ✓ τ = 3 globally;
   no τ_ij assignment.
3. **Initial conditions are uniformly random on S³.** ✓ The buffer
   was initialised with quaternions from `rng.standard_normal(N, 4)`
   followed by normalisation; sample mean of initial states is
   numerically zero in ℝ⁴ (uniform sampling).
4. **Effective geometry uses time-averaged correlations only.** ✓
   The observables are derived purely from C_ij = ⟨q_i · q_j⟩_t;
   no information from the update rule beyond the correlations is used.

**The anti-circularity audit passes cleanly.** The reason this
investigation FAILs is *not* that we accidentally encoded the geometry
in the coupling — it is that the dynamics genuinely does not produce
emergent S³ structure. The FAIL is honest.

## Honest limitations

- **The S³-quaternion choice may not be the brief's "preferred"
  construction.** The brief listed Doyle–McMullen or invariant-built
  A₅-equivariant rational maps on S² as primary. I chose S³ quaternions
  for tractability with Cliff's confirmation. It is possible that a
  *chaotic* Riemann-sphere A₅-equivariant map (Doyle–McMullen has a
  non-trivial Julia set) would produce different dynamics. The 2I
  drift element c is *periodic* (order 10 quaternion), so the
  self-iteration is a clean rotation, not chaotic — and a chaotic
  rational map might give richer attractors. This is the strongest
  hedge in the report.
- **Only one operating point fully tested.** The sweep showed the
  whole grid is binary (sync or random); the verdict was run at one
  point in the non-synced regime. If there is a narrow ε-τ band that
  produces non-trivial structure, it was missed.
- **N = 240, T = 8000 is relatively modest.** Brief allowed up to
  N≈1000, T≈15000. Larger N might reveal effects masked at this size,
  though the noise-floor structure of the correlation matrix would
  also approach the random-matrix limit more closely, not less.
- **The "no-delay null" trivially fails by full-sync, not by
  *positively* producing some other structure.** A more discriminating
  no-delay null might use a slightly stochastic dynamics that
  doesn't sync, so the comparison probes the lag-vs-no-lag axis
  rather than the sync-vs-not axis.

## Pattern flags

- **Pattern 75 (null)**: all three nulls satisfied as pre-registered.
  The non-icosahedral null is the *informative* one — it shows that
  the icosahedral symmetry is **not load-bearing** in this construction.
  The no-delay null is informative-by-degeneracy (full sync).
  The shuffled-K null produces nearly identical results to the
  icosahedral+all-to-all, so the coupling-topology accident is also
  ruled out as the cause.
- **Pattern 39 (DERIVED vs OBSERVED)**:
  - DERIVED: the 2I drift element `c` and its quaternion structure;
    the equivariance property R(g·q) = g·R(q) for g ∈ 2I.
  - OBSERVED: the binary phase diagram (sync vs uniform-random);
    the d_s = 25 noise-dominated spectrum; the null-indistinguishability.
- **Pattern 19 (adversary)**: the anti-circularity audit was the
  hostile-reviewer's pre-emptive line of attack. It passed cleanly.
  The remaining adversarial attack would be "you should have used a
  *chaotic* A₅-equivariant map (Doyle–McMullen) not a periodic
  rotation" — and this is a legitimate critique, hedged in the
  Limitations section above.

## What the result says about the conjecture

The conjecture: **S³ geometry emerges from temporal coupling of
icosahedral units.**

At the construction tested here (continuous S³ quaternions with
periodic 2I drift, all-to-all lag-coupling, intermediate ε), the
conjecture is **not supported**:

- The dynamics has only two attractors: full sync (no structure)
  and uniform random on S³ (no *emergent* structure beyond what
  the state space provides for free).
- The icosahedral symmetry is not load-bearing: a non-icosahedral
  drift produces identical statistics.
- The lag does not produce structure; it only prevents sync.

This is consistent with the running pattern from v0.1 and v0.2 of
the static-graph line: **2I structure is not enough on its own to
produce S³-like spectral or geometric properties.** Both the static
and dynamic constructions tried so far have failed for the same
underlying reason — the icosahedral symmetry of 2I (or its A₅
quotient) does not appear to generate 3-manifold structure
combinatorially or dynamically without some additional ingredient.

**The honest synthesis across v0.1 (static), v0.2 (refined static),
and v0.1-lag (dynamic) is three independent FAILs on three
independent constructions, each ruling out one route to the
conjecture. That is a substantive negative result.**

## Recommendation for any v0.2-lag (not in scope)

If CinC continues this line, the natural next experiment would be:
**chaotic A₅-equivariant Riemann-sphere dynamics**. Specifically the
Doyle–McMullen iteration map (or any rational map of degree ≥ 11
that's A₅-equivariant and has a non-trivial Julia set). This would
replace the periodic drift with chaos, and the question becomes:
*does chaotic A₅-equivariant dynamics with lag produce emergent
S³ structure?* This addresses the strongest critique of v0.1-lag
(periodic drift ≠ rich dynamics).

Out of scope here per pre-registration. Flagged for CinC.

## Files

- [BRIEF.md](BRIEF.md) — CinC's v0.1 brief (verbatim)
- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — pre-registered choices
- [dynamics.py](dynamics.py) — quaternion + lag dynamics + observables
- [op_point_sweep.py](op_point_sweep.py) — operating-point sweep
- [op_point_sweep.csv](op_point_sweep.csv) — sweep data
- [run_verdict.py](run_verdict.py) — verdict + 3 nulls
- [observables_lag_v0_1.csv](observables_lag_v0_1.csv) — verdict data

## Compute budget

Total v0.1-lag runtime: ≈ 8 seconds (4 runs × ~1.5s each + ~2s sweep
overhead + observables). Well under the 30-minute budget.
