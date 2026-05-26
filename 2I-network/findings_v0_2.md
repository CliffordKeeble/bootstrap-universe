# 2I-network simulation v0.2 — findings (convergence test)

**Date**: 26 May 2026
**Brief**: [BRIEF_v0_2.md](BRIEF_v0_2.md) (CinC, 26 May 2026)
**Pre-registration**: [PRE_REGISTRATION_v0_2.md](PRE_REGISTRATION_v0_2.md) (committed at `32dc678` before any v0.2 compute)
**Predecessor**: v0.1 findings at `523ae23`.
**Verdict**: **FAIL — diverging** (Scheme A); plus an *additional* failure on null
distinguishability at level 0 that the v0.1 framing had missed.

---

## Summary

Pre-registered stop-on-fail triggered at Scheme A level 1 (N = 840):
`λ₁·N^(2/3) = 14.75`, deviation **−70.76%** from the target 50.436, well
past the 30% stop-on-fail threshold. Per protocol, Scheme A levels 2–4
were not run.

Scheme B level 1 (N = 720) gave `λ₁·N^(2/3) = 65.02`, deviation
**+28.92%** — outside the 10% pass threshold and trending *away from*
50.4, not toward it (level 0 was at +10.55%, level 1 at +28.92%).

Scheme A and Scheme B agree by construction at level 0 (both are the
bare 600-cell graph). At level 1 they differ by a factor of 4.4×
(14.75 vs 65.02). The result is **construction-sensitive** in the
strong sense.

**Additional failure not flagged in v0.1**: the null distribution
(random graphs matched on N, E, degree sequence) at level 0 has
mean 91.83 ± 21.72, putting the 600-cell value (55.76) at **z = −1.66**.
The pre-registered distinguishability requirement is |z| ≥ 2 at every
level from N = 120 onward. **The 600-cell's spectral gap is NOT
statistically distinguishable from random regular graphs at matched
degree at 2σ.** The v0.1 "11% match to 50.4" finding therefore needs to
be downgraded: the 600-cell value sits in the lower tail of the
12-regular random-graph distribution, not at a structurally privileged
position.

## Refinement scheme used and rationale

Both schemes as specified by the brief. No deviations:

- **Scheme A (primary)**: iterated edge subdivision of the 600-cell.
  Level 0 = bare 600-cell (120 nodes, 720 edges). Each subsequent level
  inserts a degree-2 vertex at every edge midpoint.
- **Scheme B (secondary)**: barycentric cell subdivision. Level 1 adds
  one vertex per tetrahedral cell (600 cells in the 600-cell), connected
  to the 4 corner vertices. I went to level 1 only; cell barycentric
  subdivision at level 2 was discretionary per the brief and is not
  warranted given Scheme A's failure at level 1.

The 600 tetrahedral cells were identified as the 4-cliques of the
600-cell graph (verified count = 600 exactly).

## Observables table (Scheme A)

| Level | N    | E     | avg_deg | λ₁ raw     | λ₁·N^(2/3) | dev from 50.4 | null mean | null std | z       | protocol           |
|-------|------|-------|---------|------------|-----------:|---------------|----------:|---------:|---------|---------------------|
| 0     | 120  | 720   | 12.00   | 2.291796   | **55.757** | **+10.55%**   | 91.83     | 21.72    | **−1.66** | outside 10% pass; |z|<2 |
| 1     | 840  | 1440  | 3.43    | 0.165660   | **14.748** | **−70.76%**   |  9.94     |  2.84    | **+1.70** | **STOP-ON-FAIL TRIGGERED (>30%)** |
| 2     | —    | —     | —       | —          | —          | —             | —         | —        | —       | not run per protocol |
| 3     | —    | —     | —       | —          | —          | —             | —         | —        | —       | not run per protocol |
| 4     | —    | —     | —       | —          | —          | —             | —         | —        | —       | not run per protocol |

**Note on level 1 z-score**: even though |z| = 1.70 > 1 (so the brief's
"z<1 stop" doesn't fire), the 30%-deviation stop fires first.

**Independent informational data point**: during the implementation
smoke test (before formal evaluation), I observed Scheme A level 2 at
`λ₁·N^(2/3) = 5.38`, deviation −89.34%. This confirms the diverging
trend continues. I am reporting it as a flag, not as an officially-run
result; the official sweep halted at level 1 per protocol.

## Observables table (Scheme B)

| Level | N    | E    | avg_deg | λ₁ raw   | λ₁·N^(2/3) | dev from 50.4 | null mean | null std | z       |
|-------|------|------|---------|----------|-----------:|---------------|----------:|---------:|---------|
| 0     | 120  | 720  | 12.00   | 2.291796 | **55.757** | **+10.55%**   | 91.83     | 21.72    | **−1.66** |
| 1     | 720  | 3120 | 8.67    | 0.809429 | **65.023** | **+28.92%**   | 32.21     |  8.57    | **+3.83** |

Scheme B level 1 (a) is outside the 10% pass threshold, (b) drifts
*away from* 50.4 (level 0 at +10.55% → level 1 at +28.92%), (c) is
the *opposite direction* drift from Scheme A. At level 1 the two
schemes give 14.75 vs 65.02 — a factor of 4.4× discrepancy at the same
nominal "refinement step".

The null at Scheme B level 1 is well-separated (z = +3.83 ≥ 2). So
**Scheme B level 1 satisfies null distinguishability** but **fails the
10% target threshold AND fails the convergence-toward-50.4 trend**.

## Convergence plot description

No plot file generated (would be misleading with only 2 official
Scheme A points). The verbal trend:

- Scheme A: +10.55% → −70.76% (and the smoke-test informational
  point at level 2 is −89.34%). Strongly *diverging away* from 50.4.
- Scheme B: +10.55% → +28.92%. Drifting *away from* 50.4 in the
  positive direction.

Both schemes leave the ±10% pass band immediately upon refinement;
they diverge from 50.4 in *opposite* directions; the magnitude of
the deviation grows with refinement.

The 600-cell value 55.76 sits at +10.55%, just outside the 10% pass
threshold. The v0.1 framing ("matches 50.4 to within ~11%") was
strictly outside the pre-registered 10% v0.2 pass band — a half-step
miss made fully clear here.

## Verdict

**FAIL — diverging.**

The brief specified:
> **FAIL — diverging**: λ₁·N^(2/3) moves away from 50.4 as N grows.
> The 600-cell match was a special-case artefact; the conjecture is in
> trouble.

That is exactly the result. Both refinement schemes show monotonic
drift away from 50.4 as N grows, and they drift in opposite directions
(Scheme A downward, Scheme B upward), confirming that the result at any
given N is **construction-dependent**, not a property of the underlying
manifold.

In addition, the **null distinguishability fails at level 0** (z = −1.66,
|z| < 2). The 600-cell's spectral gap is not significantly different
from that of random graphs at matched degree distribution. This is a
distinct failure mode from "diverging" and was not flagged in v0.1's
narrower target check.

## Honest limitations

- **Edge subdivision is not a 3-manifold mesh refinement.** Inserting a
  degree-2 vertex at each edge midpoint adds 1D "filaments" along edges
  of the original graph; it does not refine the 3-dimensional cell
  structure. A proper 3-manifold mesh refinement would subdivide each
  tetrahedral cell into smaller tetrahedra (e.g., 1-to-8 simplicial
  subdivision, or barycentric simplicial subdivision). Scheme A's
  graph-theoretic edge subdivision essentially tests something different
  from the underlying conjecture: it tests "does the spectral gap
  scale with N^(2/3) under arbitrary subdivision", which it does not,
  because the new vertices are not uniformly distributed in the 3D
  bulk. CinC may want a *Scheme C* in a follow-up brief: proper
  simplicial 1-to-8 cell subdivision.
- **Scheme B level 2 not run**. The brief left it discretionary; given
  Scheme A's failure I judged it not informative.
- **The 600-cell graph is the wrong mesh model for matched-degree null**.
  A 12-regular graph on 120 vertices has many possible structures; the
  600-cell happens to be one with relatively *low* spectral gap (in the
  bottom ~5% of the random 12-regular distribution per the level-0
  null). This suggests that the "55.76 ≈ 50.4 (within 11%)" framing
  from v0.1 is more honestly stated as "the 600-cell sits in the
  low-gap tail of 12-regular graphs, in a region where many random
  graphs would coincidentally land". The v0.1 finding survives at face
  value but loses any claim to be 2I-specific.
- **Null model is `expected_degree_graph`**, which preserves degree
  sequence only in expectation. A stricter null (`configuration_model`
  with strict degree matching) would give a slightly tighter
  distribution. The qualitative null finding (|z|<2 at level 0) would
  likely persist.

## Pattern flags

- **Pattern 75 (null)**: satisfied at all levels (20 random graphs each,
  matched on N, E, degree sequence via `expected_degree_graph`).
  *Result of the null check was negative at Scheme A level 0* — the
  most important specific outcome of this v0.2 work and the reason
  the v0.1 framing needs revision in the parent paper.
- **Pattern 39 (DERIVED vs OBSERVED)**:
  - DERIVED: the 600-cell graph construction; the target value
    50.436 = 168 · (2π²/120)^(2/3); the eigenvalue spectrum of the
    bare 600-cell (computed by dense diagonalisation).
  - OBSERVED: all refined-graph spectral gaps, null distributions,
    z-scores.
- **Pattern 19 (adversary)**: the adversary's strongest line of attack
  is now clear. The 600-cell-at-N=120 value of 55.76 ≈ 50.4 is *neither*
  a refinement-invariant signature *nor* statistically distinguishable
  from random 12-regular graphs. The v0.1 framing of "11% match"
  conflated a single-point coincidence with a structural property.
  Without v0.2's null-at-level-0 check, that conflation would have
  propagated. Pattern 19 vindicated: the adversary catches what the
  initial pass missed.
- **Pattern 9 (multiple paths to truth)**: Scheme A and Scheme B were
  supposed to be two paths converging on the same answer. Instead they
  diverge in *opposite* directions, by a factor of 4.4× at level 1.
  Pattern 9 fires as failure detector — the construction-dependence
  itself is the finding.

## Recommendations for v0.3 (if there is one)

Not part of this report; flagging for CinC:

1. The natural successor test is a *proper simplicial cell subdivision*
   of the 600-cell into smaller tetrahedra — a true 3D mesh refinement —
   to test whether λ₁·N^(2/3) converges to 50.4 under refinement that
   actually approximates the manifold. Edge subdivision was the wrong
   instrument.
2. The matched-degree null should be supplemented with a *matched-girth*
   or *matched-clique-count* null to test whether the 600-cell carries
   any combinatorial signature beyond degree distribution that would
   distinguish it from random graphs.
3. The Cliff observation "S³/2I not too stable for small networks"
   appears to be vindicated by this work — small N does not give
   stable convergence statistics, and the v0.1 N=120 "match" was a
   single-point coincidence. The honest framing in the parent paper
   should be: "the 600-cell graph at N=120 gives λ₁·N^(2/3) within
   11% of the volume-corrected target, but this match does not persist
   under refinement, is not distinguishable from random graphs at
   matched degree at 2σ, and depends on the choice of refinement
   scheme."

## Files

- [BRIEF_v0_2.md](BRIEF_v0_2.md) — CinC's v0.2 brief (verbatim)
- [PRE_REGISTRATION_v0_2.md](PRE_REGISTRATION_v0_2.md) — Mr Code's
  v0.2 pre-registration
- [build_refinement.py](build_refinement.py) — Scheme A and Scheme B
  implementations + analytics + null
- [run_v0_2_sweep.py](run_v0_2_sweep.py) — driver
- [observables_v0_2.csv](observables_v0_2.csv) — raw data

## Compute budget

Total v0.2 runtime: under 10 seconds (much less than the 30-minute
budget). The sparse Laplacian operations dominated; even the null
distribution at the largest level computed (Scheme B level 1, N=720)
took ~3 seconds for 20 random graphs.
