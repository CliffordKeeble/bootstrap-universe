# 2I-network simulation — findings

**Date**: 25 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 25 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md) (committed at `06cbd9b` before any compute)
**Verdict**: **FAIL — construction (both primary and secondary), with a target-value caveat that needs CinC adjudication.**

---

## Summary

The pre-registered stop-on-fail clause at N = 120 triggered. None of the three
R³ geometric constructions (face-sharing, axis-sharing, vertex-sharing)
produced a λ₁ anywhere near the pre-registered target of 168 ± 17 at N = 120.
The protocol was followed: face-sharing was tried first; it failed; axis-sharing
was tried; it also failed; per protocol I stopped and wrote this report.
N = 600, N = 1200, N = 6000 runs were **not** performed.

**One finding is more important than the construction comparison**: the
600-cell graph itself, constructed exactly as a 12-regular graph on 120
vertices (the brief's known-answer reference target), has

λ₁ = 2.292, λ₁/⟨deg⟩ = 0.191, λ₁·N^(2/3) = 55.76

None of these is within 10% of 168. **The brief's pre-registered target
value cannot be reached at N = 120 by the ideal structure either.** This
is either a normalization discrepancy in the brief, or it implies the
conjecture's specific numerical prediction (168) does not survive contact
with the 600-cell.

A volume-aware reading of the conjecture would predict λ₁·N^(2/3) → 168 · V^(2/3),
where V = vol(S³/2I) = 2π²/120 ≈ 0.165. That gives a target of ≈ 50 at N = 120 —
which the 600-cell matches to within ~12%. This may be the real content of
the conjecture; the literal "168" in the brief appears to be missing a volume
factor.

## Construction chosen and rationale

**Primary**: face-sharing, R³ geometric reflection. New icosahedron = reflection
of parent through the face plane; vertex coincidences detected by spatial-hash
lookup with tolerance ε = 1×10⁻⁶.

This produced a tree (119 edges, avg degree 1.98) up to N=120, because
icosahedra do not face-pack in R³ (dihedral angle ≈ 138.19° doesn't divide 360°).
No consolidation occurred — the "attach to all adjacent nodes" rule from the
brief is geometrically vacuous in R³ since icosahedron-vertex coincidences
under reflection-orbit are measure-zero events.

**Secondary (run after primary failed at N=120)**: axis-sharing, R³ via 180°
rotation about the shared axis. Also a tree (1202 vertices, 119 edges).

**Tertiary (diagnostic per brief allowance)**: vertex-sharing, R³ via
translation by 2·direction along the shared vertex. **This was the only
construction that produced non-trivial consolidation in R³** — 257 edges at
N=120 vs. 119 for a tree. The translation-orbit happens to produce vertex
coincidences in R³ that the rotational/reflective constructions do not.

## Observables table

| construction | N | edges | avg deg | λ₁ (raw) | λ₁/⟨deg⟩ | λ₁·N^(2/3) | d | null λ₁ (mean±std) | z |
|---|---|---|---|---|---|---|---|---|---|
| face       |  12 |  11 | 1.83 | 1.0000 | 0.5455 |  5.241 | n/a  | 0.503 ± 0.376 |  +1.32 |
| face       |  60 |  59 | 1.97 | 0.0477 | 0.0243 |  0.731 | 2.55 | 0.033 ± 0.018 |  +0.80 |
| face       | 120 | 119 | 1.98 | 0.0477 | 0.0241 |  1.161 | 2.66 | 0.005 ± 0.003 | +14.54 |
| axis       |  12 |  11 | 1.83 | 0.2583 | 0.1409 |  1.354 | n/a  | 0.503 ± 0.376 |  −0.65 |
| axis       |  60 |  59 | 1.97 | 0.0464 | 0.0236 |  0.712 | 2.20 | 0.033 ± 0.018 |  +0.73 |
| axis       | 120 | 119 | 1.98 | 0.0266 | 0.0134 |  0.646 | 2.29 | 0.005 ± 0.003 |  +7.37 |
| vertex     |  12 |  11 | 1.83 | 1.0000 | 0.5455 |  5.241 | n/a  | 0.503 ± 0.376 |  +1.32 |
| vertex     |  60 | 109 | 3.63 | 0.7327 | 0.2020 | 11.230 | 1.94 | 0.700 ± 0.090 |  +0.36 |
| vertex     | 120 | 257 | 4.28 | 0.7776 | 0.1816 | 18.919 | 2.18 | 0.616 ± 0.039 |  +4.12 |
| **600-cell ref** | **120** | **720** | **12.00** | **2.2918** | **0.1910** | **55.757** | **1.42** | **5.290 ± 0.058** | **−51.6** |

Raw data in [observables.csv](observables.csv); 600-cell eigenvalues in
[reference_600cell_eigenvalues.csv](reference_600cell_eigenvalues.csv).

### Stop-on-fail thresholds (brief)

- **λ₁ within 10% of 168 (= 151.2–184.8) for at least one normalisation at N=1200**:
  At N=120, no construction or 600-cell reference achieves any value within
  10% of 168. Best is the 600-cell at 55.76 (factor of 3 off).
- **d ∈ [2.8, 3.2] at N=1200**: At N=120, face-sharing gives d=2.66
  (just below the lower bound), 600-cell gives d=1.42 (well below — small
  diameter makes the middle-60% fit window degenerate). Vertex-sharing is
  the closest at 2.18 but still well below.
- **N = 120 wildness check**: triggered. None of the three R³ constructions
  is within an order of magnitude of any reasonable interpretation of 168
  (and even the 600-cell reference is off by a factor of 3).

## Plots (not generated — runs were halted at N=120)

No plots were produced. The brief's "λ₁ vs N" trend would have used N=12,60,120,
600,1200,6000, but per stop-on-fail at N=120 I halted at three points.
Below is the verbal trend of λ₁·N^(2/3) for the three constructions:

```
N    face   axis   vertex   600-cell
12    5.24   1.35    5.24
60    0.73   0.71   11.23
120   1.16   0.65   18.92    55.76
                                             (reference target: 168)
```

The trend at three points is too sparse to extrapolate confidently. Face and
axis are monotonically decreasing or flat; vertex is monotonically increasing
(only construction trending toward the target). The 600-cell, evaluated
exactly, is the upper benchmark for this geometry.

## Heat-kernel return probability

Not reported. The brief lists this as Observable 3 (secondary, informational).
Computation was implemented but only meaningful at large N (≥ 60); given
runs halted at N=120 across three constructions, no asymptotic α fit was
performed. The exponential floor 1/N is too large to identify α reliably
at N=120 with the spectra observed here.

## Heat-kernel-equivalent for 600-cell

For diagnostic value, the 600-cell graph at N=120 has λ₁ = 2.29 with
multiplicity 4 (the four-dimensional 2I-irrep). Density of states near the
gap is governed by this and the next plateau at λ = 5.53 (multiplicity 5),
both being known 2I representation-theoretic values. This rules out an
"accidental" near-zero gap.

## Verdict

**FAIL — construction (primary and secondary)**.

The pre-registered stop-on-fail clause at N=120 fired. All three R³ geometric
constructions miss the 168 ± 17 threshold by 1–3 orders of magnitude across
all three normalisations.

**However**, the 600-cell reference also fails the same threshold at N=120
by a factor of ~3, while matching a *volume-corrected* form of the prediction
(168 · V(S³/2I)^(2/3) ≈ 50) within ~12%. This shifts the responsibility for
re-examination from "is the construction wrong?" to "is the target value
wrong (or is the volume normalisation missing from the brief's prediction)?"

Per the brief's hygiene rules, I did not silently revise the target. I
report literally against the pre-registered 168, and flag the
volume-normalisation issue separately for CinC's adjudication.

## Honest limitations

- **The R³ realisation was always the wrong embedding for closure.** The
  brief explicitly noted "4D embedding or careful combinatorial bookkeeping".
  I committed to R³ in the pre-registration. This made the consolidation rule
  geometrically vacuous for face/axis-sharing — both produced trees. The
  brief anticipated this with the stop-on-fail at N=120.
- **The 600-cell direct construction is the correct known-answer test at
  N=120** and is reported as the reference. A construction whose growth
  *intrinsically produces* the 600-cell at N=120 was not built in this
  investigation; if CinC wants that, it needs a new brief targeting the
  S³-embedded growth rule rather than R³ reflection.
- **The dimension-fit at small N is unreliable.** With diameter ~ 9 (600-cell)
  or ~ 5 (caterpillar trees) at N=120, the "middle 60%" fit window is just
  3–4 R-values. The d-values reported should not be over-interpreted at
  N=120.
- **The null model (random regular graph at matched ⟨deg⟩) is a generic null.**
  It is *not* a "random 2I-coupled network" null. The z-scores show all four
  graphs (face, axis, vertex, 600-cell) are statistically distinguishable
  from the random regular null, but this doesn't tell us whether the
  distinguishability is *because of* 2I structure or just because of the
  particular degree sequence and connectivity.
- **The brief's literal target value 168 was never within reach.** Best
  achievable at N=120, by the 600-cell itself, is 55.76 in the λ·N^(2/3)
  normalisation. A volume-corrected form (≈50) is within reach and is
  approximately matched by the 600-cell.

## Pattern flags

- **Pattern 75 (null)**: satisfied by the random regular graph null at
  matched average degree, 20 samples per N. z-scores reported. The
  2I-structured networks are statistically distinguishable from the random
  null at N=120, but distinguishability alone is not the success criterion —
  the criterion was λ₁ → 168, which no network reaches.
- **Pattern 39 (DERIVED vs OBSERVED)**:
  - DERIVED (symbolic / exact): The 600-cell vertex set was constructed
    exactly from the standard 4D coordinates; its 12-regular adjacency was
    verified analytically (720 = 120·12/2 edges); its Laplacian was computed
    by dense eigendecomposition.
  - OBSERVED (computed / numerical): All graph Laplacian eigenvalues for the
    three R³ constructions; dimension fits; null z-scores. These are
    numerically reliable but interpretation-sensitive.
  - The target value 168 is DERIVED in Paper 117 (not read here per brief
    hygiene) for the continuous Laplacian on S³/2I. The translation to a
    graph-Laplacian prediction at finite N requires a volume factor not
    spelled out in the brief — this is the central flag for CinC.
- **Pattern 19 (adversary)**: An adversary would attack the volume-normalisation
  silence in the brief most aggressively. The fact that even the ideal
  600-cell misses 168 by 3× means any construction that approached 168
  would be a *worse* approximation of the underlying continuum theory, not
  a better one — unless the conjecture's framing is shifted to match the
  volume-corrected target ≈ 50, at which point the 600-cell already
  reasonably matches and the question becomes whether *some* combinatorial
  growth rule produces the 600-cell (or similar) at N=120 from a single-icosahedron
  seed. None of my three R³ rules does this; the brief's anticipated stop-on-fail
  at N=120 is therefore exactly right for the constructions I tried.

## Files

- [BRIEF.md](BRIEF.md) — CinC's brief (verbatim)
- [PRE_REGISTRATION.md](PRE_REGISTRATION.md) — my pre-registered construction choice
- [build_network.py](build_network.py) — three R³ growth implementations + observables + null
- [reference_600cell.py](reference_600cell.py) — explicit 600-cell graph and spectrum
- [run_sweep.py](run_sweep.py) — sweep driver writing observables.csv
- [observables.csv](observables.csv) — all data
- [reference_600cell_eigenvalues.csv](reference_600cell_eigenvalues.csv) — 600-cell Laplacian spectrum (120 values)

## What I did not do (and would in a follow-up brief)

If CinC chooses to re-issue the brief with the volume-normalisation question
resolved, the natural next steps are:

1. **Build a growth rule that produces the 600-cell exactly at N=120.** This
   requires either an S³ (4D) embedding or a combinatorial-closure rule
   that forces sharing at edges (3-icosahedra-per-edge in the icosahedral
   pseudomanifold sense). Both are doable but are new constructions, not
   in the original pre-registration.
2. **Extend the 600-cell to larger N by gluing copies** along icosahedral
   sub-clusters, then test the asymptotic λ₁·N^(2/3) limit. This tests
   whether the conjectured limit is matched, not just at N=120 but in
   the large-N regime where the volume-correction may be the dominant
   issue.
3. **Test heat-kernel α at larger N**, which is a more robust 3-manifold
   diagnostic than λ₁ at fixed small N.
