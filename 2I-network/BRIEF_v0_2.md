# Mr Code brief — 2I-network simulation v0.2 (convergence test)

**Date**: 26 May 2026
**Brief author**: CinC
**Predecessor**: v0.1 brief, commit 06cbd9b; findings commit 523ae23.
**Pre-registration document**: commit to git BEFORE running any computation.
**Estimated runtime**: ~1 session for Mr Code, well under 1 hour.

---

## What v0.1 found and what v0.2 needs to test

v0.1 ran three R³ growth rules (face-sharing, axis-sharing, vertex-sharing), all of which collapsed to tree-like structures because icosahedra don't tile flat 3-space (dihedral angle ≈ 138.19° doesn't divide 360°). Stop-on-fail at N=120 triggered cleanly.

The key finding was buried in Mr Code's flag, not the headline: **the exact 600-cell graph at N=120 gives λ₁·N^(2/3) = 55.76, which matches the *volume-corrected* target 168·vol(S³/2I)^(2/3) ≈ 50.4 to ~11%.** The brief's literal target (168) was the un-normalised continuum eigenvalue from Paper 117; the right comparison for a finite graph is the volume-normalised version. CinC error in v0.1; corrected here.

**v0.2 mission**: test whether the 600-cell's spectral signature *persists under refinement*. The known-answer test at minimum resolution (N=120) basically passes with proper normalisation; the question is whether λ₁·N^(2/3) stays near 50.4 as N grows from 120 to ~10000 via successive refinement of the 600-cell. If yes, the conjecture has serious evidence beyond a single-point match. If no, the N=120 match was an artefact of being at the special minimum-resolution point.

## Task

Construct refinements of the exact 600-cell graph (the natural minimum-resolution discretisation of S³/2I) at increasing N. Compute the graph Laplacian spectrum and the volume-corrected first eigenvalue at each refinement level. Test convergence: does λ₁·N^(2/3) stay within 10% of 50.4 across refinement levels from N=120 to N≈10000?

## Construction

**Starting point**: the exact 600-cell graph from v0.1's `reference_600cell.py`. 120 vertices, 720 edges, vertex-transitive, 2I-symmetric. Use this verbatim as N=120 and as the seed for refinement.

**Refinement scheme A (primary): iterated edge subdivision.** For each edge (u,v) in the current graph, insert a new vertex w at the edge midpoint, replacing the edge (u,v) with two edges (u,w) and (w,v). All original vertices and edges-by-replacement are preserved; the new vertices sit at edge midpoints.

- Level 0: N = 120, E = 720 (original 600-cell)
- Level 1: N = 120 + 720 = 840, E = 1440
- Level 2: N = 840 + 1440 = 2280, E = 2880
- Level 3: N = 2280 + 2880 = 5160, E = 5760
- Level 4: N = 5160 + 5760 = 10920, E = 11520

**Refinement scheme B (secondary, comparison): cell barycentric subdivision.** For each 3-cell (tetrahedron in the 600-cell), add a vertex at the cell centre and connect it to the cell's 4 vertices. The 600-cell has 600 tetrahedral cells, so this adds 600 vertices.

- Level 0: N = 120 (original)
- Level 1: N = 120 + 600 = 720
- Level 2: requires barycentric subdivision of the level-1 graph — at Mr Code's discretion if level 1 is worth pursuing

**Why two schemes**: if the conjecture holds, BOTH refinement schemes should produce λ₁·N^(2/3) converging to ~50.4 as N grows. If they disagree systematically, the result is construction-dependent and the conjecture has not been cleanly tested. Pattern 9 (multiple paths to truth) applied at the construction level.

**Implementation note**: edge subdivision is a standard NetworkX operation. The 600-cell's edge set comes from `reference_600cell.py` already. Iteration is straightforward.

## Pre-registered observables and thresholds

**Observable 1: λ₁·N^(2/3), the volume-corrected first eigenvalue scaling.**

- Target value: **50.4** (= 168 · vol(S³/2I)^(2/3) = 168 · (2π²/120)^(2/3))
- Pass threshold: **at every refinement level from N=120 to N=10000, |λ₁·N^(2/3) − 50.4| / 50.4 ≤ 10%**
- This is stricter than v0.1's threshold because we're now testing *persistence*, not just minimum-N match.

**Observable 2: convergence rate.**

- Compute the trend of λ₁·N^(2/3) as N grows.
- Pass: the trend is flat-or-converging-toward-50.4. Fail: the trend diverges away from 50.4.
- This is informational; not pass/fail on its own.

**Observable 3: scheme A vs scheme B agreement.**

- At each available common N (or interpolated), compare λ₁·N^(2/3) from scheme A vs scheme B.
- Pass: agreement to within 5% at all comparable N.
- Disagreement >5% means the construction matters and the result is not robustly the 600-cell's spectral signature persisting.

**Observable 4: heat-kernel return scaling.**

- Compute p(t) at the largest N reached in scheme A.
- Pass: p(t) ~ t^(-α) with α ∈ [1.3, 1.7] in an intermediate range.
- This is informational; not pass/fail.

## Pre-registered null hypothesis

For each level in scheme A, generate **20 random graphs** matched on:
1. Same number of vertices N
2. Same number of edges E
3. Comparable degree distribution (use NetworkX's `expected_degree_graph` or equivalent)

Compute λ₁·N^(2/3) for each random graph. Report:
- Mean and standard deviation of the null distribution at each N
- z-score of the 600-cell-refinement value against the null distribution

**Distinguishability requirement**: at every level from N=120 onward, the 600-cell-refinement λ₁·N^(2/3) must be at least 2σ from the null mean. If the null distribution overlaps the refinement value substantially, the spectral signature is generic-graph behaviour, not 2I-specific structure.

## Stop-on-fail protocol

**Stop and report failure if**:
- At any N from 120 to 10000, λ₁·N^(2/3) deviates more than 30% from 50.4 in scheme A. (Tighter than v0.1's threshold; if convergence is real, 30% deviation is well outside what should occur.)
- Or: scheme A and scheme B disagree by more than 15% at any common N.
- Or: at any N from N=840 onward, the z-score against the null is less than 1 (the 2I structure is indistinguishable from a random graph at that scale).

If any of the above triggers at level k, stop. Do not run higher levels. Report the level at which failure occurred and the values.

**Stop and report success if**:
- All four pass-conditions hold cleanly through level 4 (N≈10920).
- z-score against null ≥ 2 at every level.

**Continue with caution if**:
- λ₁·N^(2/3) is between 10% and 30% off but trending toward 50.4 as N grows. Note this as PARTIAL pass and continue.
- Scheme A passes but Scheme B disagrees by 5-15%. Note as construction-sensitive and continue.

## Implementation notes

**Language**: Python, NetworkX, scipy.sparse.linalg.

**Sparse eigenvalue computation**: at N=10920 the Laplacian is 10920×10920 sparse. eigsh with sigma=0 (shift-invert) for smallest 5 eigenvalues. Should complete in seconds per call.

**Random graph generation**: use `networkx.random_regular_graph` or `networkx.expected_degree_graph` depending on which is more faithful to the 600-cell's degree structure. The 600-cell is 12-regular at level 0; later levels have mixed degree distribution. Match the degree distribution, not just the average.

**Reuse v0.1 code where possible**: `reference_600cell.py` provides the seed graph. `build_network.py` may have useful Laplacian/null computation utilities.

**Performance budget**: aim for ≤30 minutes total compute including null. If it runs longer, profile and report.

**Numerical sanity**: at every level, verify connectivity, check that λ₀ is numerically zero (within 1e-10), confirm the first non-zero eigenvalue is unambiguously identifiable.

## Deliverable

A markdown report (`findings_v0_2.md`) committed alongside the v0.2 brief, structured as:

```
## Results

### Refinement scheme used and rationale

Notes on construction choices and any deviations from the brief's primary scheme.

### Observables table (scheme A)

| Level | N | E | avg_deg | λ₁ raw | λ₁·N^(2/3) | dev from 50.4 | null mean | null std | z-score |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 120 | 720 | 12 | ... | 55.76 | +10.6% | ... | ... | ... |
| 1 | 840 | 1440 | ... | ... | ... | ... | ... | ... | ... |
| 2 | 2280 | ... | ... | ... | ... | ... | ... | ... | ... |
| 3 | 5160 | ... | ... | ... | ... | ... | ... | ... | ... |
| 4 | 10920 | ... | ... | ... | ... | ... | ... | ... | ... |

### Observables table (scheme B)

Same columns, levels 0 and 1 minimum.

### Convergence plot description

λ₁·N^(2/3) vs N (log scale), with 50.4 reference line and ±10% / ±30% bands.

### Verdict

One of:
- **PASS**: λ₁·N^(2/3) stays within 10% of 50.4 across all levels; null distinguishability holds; scheme A and scheme B agree. The 600-cell's spectral signature persists under refinement; the 2I × S³ balance conjecture has serious evidence.
- **PARTIAL**: passes some levels but drifts at others, or scheme A/B disagree mildly. Specific values reported.
- **FAIL — diverging**: λ₁·N^(2/3) moves away from 50.4 as N grows. The 600-cell match was a special-case artefact; the conjecture is in trouble.
- **FAIL — null overlap**: the null distribution overlaps the 2I-network's values at higher N. The spectral signature is generic-graph, not 2I-specific.

### Honest limitations

- What's not tested
- Where construction sensitivity might matter
- What the null comparison might be missing

### Pattern flags

- Pattern 75 (null): how was it satisfied at each level?
- Pattern 39 (DERIVED vs OBSERVED): 600-cell graph is DERIVED; refinements are constructed-by-procedure; spectrum is computed-empirically.
- Pattern 19 (adversary): if this passes, what's the first attack a hostile reviewer would mount? (Likely: "your refinement scheme was chosen to give 50.4 by construction." The scheme A vs scheme B agreement check is the answer.)
```

## Hygiene

- Pre-registration is binding. Do not adjust the target (50.4), the pass threshold (10%), or the stop-on-fail thresholds (30%, 15%, z<1) after seeing data.
- If the volume normalisation itself is wrong — i.e., if Paper 117's 168 doesn't reduce to 50.4 by the (2π²/120)^(2/3) factor — flag immediately and stop. CinC will adjudicate before continuing.
- Do not read Paper 203 draft, the conversation prompting this brief, or other programme context. Recon must be clean.
- Report negative results with the same emphasis as positive ones.
- If construction scheme A produces wildly different results from scheme B, treat as evidence that the simulation is construction-sensitive (Pattern 75 angle: the null isn't just random graphs; it's also "would any reasonable refinement scheme give the same answer?").

⌨️ over to you.

---

## Mr Code's report

**Verdict**: **FAIL — diverging** (plus null-distinguishability fail at
level 0 that v0.1 missed). Long form: [findings_v0_2.md](findings_v0_2.md).

### Refinement scheme used and rationale

Both schemes as specified by the brief, no deviations. Scheme A through
level 1 only (stop-on-fail triggered); Scheme B through level 1 (level 2
discretionary per brief and not warranted given Scheme A result).

The 600 tetrahedral cells of the 600-cell were identified as 4-cliques
of the vertex-edge graph (verified count = 600 exactly).

Pre-registration committed before any v0.2 compute at `32dc678`.

### Observables table (Scheme A)

| Level | N | E | avg_deg | λ₁ raw | λ₁·N^(2/3) | dev from 50.4 | null mean | null std | z | protocol |
|-------|---|---|---|---|---|---|---|---|---|---|
| 0 | 120 | 720 | 12.00 | 2.291796 | **55.757** | **+10.55%** | 91.83 | 21.72 | **−1.66** | outside 10%; \|z\|<2 |
| 1 | 840 | 1440 | 3.43 | 0.165660 | **14.748** | **−70.76%** | 9.94 | 2.84 | +1.70 | **STOP-ON-FAIL** |
| 2 | — | — | — | — | — | — | — | — | — | not run |
| 3 | — | — | — | — | — | — | — | — | — | not run |
| 4 | — | — | — | — | — | — | — | — | — | not run |

(Informational: implementation smoke-test gave Scheme A level 2 at
λ₁·N^(2/3) = 5.38, deviation −89.34%. Confirms the diverging trend.
Not part of the official sweep.)

### Observables table (Scheme B)

| Level | N | E | avg_deg | λ₁ raw | λ₁·N^(2/3) | dev from 50.4 | null mean | null std | z |
|-------|---|---|---|---|---|---|---|---|---|
| 0 | 120 | 720 | 12.00 | 2.291796 | **55.757** | **+10.55%** | 91.83 | 21.72 | **−1.66** |
| 1 | 720 | 3120 | 8.67 | 0.809429 | **65.023** | **+28.92%** | 32.21 | 8.57 | +3.83 |

### Convergence

- Scheme A: +10.55% → −70.76% (→ −89% per smoke test). Diverging away,
  monotonic, sharp drop.
- Scheme B: +10.55% → +28.92%. Drifting away in the *opposite*
  direction from Scheme A.
- At level 1, the two schemes differ by a factor of 4.4× (14.75 vs
  65.02). The result is strongly construction-sensitive.

### Verdict (per brief's options)

**FAIL — diverging.** λ₁·N^(2/3) moves away from 50.4 as N grows, in
*both* refinement schemes, in *opposite directions*. The 600-cell
N=120 match was a single-point coincidence, not a refinement-invariant
signature.

Compounding the failure: the **null distinguishability test fails at
level 0**: z = −1.66, |z| < 2. The 600-cell's spectral gap is not
statistically distinguishable from random graphs at matched (N, E,
degree-sequence). The v0.1 "11% match to 50.4" finding therefore
needs to be downgraded — the value sits in the low-gap tail of the
random 12-regular distribution, not at a structurally privileged
position.

### Honest limitations

- **Edge subdivision is not a true 3-manifold mesh refinement.** It
  adds 1D filaments along original edges, not 3D bulk refinement.
  A future Scheme C (proper simplicial cell subdivision into smaller
  tetrahedra) is the natural successor test.
- **Scheme B level 2 not run.** Discretionary per brief; not informative
  given Scheme A failure.
- **`expected_degree_graph` null** preserves degree sequence in expectation
  only; a stricter `configuration_model` null would give a tighter
  distribution but the qualitative |z|<2 finding at level 0 would persist.

### Pattern flags

- **Pattern 75 (null)**: satisfied at every level with 20 samples. The
  level-0 null check is the most consequential finding — the
  600-cell does not stand out at 2σ against random 12-regular graphs.
- **Pattern 39 (DERIVED vs OBSERVED)**: 600-cell graph and target
  50.436 are DERIVED; refined-graph spectra and null comparisons are
  OBSERVED. The v0.1 "11% match" was OBSERVED single-point; v0.2's
  failure of refinement-invariance and null-distinguishability is
  also OBSERVED.
- **Pattern 19 (adversary)**: the adversary's attack vector that v0.2
  was supposed to defuse — "refinement scheme chosen to give 50.4 by
  construction" — instead becomes a *positive* adversary finding:
  *no* reasonable refinement gives 50.4. The v0.1 single-point match
  was construction-dependent, not construction-robust.
- **Pattern 9 (multiple paths to truth)**: Scheme A and Scheme B were
  meant to converge; instead they diverge in *opposite* directions
  by factor 4.4×. Pattern 9 fires as failure detector.

### Compute budget

Total v0.2 runtime: under 10 seconds. Well inside the 30-minute budget.
