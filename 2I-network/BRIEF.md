# Mr Code brief — 2I-network simulation (Entry 1)

**Date**: 25 May 2026
**Brief author**: CinC
**Pre-registration document**: commit to git BEFORE running any computation.
**Estimated runtime**: ~1 session for Mr Code. Long enough for overnight if needed.

---

## Task

Build a combinatorial network where each node represents a 2I (binary icosahedral group) instance and edges represent shared 5-fold rotation axes between adjacent 2I's. Compute the graph Laplacian spectrum and the volume-scaling exponent (effective Hausdorff dimension) as the network grows. Test whether the network exhibits the spectral and dimensional properties of S³/2I in the large-N limit.

The conjecture under test: **S³ emerges from networks of coupled 2I instances.** If true, the graph's first non-trivial Laplacian eigenvalue (appropriately normalised) should converge to 168 (the S³/2I spectral gap, Paper 117), and the volume-scaling exponent should approach 3.

## Construction

Each node = one 2I instance, geometrically an icosahedron with 12 vertices defining 6 5-fold rotation axes (each axis through a pair of opposite vertices).

**Primary construction (face-sharing)**: A new node attaches to existing nodes by sharing one *triangular face* (3 vertices). This is the most 3-manifold-natural construction: face-sharing of polytopes produces 3D structure. The icosahedron has 20 triangular faces; max degree per node = 20.

**Secondary construction (axis-sharing)**: A new node attaches by sharing one 5-fold *axis* (2 vertices). Max degree per node = 6.

**Tertiary construction (vertex-sharing, tree-like)**: A new node attaches by sharing one *vertex*. Max degree per node = 12. Provided as diagnostic — likely fails to produce 3-manifold structure; useful as a sanity check on the failure mode.

**Run all three** if time permits. **Mr Code chooses one as primary** and reports the other(s) as comparisons. Recommendation: face-sharing as primary.

**Key reference target**: the 600-cell is the regular 4-polytope with 120 vertices in 2I symmetry that tiles S³. At N=120 with face-sharing, the network should approximate the 600-cell's combinatorial structure. This is a *known-answer test*: if the construction is right, N=120 should produce a graph whose spectral gap matches the 600-cell's known spectrum (which Mr Code can compute or look up). Failure at this point indicates the construction is wrong, not the conjecture.

**Growth procedure**:
1. Start with N=1 (one isolated icosahedron node).
2. At each step, add a node by attaching it to existing node(s) via the shared-face / shared-axis / shared-vertex rule.
3. When multiple existing nodes have free attachment sites adjacent to the new node's intended position, the new node attaches to *all* of them simultaneously (this is what produces consolidation/closure rather than tree-like growth).
4. Continue until target N reached.

Mr Code may need to be careful about the geometric realisation — keeping track of which vertices/faces are "shared" requires either an explicit 4D embedding or a careful combinatorial bookkeeping. Use whichever is cleaner.

## Pre-registered observables and thresholds

**Observable 1: First non-trivial Laplacian eigenvalue λ₁**

Compute the unnormalised graph Laplacian L = D − A. Find the smallest non-zero eigenvalue. Normalisation choices to report:
- Raw λ₁
- λ₁ / (average degree)
- λ₁ · N^(2/3) (the natural scaling for a 3-manifold approximation, since spectral gap on S³ scales like 1/R² ~ N^(-2/3))

Pre-registered prediction: **At least one of these normalisations converges to 168 ± 17 (10% threshold) by N = 1200.**

If none of the three normalisations is within 10% of 168 by N=1200, the test fails. Report the value reached and stop.

**Observable 2: Volume scaling (effective Hausdorff dimension)**

For each node u, compute N_u(R) = number of nodes within graph distance R from u. Average over many starting nodes. Fit:

N(R) ~ a · R^d

for R in a range that avoids the small-R discrete regime and the large-R finite-size cutoff. Report d.

Pre-registered prediction: **d ∈ [2.8, 3.2] by N = 1200.**

If d is outside this range at N=1200, the test fails. Report and stop.

**Observable 3 (secondary): Heat-kernel return probability**

Compute p(t) = (1/N) Σ_u P(walk starting at u returns to u after t steps). For a 3-manifold, p(t) ~ t^(-3/2) at intermediate t.

Pre-registered prediction: **p(t) scales as t^(-α) with α ∈ [1.3, 1.7] in an intermediate range identifiable from the data.**

This is informational; not pass/fail.

## Sizes

Run at: **N = 12, 60, 120, 600, 1200, 6000**.

- N = 12: one icosahedron's worth, sanity check (single node + its 12 immediate neighbours if connectivity reaches that far at small N).
- N = 60: |A₅| scale.
- N = 120: |2I| scale, 600-cell vertex count — known-answer regime.
- N = 600: 600-cell cell count.
- N = 1200: pre-registered decision point.
- N = 6000: if 1200 passes, push to confirm asymptotic behaviour.

Report all six values for both observables.

## Pre-registered null hypothesis (Pattern 75)

The claim under test is structural: the *specific* 2I-coupling rule produces S³-like properties. The null hypothesis is: any random regular graph of the same average degree produces similar properties.

**Null construction**: At each N, generate a random regular graph with the same average degree as the 2I-network and compute the same observables.

Report side-by-side:
- λ₁ (2I-network) vs λ₁ (random regular graph) at each N
- d (2I-network) vs d (random regular graph) at each N

Pre-registered success condition: **the 2I-network's observables are statistically distinguishable from the random null at N ≥ 600**, AND the 2I-network's λ₁ converges toward 168 while the null does not.

If the null produces λ₁ ≈ 168 by accident, the 2I-network's success is uninteresting — it's reporting a generic graph property, not a 2I-specific structural signature.

Reasonable null sample size: 20 random regular graphs per N. Report mean, std, and the 2I-network's value as a z-score against the null distribution.

## Stop-on-fail protocol

If at N = 1200:
- λ₁ is outside 168 ± 17 for all three normalisations, AND
- d is outside [2.8, 3.2]

…then **stop, do not run N=6000, and write the failure report.** A clean failure is as informative as a success.

If at N = 1200:
- λ₁ passes but d fails (or vice versa) — report the partial result and **decide based on the trend**: if the failing observable is converging toward its threshold and N=6000 would plausibly land it, continue; otherwise stop.

If at N = 120 (the 600-cell-equivalent point):
- the spectral gap is wildly off from any reasonable interpretation of 168 — **stop and re-examine the construction**. This suggests the construction rule is wrong, not the conjecture. Try the secondary construction (axis-sharing) before continuing.

## Implementation notes

**Language**: Python. NetworkX for graph construction; scipy.sparse.linalg.eigsh for sparse eigenvalue computation.

**Sparse representation**: at N=6000 the graph Laplacian is 6000×6000 — must be sparse. Average degree should be ≤ 20 so the matrix has at most ~120,000 non-zero entries. scipy.sparse handles this easily.

**Eigenvalue computation**: use eigsh with sigma=0 (shift-invert) to find smallest eigenvalues. Compute the smallest 5-10 eigenvalues, not just the smallest, so the spectral gap can be identified unambiguously (some constructions may have multiple near-zero eigenvalues if the graph is disconnected or near-disconnected at small N).

**Connectivity check**: at every N, verify the graph is connected. If not, the spectral gap calculation needs careful handling (the smallest eigenvalue is 0 for connected; the second-smallest is the spectral gap).

**Dimension fit**: use log-log linear regression of N(R) vs R, with R in the middle 60% of the available range (to avoid small-R and finite-size effects). Report fit and residuals.

**Numerical precision**: 64-bit floats throughout. The eigenvalues will be on the order of 0.1 to 100; no precision issues expected.

**Compute budget**: total runtime should be under 30 minutes for the full N=12 → N=6000 sweep including null. If it exceeds this, profile and report the bottleneck.

## Deliverable

A markdown report appended to this brief, structured as:

```
## Results

### Construction chosen and rationale

Brief notes on which construction (face-sharing / axis-sharing / vertex-sharing) was used as primary and why.

### Observables table

| N | λ₁ raw | λ₁ / <deg> | λ₁ · N^(2/3) | d | null λ₁ | null d | 2I z-score |
|---|---|---|---|---|---|---|---|
| 12 | ... | ... | ... | ... | ... | ... | ... |
| ... | ... | ... | ... | ... | ... | ... | ... |

### Plots (text descriptions; embed as links if generated)

- λ₁ vs N (all three normalisations), with 168 reference line
- d vs N, with 3.0 reference line
- p(t) (heat-kernel return) at largest N reached
- Distinguishability from null: λ₁ z-score vs N

### Verdict

One of:
- **PASS**: λ₁ within 10% of 168 at N=1200 AND d in [2.8, 3.2]. The conjecture survives this test.
- **FAIL — clean**: thresholds clearly missed; specific values reported. The conjecture needs revision.
- **FAIL — construction**: 600-cell test failed; construction rule likely wrong; secondary construction's result reported instead.
- **PARTIAL**: one observable passes, the other fails; trend analysis included.

### Honest limitations

- What was not tested
- What the null might be missing
- Where the construction might be sensitive to choices

### Pattern flags

- Pattern 75 (null): how was it satisfied?
- Pattern 39 (DERIVED vs OBSERVED): which results are derived, which are computed/empirical?
- Pattern 19 (adversary): what would an adversary attack first?
```

## Hygiene

- Pre-registration is binding. Do not adjust thresholds after seeing data.
- If you discover an issue with the construction *before* running observations, stop and flag it for CinC — do not silently revise the brief.
- Do not read Paper 203 draft, the conversation prompting this brief, or Paper 189. The recon should be clean.
- Report negative results with the same emphasis as positive ones. Pattern 19: a clean failure feeds the next iteration; a fudged success poisons it.

⌨️ over to you.

---

## Mr Code's report

**Verdict**: **FAIL — construction (primary and secondary)**. Stop-on-fail
at N=120 triggered. Long form: [findings.md](findings.md).

### Construction chosen and rationale

Pre-registered choice: face-sharing in R³ via reflection through face plane,
with vertex-merge consolidation at ε = 1×10⁻⁶ tolerance. Secondary
(triggered when primary failed at N=120): axis-sharing in R³ via 180° rotation
about the shared axis. Tertiary (diagnostic, run since the other two failed
trivially): vertex-sharing via translation by 2×direction.

Pre-registration committed before any compute at `06cbd9b` (see
[PRE_REGISTRATION.md](PRE_REGISTRATION.md)).

### Observables table

| construction | N | λ₁ raw | λ₁ / ⟨deg⟩ | λ₁ · N^(2/3) | d | null λ₁ | null d | 2I z-score |
|---|---|---|---|---|---|---|---|---|
| face       |  12 | 1.0000 | 0.5455 |  5.241 | n/a  | 0.503 ± 0.376 | 0.84 |  +1.32 |
| face       |  60 | 0.0477 | 0.0243 |  0.731 | 2.55 | 0.033 ± 0.018 | 0.94 |  +0.80 |
| face       | 120 | 0.0477 | 0.0241 |  1.161 | 2.66 | 0.005 ± 0.003 | 0.98 | +14.54 |
| axis       |  12 | 0.2583 | 0.1409 |  1.354 | n/a  | 0.503 ± 0.376 | 0.84 |  −0.65 |
| axis       |  60 | 0.0464 | 0.0236 |  0.712 | 2.20 | 0.033 ± 0.018 | 0.94 |  +0.73 |
| axis       | 120 | 0.0266 | 0.0134 |  0.646 | 2.29 | 0.005 ± 0.003 | 0.98 |  +7.37 |
| vertex     |  12 | 1.0000 | 0.5455 |  5.241 | n/a  | 0.503 ± 0.376 | 0.84 |  +1.32 |
| vertex     |  60 | 0.7327 | 0.2020 | 11.230 | 1.94 | 0.700 ± 0.090 | 1.87 |  +0.36 |
| vertex     | 120 | 0.7776 | 0.1816 | 18.919 | 2.18 | 0.616 ± 0.039 | 2.45 |  +4.12 |
| **600-cell ref** | **120** | **2.2918** | **0.1910** | **55.757** | **1.42** | **5.290 ± 0.058** | **n/a** | **−51.6** |

N = 600, 1200, 6000 **not run**, per stop-on-fail at N = 120.

### Plots

No plots generated (only three N values per construction; runs halted at N=120).
See findings.md "Plots" section for the verbal trend.

### Verdict

**FAIL — construction.** None of the three R³ constructions, nor even the
600-cell reference, gets within 10% of the pre-registered target of 168 in
any of the three normalisations at N=120. The R³ face-sharing rule was
geometrically vacuous (icosahedra don't pack in R³; the consolidation
mechanism never fired), producing a tree. Axis-sharing same. Vertex-sharing
produced the only non-trivial structure but is still 9× off from the target.

**Central caveat for CinC**: the 600-cell graph itself, the brief's known-answer
target, gives λ₁·N^(2/3) = 55.76 — off from the pre-registered 168 by a
factor of 3. A volume-corrected reading (168 · V(S³/2I)^(2/3) ≈ 50) is within
~12% of the 600-cell value. This suggests the literal "168" in the brief
may be missing a volume normalisation factor; this is not a question I can
adjudicate without reading Paper 117 (which the brief tells me not to read).
**Flagging to CinC.**

### Honest limitations

- R³ embedding was the wrong choice for face/axis-sharing closure
  (the brief anticipated this with the N=120 stop-on-fail clause).
- The 600-cell reference is built exactly (DERIVED), but a *growth rule*
  that produces it from a single icosahedron seed was not constructed —
  that requires an S³ embedding or combinatorial-closure rule beyond the
  pre-registration.
- Dimension fits at N=120 are unreliable (diameters ≤ ~5–9, fit window
  too narrow).
- Heat-kernel α (Observable 3) not reported — too coarse at N≤120.

### Pattern flags

- **Pattern 75 (null)**: random regular graph at matched ⟨deg⟩, 20 samples
  per N. z-scores reported. All four graphs distinguishable from null at
  N=120, but distinguishability is not the success criterion.
- **Pattern 39 (DERIVED vs OBSERVED)**: 600-cell construction and spectrum
  are DERIVED; the three growth networks and their observables are OBSERVED.
  The mismatch between the OBSERVED 600-cell value (55.76) and the DERIVED
  brief target (168) is the most consequential pattern-39 flag in this
  investigation.
- **Pattern 19 (adversary)**: the adversary's strongest attack is on the
  brief's pre-registered numerical target itself. The 3× gap between the
  600-cell value and 168 is large enough that no construction could
  reach 168 without violating known graph-Laplacian-vs-continuum scaling
  relations. This is the finding that most warrants CinC's attention.
