# Pre-registration — icosahedral lag simulation v0.1

**Date**: 28 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 28 May 2026)
**Binding**: this document is committed to git **before any computation has been run**.
**Authority**: Cliff explicitly confirmed the two delegated choices below
via AskUserQuestion at session start (S³ quaternions, N=240 / T=8000).

---

## Construction choice (the brief delegates)

**Continuous S³ quaternion states with 2I drift.** Each unit i carries a state
`q_i ∈ S³` (a unit quaternion in ℝ⁴). This is closer to the brief's
preferred Riemann-sphere intent than the discrete-2I fallback, but uses
S³ rather than S² because 2I lives naturally on S³ (the 120 elements of
2I are unit quaternions / vertices of the 600-cell). Continuous state →
the emergent geometry can be richer than the bare 120-point Cayley
structure.

The brief's two listed defaults (Doyle-McMullen on S², or invariants
f₁₂/H₂₀/T₃₀) were both judged too expensive to implement-and-debug
within this session's budget. The S³ construction is **A₅/2I-equivariant
by construction** (the drift element c is in 2I, so `R(q) = c·q`
commutes with the right 2I action), so it satisfies the brief's
"commutes with the A₅ action" requirement.

## Dynamics

Single fixed-form update rule applied at every step `t` to every unit `i`:

```
q_drift  = c · q_i(t)                                 # self-iteration, fixed 5-fold
q_neigh  = mean_{j ≠ i} K_ij · q_j(t − τ)             # lagged neighbour mean in R^4
q_i(t+1) = normalize( (1 − ε) · q_drift  +  ε · q_neigh )
```

Where:
- `c ∈ 2I` is a fixed 5-fold rotation quaternion. Specifically
  `c = (cos(π/5), sin(π/5)·n̂)` with `n̂` aligned to a 600-cell vertex
  axis. Concretely I'll use the quaternion `(φ/2, 1/(2φ), 1/2, 0)`,
  which is in 2I (it's a 600-cell vertex) and represents a 2π/5
  rotation about the axis (1/(2φ), 1/2, 0)/||·||.
- `ε ∈ (0, 1)` is the coupling strength. Operating point fixed via
  the pre-data sweep described below.
- `τ ≥ 0` is the global integer lag in time-steps. Operating point
  fixed via sweep.
- `K_ij` is **all-to-all uniform**: `K_ij = 1/(N-1)` for `i ≠ j`,
  `K_ii = 0`. (For the shuffled-coupling null, K becomes a sparse
  random matrix; see nulls below.)
- The `mean_{j ≠ i} K_ij · q_j(t − τ)` is computed in R⁴ (treating
  quaternions as 4-vectors), so it's a uniform weighted average.
- `normalize` projects back to the unit sphere S³: `q / ||q||`.
- Quaternion multiplication `c · q` is standard Hamilton product.

## Sign convention for antipodes

Quaternions q and -q represent the same SO(3) rotation. For the
emergent-geometry test we treat them as **distinct points on S³**
(don't fold). Rationale: the brief tests for S³ structure, not
SO(3). If folding were intended, the brief would say S³/{±1} = SO(3).

For dimension and curvature on the symmetric-space SO(3), the right
fold could be revisited if results are ambiguous; pre-reg uses the
S³ convention.

## Effective geometry reconstruction

1. **Correlation**: `C_ij = ⟨ q_i(t) · q_j(t) ⟩_t` averaged over the
   steady-state recording window. Quaternion dot product ∈ [-1, 1]
   (cos of half-angle between rotations).
2. **Distance**: `d_ij = arccos(C_ij)` ∈ [0, π]. Brief offered
   `arccos` or `-log`; I pick `arccos` because it has the natural
   S³ geodesic-distance interpretation under the spherical embedding.
3. **Weighted graph**: edge weight `W_ij = max(0, C_ij)` (drop
   anti-correlated pairs for the random-walk graph). Sanity: at
   the operating point I'll confirm the graph is connected when
   restricted to positive correlations; otherwise I'll use
   `W_ij = exp(-d_ij²/σ²)` with σ = median(d).

## Observable computations

- **Observable 1 (spectral dimension d_s)**: graph Laplacian
  `L = D − W`. Compute the smallest `k=20` eigenvalues with eigsh
  (sigma=0); approximate the heat kernel return
  `p(t) ≈ (1/N) Σ_k exp(-t·λ_k)`; fit `log p(t) ~ -(d_s/2)·log(t) +
  const` over the intermediate range (middle 60% of useful log-t,
  same convention as v0.1).
- **Observable 2 (saturation)**: for each unit u, compute
  `N_u(r) = #{j : d_ij ≤ r}` averaged over all u. Plot `N(r)` vs `r`.
  Check (a) small-r power-law slope (3 expected if d_s = 3), (b)
  saturation at `N` (compact).
- **Observable 3 (curvature)**: triangle comparison.
  For each triple (i, j, k) sampled from the largest connected
  component, compute the triangle-defect
  `δ = d_ij + d_jk + d_ki − π`. On S³ of geodesic radius ≤ π,
  a "spherical" triangle has positive defect.
  Report sign of mean δ over a random sample of triples.
- **Observable 4 (only if O1-O3 pass)**: graph Laplacian λ₁,
  normalised as `λ₁·N^(2/3)`, compared to v0.2 target 50.4.

## Operating-point sweep

Pre-data sweep at small N=30, T=2000 over a small grid:
- ε ∈ {0.05, 0.1, 0.2, 0.4}
- τ ∈ {0, 1, 5, 20}

Pick the (ε, τ) pair where:
(a) units do NOT trivially synchronise (mean inter-unit correlation < 0.95)
(b) units do NOT trivially decorrelate (mean inter-unit correlation > 0.05)
(c) variance of the inter-unit correlation distribution is highest
    (most differentiated structure)
(d) the dynamics is stable in the second half of the run (no
    divergence; first-half-vs-second-half mean correlation difference
    < 0.05)

If no operating point satisfies all four, expand the grid up to
ε ∈ {0.025, 0.05, …, 0.6}, τ ∈ {0, 1, 2, 5, 10, 20, 50}.

If the dynamics still does not settle, report **FAIL–no-attractor**
and stop. This is one of the brief's named failure modes.

## Verdict run

After operating point fixed:
- N = 240 (Cliff confirmed)
- T = 8000 total steps; first 2000 discarded as transient; correlations
  computed over the last 6000.
- Same (ε, τ) as the operating point.
- Seed 20260528 for all randomness.

## Three nulls

Run the identical pipeline (same N, T, ε, τ, seed) with one change each:

1. **Non-icosahedral null**: replace the 2I drift element `c` with a
   generic non-icosahedral rotation: `c_null = (cos(π/3.7), sin(π/3.7)·n̂)`
   for an irrational angle and a fixed axis `n̂ = (1, 0, 0, 0)/||·||`.
   No 2I commutation property.
2. **No-delay null**: `τ = 0`. Same icosahedral `c`. Same ε.
3. **Shuffled-coupling null**: same icosahedral `c`, same τ, but
   `K_ij` is sparse random with the same mean and same total weight
   as the all-to-all case. Specifically: pick the `N·(N-1)/2`
   pairs at random with probability `p = 1.0` (i.e., the all-to-all
   is full anyway), then shuffle the weights to random unit
   quaternions on the matrix. The brief notes this null is
   degenerate for fully-uniform all-to-all; I commit to
   `K_ij ~ Uniform[0, 2]` (mean 1) instead of constant, then shuffle
   per-pair.

Distinguishability requirement (pre-registered): at the verdict N,
the icosahedral+lag d_s must be ≥ 0.5 closer to 3.0 than any null's
d_s; AND saturation must be present in icosahedral+lag and absent in
at least one null.

## Anti-circularity self-audit (pre-registered)

I list now the checks I will run AFTER the verdict to confirm no
hidden circularity:

1. The coupling matrix `K` is N×N, all-to-all, *no spatial labelling*.
   I will export the K and confirm no metric structure was encoded.
2. The lag τ is a single integer, *not per-pair*. I will confirm
   `τ_ij` was not assigned.
3. The state space S³ has a natural metric, but the *initial conditions*
   are randomly sampled uniformly on S³ — no clusters pre-placed.
   I will verify init distribution is uniform.
4. The emergent geometry uses time-averaged correlations only — no
   dependence on the dynamics' update rule beyond the correlations.

If any of these checks reveal hidden geometry, the PASS is downgraded
to a circular-construction FAIL regardless of observables.

## What I will NOT do

- Re-tune ε or τ after seeing verdict observables.
- Switch dynamics or coupling after seeing the verdict.
- Read Paper 203, the brief-generating conversation, or 2I-network
  findings during this work.
- Fold antipodes (q ↔ -q) after the fact to rescue a result.

## Files

- `dynamics.py` — quaternion + lag implementation
- `run_v0_1.py` — driver: sweep → verdict → nulls
- `observables_lag_v0_1.csv` — data
- `findings_icosahedral_lag_v0_1.md` — long-form report
