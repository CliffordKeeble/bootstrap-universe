# Pre-registration v0.2 — 2I-network convergence test

**Date**: 26 May 2026
**Author**: Mr Code
**Brief**: [BRIEF_v0_2.md](BRIEF_v0_2.md) (CinC, 26 May 2026)
**Binding**: this document is committed to git **before any v0.2 computation has been run**.

## Volume-correction sanity check (not data-dependent)

The brief's target value 50.4 is `168 · vol(S³/2I)^(2/3)`. I verified this
arithmetically before pre-registering:

- `vol(S³)` (unit 3-sphere in R⁴) = `2π²` ≈ 19.7392
- `vol(S³/2I)` = `2π²/120` ≈ 0.164493
- `(2π²/120)^(2/3)` ≈ 0.29978
- `168 × 0.29978` ≈ **50.36**

This matches CinC's stated 50.4 to better than 0.1%. The volume-correction
formula in the brief is consistent with the spectral-geometry paper read
alongside the brief (which gives λ_1(S³/2I) = 168 in the unit-radius S³
convention with the Ikeda 1980 calculation as the source). I have **no
flag** on the target value itself, so I proceed.

This brief reads paper 117/spectral-geometry and paper 119 (Fibonacci ladder)
under explicit extension of trust by Cliff alongside the brief. The
forbidden papers (203 draft, conversation prompting brief) remain unread.

## Refinement-scheme commitments

I commit to the brief's two refinement schemes exactly as specified:

**Scheme A (primary, edge subdivision)**: for each edge `(u, v)` insert
midpoint `w` and replace `(u, v)` with `(u, w)` and `(w, v)`. Iterate
through levels 0..4. The level table from the brief:

| Level | N     | E     |
|-------|-------|-------|
| 0     | 120   | 720   |
| 1     | 840   | 1440  |
| 2     | 2280  | 2880  |
| 3     | 5160  | 5760  |
| 4     | 10920 | 11520 |

**Scheme B (secondary, barycentric cell subdivision)**: for each
tetrahedral cell of the 600-cell, add a vertex at the cell centre and
connect it to the 4 cell vertices. Level 1 gives N = 720, E = 720 +
4·600 = 720 + 2400 = 3120 (original 720 edges + 2400 new cell-to-vertex
edges). I will go to level 1 only unless time/budget permits more.

## Null model

20 random graphs per level matched on:
1. **Exact N** (vertex count)
2. **Exact E** (edge count)
3. **Degree sequence** — to match the heterogeneous degree distribution
   of refined networks (level ≥ 1 has degree-12 originals + degree-2
   midpoints in Scheme A; degree-32 originals + degree-4 cell-centres
   in Scheme B level 1)

I commit to `nx.expected_degree_graph(degree_sequence, selfloops=False)`
as the primary generator, **fallback** to `nx.configuration_model` with
self-loop / multi-edge removal if `expected_degree_graph` produces too
many disconnected components.

The brief allows either; I commit to expected_degree_graph as primary
because it's simpler and more standard.

Random seed: `20260526`.

## Eigenvalue computation

- `scipy.sparse.linalg.eigsh` with `sigma=0` (shift-invert) for the
  smallest 5 eigenvalues at each level.
- Dense fall-back `numpy.linalg.eigvalsh` for N ≤ 1500 as a sanity
  cross-check; otherwise sparse only.
- Spectral gap = smallest eigenvalue > 1e-9 (numerical zero threshold).
- Connectivity check at every level: if disconnected, work on the
  largest component and flag.

## Heat-kernel α (Observable 4)

Informational, not pass/fail per brief. I will compute it only at the
largest N reached in Scheme A, using the smallest ~200 eigenvalues
sparsely if N > 1500. Fit α over the intermediate t range where p(t)
- 1/N > 1e-12. Fitting window: middle 60% of the log(t) range, same
as the v0.1 dimension fit.

## What I will NOT do

- I will not adjust the target value 50.4 after seeing data.
- I will not adjust the 10%, 15%, 30%, or z-thresholds.
- I will not switch refinement schemes if Scheme A fails — Scheme A
  is the primary; if it fails I report and stop per the brief's
  stop-on-fail protocol.
- I will not silently revise the brief; if I find an issue I flag and
  stop.

## Files

- `build_refinement.py` — Scheme A and Scheme B implementations.
- `run_v0_2_sweep.py` — driver, observables, null, CSV output.
- `observables_v0_2.csv` — all data.
- `findings_v0_2.md` — long-form report.
- Existing `reference_600cell.py` is reused verbatim as the level-0 graph.
