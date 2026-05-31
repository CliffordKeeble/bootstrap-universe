# Pre-registration — Sub E (mechanism verification for Paper 150 v2.1)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Binding**: this document is committed to git **before any probe is run and before any minima are identified**.

## Operational definitions (frozen at pre-registration)

**Probe**: `P_d(t) = sqrt(|Re²(Z_φ(t)) − d · Im²(Z_φ(t))|)` where `Z_φ(t)`
is Paper 150 v2.0's Fejér-weighted golden-angle Dirichlet sum
(`compute_zeta_gen` with `alpha = phi`, `norm_type = 'golden'` modified
to use parameterised `d`).

I use the existing `paper-203-sub-c/probe.py:compute_probe(t, d, N)`
which is exactly this with d parameterised.

**Grid**: t ∈ (1, 1000] at Δt = 0.01, so the grid has 99,900 points.
N_terms = 5000 (Paper 150 baseline).

**Local minimum at index i**: P_d[i−1] ≥ P_d[i] ≤ P_d[i+1].

**Non-overlap window** for the per-d top-100 selection: window width
**2** as the brief specifies (i.e., each chosen minimum must be ≥ 2 t-units
from every previously-chosen one).

**Per-d top-100**: I will select the 100 lowest-value local minima
greedily — pick lowest, mark its ±1 t-neighborhood as occupied, pick
next lowest outside occupied regions, repeat until 100.

**Sharing**: for (d, d') with d ≠ d', for each minimum t* of d, count
as "shared" if there exists a minimum t' of d' with |t* − t'| ≤ 2·Δt = 0.02.
**This is a tight test**: minima must align to within 0.02 of each
other to be deemed "shared". Sharing percentage = count / 100.

## Anti-circularity

- The minima are computed from P_d(t) directly, not by comparison with
  any L-zero list. The test is about the probe's internal structure
  across different d values, not its correlation with any specific
  L-function.
- My Flag 1 reading was made in `findings_paper_203_sub_c.md` and
  committed at `4ca0f57` (before the v0.4 brief, before Sub E was
  designed). The test is straightforward verification.

## Pre-committed thresholds (re-stated from brief)

- ≥80% mean off-diagonal sharing → **DECORATION-CONFIRMED**
- 40–80% → **DECORATION-PARTIAL**
- <40% → **DECORATION-REFUTED** (stop and flag)

## What I will NOT do

- Adjust Δt or window width after seeing data.
- Pick a different number of minima than 100 per d.
- Reframe the verdict thresholds.

## Files

- `sub_e.py` — implementation
- `sharing_matrix.csv` — 4×4 results
- `findings_paper_150_sub_e.md` — long-form report
