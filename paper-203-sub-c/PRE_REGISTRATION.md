# Pre-registration — Paper 203 Sub C (discriminant variation)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Binding**: this document is committed to git **before any zero data is pulled and before any computation has been run**.
**Pattern 97 note**: Papers 150 v2.0, 196 §5, 164 §2 are NOT in the local
repo. I am working from CinC's structural summary in the brief, and from
the Paper 150 v2.0 companion code at `golden-zeta/` (which encodes the
"what worked and what was pre-registered" baseline directly). If CinC
finds my interpretation of Paper 196 §5's prediction inconsistent with
the paper text, the diagonal-vs-matrix logic should be re-clarified
before findings are written.

---

## Construction commitments

### Matrix dimensions and discriminants

I commit to the brief's specified 4×4 matrix:
- **Norm rows (probe)**: d ∈ {3, 5, 7, 13} (squarefree discriminants).
- **Character columns (target)**: same set d' ∈ {3, 5, 7, 13}.
- **Conductors**: d ≡ 1 (mod 4) → conductor = d; d ≡ 3 (mod 4) → conductor = 4d.
  - (d = 5, cond = 5), (d = 13, cond = 13), (d = 3, cond = 12), (d = 7, cond = 28).

### Character tables (computed by hand from Kronecker symbol; reproduced in code)

For each character χ_D, the values on residues mod D are:

**χ_5 (D=5)**, primitive real character mod 5:
- residue 0: 0
- 1: +1, 2: −1, 3: −1, 4: +1

**χ_13 (D=13)**, Legendre symbol mod 13:
- residue 0: 0
- 1: +1, 2: −1, 3: +1, 4: +1, 5: −1, 6: −1,
- 7: −1, 8: −1, 9: +1, 10: +1, 11: −1, 12: +1

**χ_12 (D=12)** for ℚ(√3):
- gcd(n, 12) > 1: 0
- 1: +1, 5: −1, 7: −1, 11: +1

**χ_28 (D=28)** for ℚ(√7):
- gcd(n, 28) > 1: 0
- 1: +1, 3: +1, 5: −1, 9: +1, 11: −1, 13: −1,
- 15: −1, 17: −1, 19: +1, 23: −1, 25: +1, 27: +1

Cross-checks in code:
- Sum of χ over a full period = 0 (primitive real characters are balanced).
- Multiplicativity: χ(ab) = χ(a)χ(b).
- These will be verified before zero finding.

### L-function zero source

The brief allows "LMFDB or in-script computation". I have no network
access for LMFDB pulls. I commit to **in-script computation using
mpmath**:

```
L(s, χ_q) = q^(−s) · Σ_{a=1}^{q−1} χ(a) · ζ(s, a/q)
```

where ζ(s, a/q) is the Hurwitz zeta (`mpmath.zeta(s, a/q)`). This is
exact analytic continuation; no truncation.

**mpmath precision (dps)**: 15 (default) for initial zero finding;
20 for the final refined zero list. The matching window is W = 1.0,
so zero precision better than 10⁻³ is more than sufficient.

**Zero finding procedure**:
1. Compute |L(1/2 + it, χ)| on a fine grid t ∈ [0.5, T_max] with dt = 0.1.
2. Identify local minima below a fraction (10%) of the median magnitude.
3. For each minimum candidate, refine with `mpmath.findroot` against
   the real and imaginary parts of L.
4. Save (n, t_n) to per-character CSV files.

**Coverage**: T_max chosen so that the number of L-zeros in [0, T_max]
is ~ 100 (matching the Paper 150 v2.0 N=100 reference). Specifically
T_max set per character so that ~200 zeros are found. For χ_5, that's
~ T_max ≈ 250 (matching Paper 150's t range). For higher conductors,
zeros are denser, so T_max smaller (~ 200 for χ_13).

For the N = 10⁴ extension on diagonal cells, T_max correspondingly
larger so that ~ 10⁴ zeros are within range.

### Probe construction

Per Paper 150 v2.0 `zeta_gen.py`, with the brief's discriminant
parameterisation:

```python
re = Σ_{n=1}^N σ_n · n^(−1/2) · cos(t·log(n) + θ_n)
im = Σ_{n=1}^N σ_n · n^(−1/2) · sin(t·log(n) + θ_n)
probe_d(t) = sqrt(|Re²(t) − d · Im²(t)|)
```

with `θ_n = 2π · {n · φ}` (golden-angle phase, Paper 150 baseline)
and `σ_n = 1 − n/N` (Fejér smoothing).

### Minima detection and matching

Per Paper 150 v2.0 `zeta_phi.py`:
- 5th percentile threshold (the probe's lowest 5% of values).
- Three-point local-minima detection.
- Deduplication window 0.3 (probe-local consolidation).
- Matching window W = 1.0 (the bound for "matched" classification).
- Mean Δ = average of distances from each L-zero to the nearest
  probe minimum.

### Null

Per Paper 150 v2.0:
- N_null = 1000 Monte-Carlo trials.
- Each trial: draw fake "minima" uniformly on [t_min, t_max] at the
  same density as the actual probe minima.
- Compute mean Δ for each trial.
- Empirical z = (signal_mean − null_mean) / null_std.

**Fixed random seed**: 42 (matching Paper 150 v2.0 baseline).

### Pipeline-soundness reproduction

The (d, d') = (5, 5) cell is the Paper 150 v2.0 baseline. I commit to:
1. Computing L(s, χ_5) zeros and comparing to the Riemann zeros at the
   matched indices to confirm the L_5 zero list is well-formed.
2. Running the (5, 5) cell at N = 1000 and confirming z ≈ −8.14 (Paper
   150 v2.0 baseline).
3. If the (5, 5) cell does NOT reproduce, STOP per the brief's
   pipeline-stop clause.

**Important reproduction note**: Paper 150 v2.0 baseline targets the
**Riemann zeta zeros** (not L(χ_5) zeros). The brief's (5, 5) cell
generalises to L(χ_5) zeros; these differ from Riemann zeros. I will
run BOTH:
- (5, 5)_Riemann: probe |Re² − 5·Im²| matched against Riemann zeros.
  This is the Paper 150 v2.0 reproduction (target z ≈ −8.14).
- (5, 5)_chi: probe |Re² − 5·Im²| matched against L(χ_5) zeros.
  This is Sub C's diagonal cell at (5, 5).
The first is the pipeline-soundness test. The second is the matrix
diagonal.

### Effect size

For Sub C-2:
```
ε(N) = (mean Δ_null − mean Δ_observed) / mean Δ_null
```
The "tightening percentage" of the signal relative to the random null.
ε = 0 means no detection; ε = 1 means perfect alignment (all minima
exactly on zeros).

### Stop-on-fail

I commit to the brief's two stop conditions:
1. If (5, 5)_Riemann reproduction at N=1000 gives |z| < 5, STOP.
2. If 3 or more diagonal cells of Sub C-1 give |z| < 5, STOP (do not
   run Sub C-2).

### Budget

Estimated total: ~30 minutes for L-function zero finding (4 characters,
~ 250 zeros each at mpmath dps=15), plus probe + matching for all 16
cells at N=1000 (~ 1-2 minutes per cell, so ~ 20-30 minutes), plus
diagonal N=10⁴ extension and Sub C-2 N-sweep. Aim ~ 2 hours total.

If L-function zero finding overruns, fall back to fewer zeros per
character (~ 100 instead of 250) and document the reduced range.

## What I will NOT do

- Adjust the discriminant set after seeing data.
- Change the classification window W after seeing data.
- Run more cells than specified or relax the Bonferroni correction.
- Switch zero source mid-run.
- Read Mr Adversary's reviews of Papers 150/203.

## Files

- `chars.py` — character tables and verification.
- `l_zeros.py` — L(s, χ) zero finder using mpmath Hurwitz zeta.
- `probe.py` — parameterised probe |Re² − d·Im²|.
- `matrix_runner.py` — 4×4 matrix sweep + nulls.
- `effect_size.py` — Sub C-2 ε(N) sweep.
- `findings_paper_203_sub_c.md` — long-form report.
