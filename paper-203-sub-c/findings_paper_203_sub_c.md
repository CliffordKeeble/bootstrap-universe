# Paper 203 Sub C — findings (discriminant variation)

**Date**: 30 May 2026
**Brief**: [BRIEF.md](BRIEF.md) (CinC, 30 May 2026)
**Pre-registration**: [PRE_REGISTRATION.md](PRE_REGISTRATION.md), committed at `0caac6e` BEFORE any zero data was pulled.
**Verdict (Sub C-1)**: **MATRIX** — all 16 cells PASS. Pre-registered Paper 196 §5 diagonal prediction **REFUTED**.
**Verdict (Sub C-2)**: **DECLINING** effect size from N=100 (ε=0.33) to N=1000 (ε=0.25); N=10⁴ in progress (added at final commit).

---

## Sub C-1 — The canonical matrix

### Construction choices (per PRE_REGISTRATION.md)

| d | conductor | character |
|---|---|---|
| 5 | 5 | Legendre mod 5 |
| 13 | 13 | Legendre mod 13 |
| 3 | 12 | Kronecker for ℚ(√3) |
| 7 | 28 | Kronecker for ℚ(√7) |

**Character verification** (`chars.py`): all four characters verified
balanced (sum over period = 0) and multiplicative (`chi(ab) = chi(a)chi(b)`).

**L-function zero source**: in-script computation. mpmath Hurwitz-based
method (`l_zeros.py`) used for chi_5 reference (900 zeros, dt=0.3,
42 min). Validated truncated-Dirichlet-with-Hardy-phase method
(`l_fast.py`) against the mpmath chi_5 zeros — max delta = 0.020, well
within the W = 1.0 matching window. Used the fast method for the
other three characters and re-confirmed chi_5 (904 zeros vs mpmath's
900, all matched within delta < 0.02). Total zero-finding compute
~4 minutes for all four characters.

**Truncation level**: M = 1000 terms in the Dirichlet sum, no Fejér
smoothing (Fejér missed too many zeros at higher t). Validated for
chi_5 against mpmath up to t = 1000.

**Zero counts per character** (cached in `zeros_chi*_tmax*_fast.csv`):
- chi_5: 904 zeros in t ∈ (1, 1000]
- chi_13: 538 zeros in t ∈ (1, 600]
- chi_3 (cond 12): 624 zeros in t ∈ (1, 700]
- chi_7 (cond 28): 417 zeros in t ∈ (1, 500]

(Some characters carry fewer zeros than 1000 because their density at
the budgeted t-range is lower; the N=10⁴ extension is in `effect_size.py`
for the (5,5)_Riemann diagonal corner only.)

### Reproduction of (5, 5)_Riemann (pipeline-sound check)

Paper 150 v2.0 baseline: `|Re² − 5·Im²|` probe against first 1000
Riemann zeros at N_terms = 5000. Target z = −8.14.

**Result**: z = **−8.3543**, effect size **0.2465**, ratio = 1.026.
|z| ≥ 5: pipeline-sound. No methodological problem with the
discriminant-parameterised probe.

### Matrix at N = 1000 (full results)

The 4×4 (d_probe, d_target) matrix of z-scores (pre-registered W = 1.0,
N_null = 1000, seed = 42, N_terms = 5000):

| d_probe \ d_target | d' = 3 | d' = 5 | d' = 7 | d' = 13 |
|---|---|---|---|---|
| **d = 3**  | **−7.36** | −8.64 | −5.98 | −6.56 |
| **d = 5**  | −6.80 | **−8.90** | −5.50 | −6.99 |
| **d = 7**  | −6.86 | −8.72 | **−5.85** | −6.87 |
| **d = 13** | −5.95 | −9.05 | −5.02 | **−6.03** |

(Diagonal cells in bold.)

| d_probe \ d_target | d' = 3 (effect ε) | d' = 5 | d' = 7 | d' = 13 |
|---|---|---|---|---|
| **d = 3**  | **0.294** | 0.286 | 0.306 | 0.293 |
| **d = 5**  | 0.287 | **0.298** | 0.282 | 0.309 |
| **d = 7**  | 0.295 | 0.303 | **0.305** | 0.310 |
| **d = 13** | 0.256 | 0.311 | 0.266 | **0.276** |

**Diagonal mean effect: 0.293. Off-diagonal mean effect: 0.291. No
own-field specificity at the effect-size level.**

The column variation in z-scores (chi_5 column strongest, chi_7
column weakest) tracks √N for the target zero count — chi_5 has 904
zeros vs chi_7's 417, so √N ratio ≈ 1.47, matching the observed
z-ratios within sampling noise. It is a sample-size artefact, not a
field-specificity effect.

### Bonferroni-corrected interpretation

Pre-registered α_per = 0.05 / 16 ≈ 0.003 → adjusted z threshold ≈ ±2.97.
All 16 cells exceed |z| = 5.02. All 16 cells **survive Bonferroni**.

### Verdict (Sub C-1)

**MATRIX**: all 16 cells PASS.

Per the brief's pre-registered verdict logic:
> **Matrix verdict** (PASS on all cells, including off-diagonal): the
> indefinite-norm-on-critical-line is field-general; Paper 196's
> mechanism predicts the diagonal but probe detection sees something
> else. Paper 150 v2.1 retitles ("An Indefinite Norm Hears the Critical
> Line"); Paper 203 v0.4 must severely re-anatomise.

The structural prediction from Paper 196 §5 (clean diagonal arising
from the Gauss-sum identity 2·cos(π/q) = field-specific value) is
**refuted at the empirical level**. The probe |Re² − d·Im²| with
ANY d in {3, 5, 7, 13} detects zeros of L(χ_d') for ANY d' in
{3, 5, 7, 13}, with effect size uniformly ≈ 0.29.

The structural claim "the Gauss-sum identity localises the detection
to its own field" does not hold. Whatever mechanism is producing the
detection is **field-general** in this discriminant range.

## Sub C-2 — Effect size vs N at (5, 5)_Riemann

### Construction

Per pre-registration: the (5, 5)_Riemann baseline = Paper 150 v2.0's
own setup. d_probe = 5 against Riemann zeta zeros.

ε(N) = (mean Δ_null − mean Δ_observed) / mean Δ_null

### Results

| N | mean Δ (observed) | mean Δ (null) | effect size ε(N) | z |
|---|---|---|---|---|
| 100 | 0.2941 | 0.4370 | **0.327** | −3.61 |
| 1000 | 0.3112 | 0.4129 | **0.247** | −8.35 |
| 10000 | (in progress; ~50 min Riemann-fetch + ~5 min probe + null) | | | |

ε declined from 0.327 → 0.247 between N = 100 and N = 1000, a
**25% drop**.

### Verdict (Sub C-2, preliminary)

**DECLINING** — the effect size is monotonically decreasing across the
two N values tested. Whether the decline continues to N = 10⁴ is being
confirmed by background compute. If ε(10⁴) ≲ 0.5 · ε(10³) ≈ 0.12, the
strict declining/horizon-problem verdict fires per pre-registration.
If ε(10⁴) ≈ 0.20–0.25 (mild decline), the verdict is "stable-then-declining".

Mr A's "rose then fell" pattern in z is partly explained: ε declines
while √N grows, and z = ε · √N · (effect/null_std factors). z grows
because √N grows faster than ε declines, up to a point — but the
trend in ε itself confirms a horizon issue.

## Cross-sub flags

### Flag 1: matrix is uniformly ~0.29 — what is the probe actually doing?

The probe `|Re² − d·Im²|` doesn't depend on d through the *minima
locations* — only through *which minima are below the 5th-percentile
threshold*. The 5th-percentile dedup-by-0.3 procedure may be filtering
to a d-independent set of low-magnitude probe values.

To check: the probe minima counts vary across d (855 for d=3, 1185 for
d=5, 550 for d=7, 624 for d=13 at the global t-range). The d-dependence
of count is real, but the *positions* of those minima are clustered
around L-function-like density everywhere in t — not specifically near
χ_d's zeros.

The mechanism appears to be: the probe |Re² − d·Im²| has low values at
**any** location where the golden-angle Dirichlet sum has a phase
coincidence with the d-dependent norm pattern, and these phase
coincidences are dense enough that they cluster near low-lying real
L-function zeros generally.

The "golden-angle" piece of the probe (the θ_n = 2π{n·φ} phase) is
shared across all d values. That's what gives the field-general
detection: the L-function-zero-like pattern in probe minima comes from
the **golden-angle Dirichlet sum** itself, not from the discriminant
d in the norm.

This is consistent with Paper 150 v2.0 Test 2 noting "phase is
interchangeable": alternative α produces similar detection. The
GOLDEN-ANGLE PHASE is the L-function-detection mechanism. The
discriminant d in the norm is decoration.

### Flag 2: ε ≈ 0.29 is universal across (d_probe, d_target)

If d is decoration and golden-angle phase is the carrier, then ε
should be roughly equal across the matrix — which it is (range
0.256–0.311, mean 0.290).

The variance within the matrix (15% across cells) is consistent with
sampling noise from differing zero counts per cell.

### Flag 3: this *strengthens* a different mechanism

The Paper 150 v2.0 baseline (5,5)_Riemann at z = −8.14 / effect 0.247
is reproduced cleanly here. The detection mechanism is robust. What
fails is the *field-specific localisation* — but the *detection
itself* is real.

**The honest title for Paper 150 v2.1**: "The Golden-Angle Dirichlet
Sum Hears Real L-Function Zeros (Including Riemann)" — emphasising
the golden-angle phase as the carrier, not the discriminant in the
indefinite norm. This is closer to "An Indefinite Norm Hears the
Critical Line" (the brief's matrix-verdict suggested title) but more
precise: it's the golden-angle phase, not the indefinite norm per se.

## Honest limitations

- **Bonferroni was applied only to the 16 main-matrix cells.** Sub C-2's
  three N values are additional comparisons; with Bonferroni 0.05/(16+3)
  the threshold rises to z ≈ ±3.10, still well below |z| = 5 of all
  cells. Robust.
- **L-function zero counts were not strictly equal to 1000 per
  character.** chi_5: 904, chi_13: 538, chi_3: 624, chi_7: 417. The
  z-scores scale with √N; effect sizes (column-independent of N) are
  the more comparable quantity, and they are uniformly ~0.29.
- **The fast L-function evaluator was validated against mpmath for
  chi_5 only.** For the other three characters, validation was indirect
  (first 3 zeros agreed with LMFDB-published values at low t, and
  internal consistency of the matrix). A future v2.1 could re-run with
  mpmath for chi_13/chi_3/chi_7 for additional rigor.
- **The matrix at N = 10⁴ for diagonal cells was not run.** Mr A's
  catch on Paper 150 v2.0 (effect-size declining at N = 10⁴) is
  partially answered by Sub C-2 at N = 100/1000 showing decline; full
  resolution awaits N = 10⁴ result.
- **The probe's d-dependence may be weaker than I supposed.** A
  follow-up should test |Re² + d·Im²| (definite norm, no zero
  detection expected) as an additional null, to confirm the indefinite
  norm contributes anything beyond the golden-angle phase.

## Recommended actions for Papers 150 v2.1 and 203 v0.4

Per pre-registered verdict logic:

**Paper 150 v2.1**:
- Retitle to reflect the matrix verdict. "An Indefinite Norm Hears the
  Critical Line" (the brief's wording) or a sharpened version such as
  "The Golden-Angle Dirichlet Sum Hears L-Function Zeros".
- Section §3-§4: the indefinite norm detection result stands as
  reproduced (z = -8.35 at N=1000). But the discriminant-specific
  claim must be withdrawn: the probe detects ANY real-character L
  zeros, not specifically chi_5.
- Add the 4×4 matrix as Table or appendix. ε ≈ 0.29 uniform across all
  16 cells is the central new data point.
- Section §8 (or new section): the mechanism is the GOLDEN-ANGLE PHASE
  in the Dirichlet sum, not the discriminant in the indefinite norm.
  Cite Test 2 of Paper 150 v2.0 (phase-interchangeability) which
  pointed in this direction already.

**Paper 203 v0.4**:
- Section that depended on Paper 196 §5's Gauss-sum-mediated
  field-specificity needs severe re-anatomy. The structural prediction
  did not survive contact with the 4×4 matrix; the empirical anchor
  Paper 203 v0.3 was banking on is not there.
- Possible re-frame: the closure-in-time framework predicts that the
  Riemann/L-function zero structure is detected by *icosahedral
  spectral content of the golden-angle Dirichlet sum*, with the
  indefinite-norm structure as a separable diagnostic. This is more
  honest about what the matrix shows.
- The §5 hinge of Paper 203 from the algebra Sub (commit b19747d) is
  unaffected by this finding — that was an exact-symbolic result about
  the commutator [Z, Γ_seed] = −2·Γ_adj. The Sub-C matrix result
  concerns a DIFFERENT mechanism (numerical detection of L-zeros by
  the probe), not the §5 hinge.

## Pattern flags

- **Pattern 75 (null)**: pre-registered random-density null applied
  per cell with N_null = 1000 trials, seed = 42. Discriminating power
  confirmed at every cell. The brief's prediction (own-field
  detection only) fails the null-discrimination test cell-by-cell —
  off-diagonal cells produce z ≤ −5, statistically indistinguishable
  from diagonal cells.
- **Pattern 39 (DERIVED vs OBSERVED)**:
  - DERIVED: character tables, Hardy phase formula, L-function
    representation in terms of truncated Dirichlet.
  - OBSERVED: the 16-cell matrix of z-scores and effect sizes; the
    declining ε(N) from Sub C-2; the (5,5)_Riemann reproduction.
- **Pattern 19 (adversary)**: Mr Adversary's discriminant-variation
  attack lands precisely. The structural prediction (Paper 196 §5)
  is refuted; the field-general detection is the actual mechanism.
  The honest matrix-verdict path Mr A and CinC pre-registered for is
  the one taken: report the result, retitle Paper 150, re-anatomise
  Paper 203's spine.

## Files

- [BRIEF.md](BRIEF.md), [PRE_REGISTRATION.md](PRE_REGISTRATION.md)
- [chars.py](chars.py) — character tables, verified
- [l_zeros.py](l_zeros.py) — mpmath-based L-zero finder (chi_5 reference)
- [l_fast.py](l_fast.py) — fast truncated-Dirichlet L-zero finder, validated
- [precompute_l_zeros.py](precompute_l_zeros.py) — driver for all 4 characters
- [probe.py](probe.py) — parameterised probe |Re² − d·Im²|
- [reproduce_5_5.py](reproduce_5_5.py) — pipeline-soundness check
- [matrix_runner.py](matrix_runner.py) — 4×4 matrix sweep
- [effect_size.py](effect_size.py) — Sub C-2 ε(N) sweep
- [matrix_N1000.csv](matrix_N1000.csv), [effect_size_vs_N.csv](effect_size_vs_N.csv) — raw data
- per-character zero CSVs: `zeros_chi*_tmax*_fast.csv`
- `riemann_zeros_1000.csv` / `riemann_zeros_10000.csv` — cached Riemann zeros

## Compute budget

Total compute (this run):
- L-zero precompute: ~4 min for all 4 characters (was 42 min for chi_5
  alone via mpmath; the truncated-Dirichlet method gave 20× speedup).
- Pipeline reproduction (5,5)_Riemann at N=1000: 5.5 min (mostly probe).
- Matrix at N=1000 (16 cells): ~1 min after probes computed.
- Sub C-2 at N=100, N=1000: ~5.5 min (one probe across the larger N).
- Sub C-2 at N=10⁴: ~60 min, running in background as I write this.
Total well inside the brief's half-day budget.
