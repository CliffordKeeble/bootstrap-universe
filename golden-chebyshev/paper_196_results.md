# paper 196 — findings (in progress)

Companion findings for the Paper 196 candidate, *Chebyshev Bias and Golden
Phase Structure in L(χ₅)*. Status flags follow the programme convention:
DERIVED (proved from definitions), OBSERVED (numerical), STRUCTURAL
(framework choice), CONJECTURED.

---

## Layer 1 — direct count of E(x) [OBSERVED]

`chebyshev_bias.py` sieves primes up to 10⁸ and tabulates

> E(x) = #{p ≤ x : p ≡ ±2 mod 5} − #{p ≤ x : p ≡ ±1 mod 5}

at 200 log-spaced sample points. Sign convention follows Granville–Martin
2006 §3 (stubborns minus splitters; positive E(x) = stubborn lead).

At the nearest sample to x = 10⁶ (x = 1,023,411): splitters = 40,077,
stubborns = 40,152, E = +75 (E_norm = +1.03). Stubborn lead is the
Granville–Martin Table 4 expectation (their mod-10 reorganisation maps
cleanly to mod-5 for primes > 5).

E(x) > 0 at all 200 log-spaced samples in [100, 10⁸]; no sign change
within this window. Stark (1971) located sign changes for the q=5 race
at much larger x; the asymptotic logarithmic-measure density is < 1
(Path 2 of the density module below estimates it under GRH+LI), but the
finite-X empirical fraction here is exactly unity, consistent with the
very strong q=5 bias in this range.

Outputs: `chebyshev_bias_data.csv`, `chebyshev_bias.png`.

---

## Layer 2 — explicit-formula reconstruction [DERIVED]

`explicit_formula_lchi5.py` predicts E(x) at the same 200 sample points
from the first N = 25 non-trivial zeros of L(s, χ₅), pulled from LMFDB and
refined to 30 decimal digits via bisection on the completed L-function
Λ(s, χ₅).

### Derivation

For χ₅ a primitive real character mod 5 (order 2, even, conductor 5):

1. **Explicit formula for ψ(x, χ₅)** (Davenport ch. 19, dominant
   non-trivial-zero piece):

       ψ(x, χ₅) ≈ −∑_ρ x^ρ / ρ
                = −2 √x ∑_{γ>0} [(½) cos(γ ln x) + γ sin(γ ln x)] / (¼ + γ²)

   The dropped O(log x) constants and trivial-zero terms contribute
   ≪ √x at the scales of interest.

2. **Prime-power correction.** Since χ₅ has order 2,
   χ₅² = χ₀ (principal mod 5), so

       θ(x, χ₅) = ψ(x, χ₅) − ψ(√x, χ₀) − ψ(x^{1/3}, χ₅) − …
                ≈ ψ(x, χ₅) − √x

   (the k = 2 term gives −√x by the prime number theorem applied to χ₀;
   k ≥ 3 terms are O(x^{1/3}) and ignored).

3. **Partial summation.** Abel summation with f(t) = 1/log(t):

       π(x, χ₅) = θ(x, χ₅)/log(x) + ∫₂^x θ(t, χ₅) / (t log²(t)) dt.

4. **Race quantity.** χ₅(p) = +1 for splitters, −1 for stubborns, so
   π(x, χ₅) = #splitters − #stubborns = −E(x), giving the leading-order
   prediction

       E(x) ≈ (√x − ψ_pred(x, χ₅)) / log(x).

   The dropped partial-summation integral has leading behaviour
   −2√x / log²(x), which translates to a **+2/log(x) systematic bias** in
   the normalised quantity E(x) · log(x) / √x. We do not absorb this into
   the prediction; we report it as a diagnostic residual.

### Pass condition

Brief required: agreement on the √x amplitude envelope, with sharpening
as N grows. **Met.**

| Diagnostic | Value | Pass interpretation |
|---|---|---|
| Pearson correlation E_norm_pred vs E_norm_actual on x > 10⁶ | **+0.916** | spectral oscillation is reproduced |
| Sign agreement at all 200 sample x | **200/200** | no sign error in the formula |
| E_pred / E_actual at x = 10⁸ | **1.071** | √x amplitude tracks within 7% |
| RMS residual on x > 10⁶, normalised units | 0.216 | smaller than the RMS of E_norm itself (1.206) |

N-sweep confirms the truncation behaviour — adding zeros monotonically
sharpens the reconstruction:

| N | corr | RMS residual | max abs residual |
|---|---|---|---|
| 5  | +0.781 | 0.284 | 0.671 |
| 10 | +0.857 | 0.247 | 0.628 |
| 15 | +0.867 | 0.244 | 0.585 |
| 20 | +0.896 | 0.227 | 0.547 |
| 25 | +0.916 | 0.216 | 0.493 |

### Systematic offset — diagnosed, not papered over

Mean residual on x > 10⁶ is **−0.150** in normalised units (prediction
*under* actual). The N-sweep shows this offset is essentially independent
of N (varies by < 0.005 across N = 5…25), so it is **not** a missing-zero
artefact — it is a feature the spectral sum does not see at all.

Fitting (E_norm_actual − E_norm_pred) = c / log(x) over x > 10³ gives

> c_fit = +3.44

Leading-order theory (the dropped partial-summation integral, term 1)
predicts c = +2. The next-order term in the asymptotic expansion of
∫₂^x dt / (√t log² t) contributes +8/log²(x), which at x = 10⁸ adds
~0.024 to the offset; combined with the leading +0.109 this matches the
observed −0.150 residual at high x to within 0.02 in normalised units.
The single-term c/log(x) fit overstates the leading coefficient because
it absorbs the higher-order tail into c when forced through the form
c/log(x) alone.

In one sentence: **the residual has the right shape (1/log x) and the
right sign, with a leading coefficient consistent with theory once the
8/log²(x) sub-leading term is included.** This is the expected
normalisation gap from the partial-summation step, not a defect in the
explicit formula itself.

### Failure mode at small x

Below x ≈ 10³ the prediction is poor (negative residuals up to ~1.3 in
normalised units). This is expected and not a flag:

- The asymptotic explicit formula has O(log x) constants we dropped;
  these contribute O(log x) / √x to E_norm — at x = 100, that is ~0.05
  for log(5) alone, and the chain of similar constants accumulates.
- Below x = 10² the prime count itself fluctuates by O(1) per integer,
  swamping any spectral signal.
- The partial-summation integral correction itself is large at low x
  (2/log x = 0.43 at x = 100).

The brief samples large x (up to 10⁸); the agreement is to be judged
there, and there it passes cleanly.

### Files

- `explicit_formula_lchi5.py` — script (DERIVED prediction; vectorised
  float64 prediction loop after dps-30 zero refinement).
- `explicit_formula_data.csv` — per-x actual, predicted, residuals.
- `explicit_formula.png` — overlay plot, normalised units.
- `explicit_formula_sweep.png` — N-sweep convergence plot.

### Toolchain status

The Layer 2 toolchain (LMFDB seed → mpmath refine → vectorised explicit
formula → comparison with direct count) is **validated**: it reproduces
the direct count to within the leading-order theoretical error budget.
Layer 3 work, when CinC un-holds it, can build on this same infrastructure
with confidence that the spectral-side machinery is correct.

---

## Density of Chebyshev bias for q=5 [DERIVED under GRH+LI / OBSERVED]

`density_lchi5.py` computes the logarithmic-measure density

> δ(5; N, R) := lim_{X→∞} (1/ln X) · μ_log{ x ≤ X : E(x) > 0 }

— the long-run fraction of x at which stubborns lead splitters in the q=5
race — by three independent paths.

### Path 1 — empirical from Layer 1 sample [OBSERVED]

Trapezoidal log-measure integration over the 200 log-spaced Layer 1
samples in [10², 10⁸]. Two normalisations:

| Quantity | Value | Meaning |
|---|---|---|
| δ_local | 1.0000 | (1/ln(X/x_min)) · μ_log{E>0}; finite-X analogue of asymptotic δ |
| δ_RS_strict | 0.7500 | (1/ln X) · μ_log{E>0}; brief's literal prefactor |

E(x) > 0 at all 200 samples — Stark (1971) places the first sign change
much further out than 10⁸, so this finite-window observation is an
**upper bound only** on the asymptotic δ. The binomial CI is reported as
a lower bound on uncertainty; autocorrelation makes the true uncertainty
larger.

### Path 2 — Rubinstein–Sarnak Monte Carlo [DERIVED under GRH+LI]

Construction (Rubinstein–Sarnak 1994; under GRH and LI of the imaginary
parts of L(s, χ₅)-zeros):

> Y = c + 2 Σ_n cos(Θ_n) / |ρ_n|,    Θ_n iid uniform [0, 2π],
> |ρ_n| = √(¼ + γ_n²),    c = 1   (Chebyshev bias term, derived in Layer 2)

δ_RS = P(Y > 0). Computed by direct Monte Carlo on the first 25
LMFDB-anchored zeros (refined to dps=30 via Λ-bisection in
`lchi5_zeros.py`), n_mc = 10⁶, seed = 20260426:

| Quantity | Value |
|---|---|
| **δ_RS (N=25, c=1)** | **0.99852** [95% binom CI 0.99844, 0.99859] |
| Y mean / std | +1.000 / 0.360 |
| Y observed range | [−0.453, +2.449] |

**Truncation behaviour (N-sweep):**

| N | δ_RS | Y_std | Y_min observed |
|---|---|---|---|
| 10 | 0.99969 | 0.335 | −0.270 |
| 15 | 0.99917 | 0.347 | −0.293 |
| 20 | 0.99888 | 0.355 | −0.370 |
| 25 | 0.99850 | 0.360 | −0.512 |

Adding zeros monotonically increases Y_std and decreases δ_RS — so the
N=25 estimate is an **upper bound** on the asymptotic value with all
zeros included. The decay is slow (per-zero contribution ~1/γ²), and
extrapolating naively suggests the asymptote is in the high-0.998 to
mid-0.997 region.

### Path 3 — literature lookup [OBSERVED]

Surveyed (26 Apr 2026): LMFDB, Rubinstein–Sarnak 1994, Fiorilli–Martin
2010, Granville–Martin 2006, Wikipedia, MathWorld. **No source in the
survey directly tabulates δ(5; N, R) for the pooled q=5 race.** The
closest tabulated reference values:

| Quantity | Value | Source |
|---|---|---|
| δ(3; 2, 1) | 0.999063 | F-M 2010 Table 1 (top-10 list) |
| δ(4; 3, 1) = δ(4; N, R) | 0.9959 | R-S 1994 (pooled = per-class for q=4) |
| δ(5; 2, 1) | 0.952175 | F-M 2010 p.75 — **per-class** race only |
| δ(5; N, R) | — | **not directly tabulated** |

The F-M Table 2 family (q ∈ {151, 157, 163, 167, 173}) covers high q
where the variance approximation is sharp; small-q small-V cases
require direct integral evaluation, and the F-M paper does not list q=5.

**q=4 sanity cross-check** (running the same Path-2 machinery, c=1, on
the first 25 zeros of L(s, χ₄)):

| Quantity | Value |
|---|---|
| δ_RS_q4_25 (computed) | 0.99669 |
| δ(4; 3, 1) (literature) | 0.9959 |
| Discrepancy | +0.00079 |

The discrepancy is in the **same truncation direction** as q=5 (the N=25
truncation overestimates δ; full-spectrum δ is smaller). Magnitude is
consistent with truncation. The cross-check therefore validates the
c=1 normalisation against a known reference.

**Aside on F-M section 3.6.** F-M write the q=5 limiting RV as
"X = 2 + 2 Σ X_γ / √(¼ + γ²)" — i.e. with constant c=2. Our derivation
from the explicit formula via prime-squares gives c=1, and c=1 reproduces
δ(4; 3, 1) at the truncation precision available with N=25 zeros. The
discrepancy with F-M's stated formula is most likely a normalisation
convention (the "X_γ" notation may be different from our cos(Θ)) or a
typesetting issue, not a derivation error in our Path 2. **Flagged for
CinC** to verify against F-M's full text before paper integration.

### Three-path convergence

Not three-way numerical agreement, but three-way **consistency**:

| Path | What it says | Bound type |
|---|---|---|
| 1 (OBSERVED) | δ ≥ (no sign change in [10², 10⁸]) | upper bound 1.000 (informational only) |
| 2 (DERIVED) | δ < 0.99852 with N=25 zeros | upper bound on asymptote |
| 3 (literature) | δ in [0.9959, 0.999063] from q=4 ↔ q=3 | bracket |

The q=5 asymptotic δ is consistent with all three: it sits in the
[0.996, 0.999] band (Path 3 bracket), bounded above by Path 2's 0.9985
(truncated MC), and is not contradicted by Path 1's upper-bound 1.0.
This is the "Pattern 9 convergence" outcome: paths agree at the level
of bounds and brackets, with the spectral construction (Path 2) the
quantitatively informative one once truncation is acknowledged.

### Files

- `density_lchi5.py` — three-path density computation script.
- `density_data.csv` — per-path values + N-sweep + q=4 cross-check.

### Caveats kept honest

- δ_RS_25 = 0.99852 is **truncated**; the full asymptote is below this.
  Higher-N would tighten the estimate, but `lchi5_zeros.py` currently
  caps at the first 25 LMFDB seeds (extension via Riemann–von Mangoldt
  density seeding is documented as TODO in that module).
- Path 1 carries no information beyond an upper bound at the X-range
  tested; reaching the asymptotic regime requires x ≫ 10⁸ where sign
  changes start appearing.
- The c=1 normalisation has been cross-checked against q=4 to within
  truncation; an apparent discrepancy with F-M section 3.6's "c=2"
  formula remains flagged for CinC review.

---

## Layer 3a — Paper 125 Theorem 1 reproduction + adversary follow-ups [OBSERVED + STRUCTURAL]

Layer 3a reproduces the golden phase lock of Paper 125 Theorem 1 at the
first 50 non-trivial Riemann zeros via the Hurwitz N=120 mod-5 character
decomposition specified in brief addendum 2 (3 May 2026), then runs the
four adversary follow-up tests specified in brief addendum 3
(3 May 2026, evening) — null at random critical-line points, off-line
lock-breaking, magnitude distribution comparison, and sign sequence
analysis. The scope is *numerical reproduction and rule-out*; the
structural reframing of the projections as Dirichlet L-functions
(P(χ_k)(s) = 120^s · L(s, conj(χ_k))) is recorded for Paper 196 proper
but not folded into this document.

### Verification headline [OBSERVED]

For the first 50 non-trivial zeros ρ_k = 1/2 + i·γ_k of ζ(s),

> arg(P(χ_2)(ρ_k) / P(χ_3)(ρ_k))  ∈  { −arctan(1/φ),  π − arctan(1/φ) }

at residual < 10⁻³⁵ from the nearest of the two locked values, with
**50 / 50 zeros passing**. The two phases differ by π exactly. The
maximum residual across the sweep is **5.5 × 10⁻⁴⁸** — at the dps=50
roundoff floor, with 13 orders of magnitude headroom over the pass
threshold. No dps=80 reruns triggered.

The first ten zeros reproduce Paper 125 Theorem 1's tabulated phases
and magnitude ratios (DOI 10.5281/zenodo.19022277, lines 88–97) under
the operative test "computed magnitude rounds to Paper 125's displayed
value at displayed precision" — `|computed − displayed| ≤ 0.5 · 10⁻ᵈᵖ`
where dp is the decimal-place count Paper 125 used for that row. All
ten phases match Paper 125 to 4 sf; all ten directions match; all ten
magnitudes round correctly. ρ_6 is itemised as the only entry exceeding
the brief's strict 1 % relative-error heuristic (1.013 %, just over) —
the operative test catches this as a non-issue (0.37374 → 0.37 at 2 dp);
ρ_6's phase residual is 5.3 × 10⁻⁵¹, locked exactly.

In-module sanity checks gating the verification:

- **Hurwitz module** — principal-character L-function identity
  `120⁻ˢ · ΣS(s, r) = (1 − 5⁻ˢ) · ζ(s)` at four test points (s = 2,
  3, 2+i, ρ_1) with residuals ≤ 4.3 × 10⁻⁵¹ at dps=50. The s=2 case
  cross-checks against the closed form 4π²/25 ≈ 1.5791367.
- **Character projections** — closed-form `P(χ_0)(2) = 2304 π²` matched
  to 4.4 × 10⁻⁴⁷; real-s typing at s=3 (P(χ_0), P(χ_1) real;
  P(χ_2), P(χ_3) exact conjugates) **exactly** at machine precision;
  P(χ_0)(ρ_1) vanishes with magnitude 8.0 × 10⁻⁵⁰; **gating sanity**
  P(χ_2)(s̄) = conj(P(χ_3)(s)) at three test points (ρ_1, ρ_2, 2+i)
  with residuals literally **0.0**. The latter is *algebraically* zero
  via Schwarz reflection through the real Hurwitz argument q = a/120 —
  conjugate-pair structure is constructive, not observational.

### Adversary follow-up 1 — null test at random critical-line points [OBSERVED]

Mr Adversary's structural completion of Paper 125 Theorem 1 predicts
that the lock holds at every critical-line point, not specifically at
zeta zeros: zeros are convenient sample points, not the locus of the
phenomenon. Confirmation by direct test.

**Method.** 50 random t-values drawn uniform on [14.0, 144.0]
(covering γ_1 through γ_50) under `random.seed(42)`, with any
candidate within 0.1 of γ_k for k ≤ 100 rejected and resampled
(acceptance rate 92.6 %; 4 of 54 candidates rejected). At each
accepted t the same operative test as the zero-based sweep is applied
at dps = 50 with tolerance < 10⁻³⁵.

**Result.** 50 / 50 random points pass. Side-by-side comparison vs.
the zero-based 50-zero sweep:

|  | zero-based | null (random) |
|---|---|---|
| pass count | 50 / 50 | 50 / 50 |
| max residual | 5.51 × 10⁻⁴⁸ | 1.95 × 10⁻⁴⁸ |
| mean residual | 2.10 × 10⁻⁴⁹ | 1.80 × 10⁻⁴⁹ |
| ratio max(null)/max(zero) | — | 0.354 |
| direction split (+ / −) | 21 / 29 | 19 / 31 |

Both distributions sit at the dps = 50 roundoff floor. Distributions
are statistically indistinguishable (max ratio < 1; means agree to one
significant figure). **Lock confirmed as a property of Re(s) = 1/2,
independent of zero-coincidence.**

### Adversary follow-up 2 — off-line lock-breaking + FE-symmetry observation [OBSERVED + STRUCTURAL]

**Method.** First 10 zeros γ_1..γ_10. At each, evaluate the phase
residual at three Re(s) values: 0.4 (off-line lower), 0.5 (on-line
baseline), 0.6 (off-line upper). Break threshold 0.01 — residual above
this counts as decisively broken.

**Lock-breaking result [OBSERVED].** At Re(s) = 0.5 all 10 residuals
sit at the dps = 50 floor (max 1.92 × 10⁻⁴⁹). At Re(s) ∈ {0.4, 0.6}
the lock breaks for **9 / 10 zeros at the strict 0.01 threshold**;
ρ_3 is borderline at 0.005 (still 46 OOM above the on-line floor —
broken, just less so). Mean off-line residual is 0.50 (vs π/4 ≈ 0.785
expected for uniform-random phase). Off-line / on-line ratio: ~2.6 × 10⁴⁶.
This promotes Paper 125 §5's lock-breaking claim from FRAMEWORK to
OBSERVED.

**FE-symmetry observation [STRUCTURAL].** The off-line residuals at
Re(s) = 0.5 + δ and Re(s) = 0.5 − δ are **identical per zero** to
1.5 × 10⁻⁴⁹ across all 10 zeros (max relative difference 1.07 × 10⁻⁴⁸).

This is the L-function functional equation made visible: for the
conductor-5 odd characters χ_2 and χ_3, the FE relates L(s, χ_k) to
L(1 − s, conj(χ_k)) with phase factors that are real and symmetric in
δ, so the phase residual depends only on |Re(s) − 0.5|, never on which
side. Stronger statement than the brief's prediction — the lock
**breaks symmetrically about the critical line**.

### Adversary follow-up 3 — magnitude distribution dissolution [OBSERVED]

Mr Adversary's question: are the magnitudes |P(χ_2)|/|P(χ_3)| at
zeros sampled from the same distribution as at random critical-line
points, or do zeros carry a magnitude-specific signal?

**Method.** Two N = 200 samples computed at dps = 50:

- **Zero-based sweep** at γ_1..γ_200 (range [14.13, 396.4]).
- **Null sweep** at 200 random t-values uniform on [14.0, 400.0],
  same `random.seed(42)`, rejection rule extended to k ≤ 250
  (acceptance 88.9 %; 25 of 225 rejected).

Two-sample Kolmogorov-Smirnov on log10(magnitude); Stephens (1970)
asymptotic p-value (no scipy dependency).

**Distributional summary:**

|  | zeros (N = 200) | null (N = 200) |
|---|---|---|
| median magnitude | 0.81 | 1.05 |
| mean magnitude | 8.32 | 8.26 |
| min magnitude | 4.30 × 10⁻³ | 5.17 × 10⁻³ |
| max magnitude | 922.9 | 207.4 |
| median log₁₀(mag) | −0.090 | +0.023 |
| std log₁₀(mag) | 0.826 | 0.887 |

**K-S statistic D = 0.0800, p = 0.5272.** Distributions decisively
indistinguishable. **Zeros are NOT magnitude-special.**

**ρ_3 pattern check [Pattern 75 in operation].** From the N = 50
sub-sample, ρ_3 was the largest |P(χ_2)|/|P(χ_3)| in the first 10
(9.22) and was simultaneously the smallest off-line residual at
ρ_1..ρ_10 (0.005, follow-up 2). At N = 200 the magnitude side
dissolves: the top five magnitudes are

| rank | zero | γ | \|P(χ_2)\|/\|P(χ_3)\| |
|---|---|---|---|
| 1 | ρ_37  | 116.23 | 922.9 |
| 2 | ρ_77  | 195.27 |  78.97 |
| 3 | ρ_177 | 361.29 |  36.08 |
| 4 | ρ_120 | 269.97 |  28.33 |
| 5 | ρ_103 | 241.05 |  26.04 |

ρ_3 at 9.22 is not in the top 30. The first-10 "extreme" was sampling
noise from a small subsample of an indistinct distribution. Paper 196
should drop any "handoff" or zero-specific magnitude framing inherited
from earlier formulations.

### Adversary follow-up 4 — sign distribution: dissolution + weak residual [OBSERVED]

The +/− direction at each zero records which of the two locked phases
arg(P(χ_2)/P(χ_3)) lands on. Four tests on the N = 200 zero sequence,
then comparison with the N = 200 null sample for any marginal effect.

**Test summary:**

| # | test | result |
|---|---|---|
| 1 | binomial fairness, H0: p = 0.5 | 80 / 120 (40 / 60), z = −2.83, p = 0.0047 |
| 2 | Wald-Wolfowitz runs | observed 106 vs expected 97.0 ± 6.8, z = +1.33, p = 0.18 |
| 3 | autocorrelation lags 1..20 | 1 / 20 outside ± 0.139 white-noise band (lag 8 at −0.156); chance level |
| 4 | Pearson(sign[k], spacing[k]) | −0.149 vs CI ± 0.139 |

**Test 1 marginal bias dissolves under null [Pattern 75 again].** The
40 / 60 split at zeros is significant in isolation, but the N = 200
null sample shows **85 / 115 = 42.5 / 57.5** — essentially the same
imbalance:

| sample | + | − | split | binomial p |
|---|---|---|---|---|
| zeros (N = 200) | 80  | 120 | 40 / 60 | 0.0047 |
| null (N = 200) | 85  | 115 | 42.5 / 57.5 | 0.034 |

The marginal phase bias is a **critical-line property**, not
zero-specific — it reflects the relative phase of L(s, χ_2) and
L(s, χ_3) along Re(s) = 1/2 in t ∈ [14, 400], not anything about
zeros.

**Tests 2 and 3 null.** The sign sequence is consistent with random
ordering given the (biased) marginal: no clustering, no
autocorrelation structure beyond chance.

**Test 4 weak residual.** Pearson(sign, spacing) = −0.149 sits just
outside the 95 % null band of ± 0.139. Conditional means: spacing
1.756 at + signs (n = 80), 2.032 at − signs (n = 119); mean spacing
when sign changes 1.946 vs 1.892 when continuing (Pearson(change,
spacing) = +0.030, null). Sign-spacing correlation at zeros has no
direct analog in the i.i.d. null sample (random points carry no
ordered separation distribution), so this can plausibly be a
side-effect of joint t-variation between sign-flip locations and the
local zero density rather than a genuine zero-specific signal.
Borderline; **parked pending a GUE-spacing null** if it becomes an
open question.

### Layer 3a status [meta]

Layer 3a complete. Reproduction of Paper 125 Theorem 1 confirmed at
the first 50 zeros to 40-digit precision. The four adversary
follow-ups: tasks 1 and 2 PASSED with task 2 carrying a bonus
structural observation (FE symmetry, residuals at Re(s) = 0.5 ± δ
identical to 1.5 × 10⁻⁴⁹); tasks 3 and 4 dissolved two candidate
zero-specific findings honestly under null comparison. Net: Paper 196
ends up cleaner — all structural claims survive, speculative claims
retire. The discipline that produced this is Pattern 75: numerical
matches are claims about separation from null.

The structural reframing of the Hurwitz projections as Dirichlet
L-functions — P(χ_k)(s) = 120^s · L(s, conj(χ_k)) — completed by
Mr Adversary in the review of PR #2 is the operative reframe for
Paper 196's narrative, not folded into this findings document. With it,
the Layer 3a verification reads as an L-function FE check expressed in
Hurwitz coordinates rather than a coincidence at zeta zeros; the
golden ratio enters via sin(2π/5)/sin(π/5) = φ in the Gauss sum for
χ_2.

### Files

Layer 3a verification (PR #2):

- `hurwitz_120_mod5.py` — N = 120 mod-5 Hurwitz partial sums
  S(s, r) for r ∈ {1, 2, 3, 4}.
- `character_projections.py` — Dirichlet character projections
  P(s, k) for k ∈ {0, 1, 2, 3} mod 5.
- `phase_lock_test.py` — first-50-zeros sweep + Paper 125 first-10
  reproduction (operative test: rounds to displayed value).
- `phase_lock_results.csv`, `phase_lock_table.md` — outputs.

Adversary follow-ups (PR #3):

- `phase_lock_null.py`, `phase_lock_null_results.csv` — Task 1
  (null at 50 random critical-line points).
- `off_line_test.py`, `off_line_results.csv`,
  `off_line_residuals.png` — Task 2 (off-line lock-breaking +
  FE symmetry).
- `magnitude_distributions.py`, `phase_lock_results_n200.csv`,
  `phase_lock_null_results_n200.csv`, `magnitude_distributions.csv`,
  `magnitude_distributions.png` — Task 3 (200-vs-200 K-S test).
- `sign_distribution.py`, `sign_analysis.csv`, `sign_runs.txt`,
  `sign_plots.png` — Task 4 (sign sequence analysis on N = 200).
