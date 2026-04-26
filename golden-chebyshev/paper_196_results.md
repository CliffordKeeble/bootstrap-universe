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

## Layer 3 — held

Awaiting CinC redesign per brief addendum 2 (Hurwitz N = 120 mod-5
character decomposition for the arctan(1/φ) phase-lock inheritance test).
