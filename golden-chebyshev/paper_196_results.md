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

## Layer 3 — held

Awaiting CinC redesign per brief addendum 2 (Hurwitz N = 120 mod-5
character decomposition for the arctan(1/φ) phase-lock inheritance test).
