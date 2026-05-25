# Paper 190 Phase 1 v0.2 — Reconnaissance Findings Report

**Bootstrap Universe Programme**
**Role:** Mr Code (computational, scope-limited per Pattern 97)
**Briefed by:** CinC
**Date:** 2026-05-20
**Pre-registration commit:** `971938b` (v0.2 brief committed unchanged before any computation)
**Verdict:** **HALT-NULL**

---

## 1. Pre-registration record

The v0.2 brief `Mr_Code_Brief_Paper_190_Phase_1_v0_2.md` was committed
unchanged to git at hash `971938b` before any code ran. The brief supersedes
v0.1 ("γ as harmonic-basis shadow of 1/√3"), which is explicitly shelved.

**Hypothesis space (locked target):**

- **H1 (STRONG):** One or more of `{γ, γ_0(ℚ(√5)), γ_0(ℚ(√−3)), γ_{2I-spectral}}`
  shows a structurally meaningful CF relationship to the others or to
  `{1/√3, log φ, π/(3√3), 2 log φ/√5}`. Equivalently: γ's CF prefix carries
  structure predictable from Euler–Maclaurin / Bernoulli machinery beyond
  what a random Stieltjes-like constant of comparable irrationality measure
  would show.
- **H0:** γ-CF behaves indistinguishably from a random K-L Stieltjes-like
  constant; Stieltjes constants of the candidate fields show no CF-level
  relationship to γ beyond chance; Bernoulli-derived prefix prediction for
  γ-CF does not exceed null at p < 0.05.

**Pre-registration acknowledgement:** Mr Code read v0.1 + the surrounding
CinC framing (the conversation that produced the brief) before this v0.2
re-brief. The hypothesis space above was locked before any computation
ran. v0.1's 1/√3-specific test is explicitly shelved, not folded in.

**Sign convention (Stieltjes constants):** `γ_0(K) =` constant term of
`[ζ_K(s) − R/(s−1)]` at `s = 1`, no sign flip (Mathar convention).

---

## 2. Task 1 — Bernoulli baseline (verification)

**PASS** at the ≥ 50-digit verification target.

| metric | value |
|---|---|
| computed γ | `0.5772156649015328606065120900824024310421593359399235988` |
| mpmath reference | `0.5772156649015328606065120900824024310421593359399235988` |
| signed error | 3.89 × 10⁻⁶¹ |
| correct digits | ~60.4 |

n = 10⁶, K = 10, 60 dps working precision. EM truncation residual at K = 10
is bounded by `|B_22/(22·n^22)| ≈ 2.8 × 10⁻¹³⁰`, so the observed error is
working-precision rounding accumulation in H_n, not EM truncation. First
10 asymptotic terms `−B_{2k}/(2k·n^{2k})` reproduce the brief §1 expansion
exactly. See [task_1_baseline.md](task_1_baseline.md) for the full table.

---

## 3. Task 2 — Continued-fraction table for 8 candidates

All rows reached the target 300 partial quotients with two-precision
stability (2000 dps primary / 1500 dps cross-check / 2200 dps iteration),
except row (h) where precision is bounded by the partial-sum convergence
rate of the spectral series.

| row | candidate | value (first 30 sig fig) | CF length | n_fail |
|-----|-----------|--------------------------|-----------|--------|
| (a) | γ | `0.577215664901532860606512090082...` | 300 | 300 |
| (b) | 1/√3 | `0.577350269189625764509148780501...` | 300 | 300 |
| (c) | log φ | `0.481211825059603447497758913424...` | 300 | 300 |
| (d) | π/(3√3) = L(1, χ₋₃) | `0.604599788078072616864692752547...` | 300 | 300 |
| (e) | 2 log(φ)/√5 = L(1, χ₅) | `0.430408940964004038889433232950...` | 300 | 300 |
| (f) | γ₀(ℚ(√5)) | `0.604679472735033350795377487236...` | 300 | 300 |
| (g) | γ₀(ℚ(√−3)) | `0.571647606368562770069214424980...` | 300 | 300 |
| (h) | γ_{2I-spectral} | `-0.01796574694...` | 6 | 6 |

**Methodology for rows (f), (g):** Method B (Dedekind factorisation,
closed-form Hurwitz expansion of `L′(1, χ)`) for high precision;
Method A (direct ideal-counting partial sum + Richardson extrapolation)
as cross-check. Methods A and B agree to ~10⁻⁵–10⁻⁶ across **four**
real-quadratic fields (ℚ(√2), ℚ(√3), ℚ(√5), ℚ(√−3)) — see
[smoke_test_modules.py](smoke_test_modules.py) for the 4-field
validation.

**Methodology for row (h):** Molien-derived analytical slope `c = 1/120`
empirically confirmed to 9 decimal digits via partial sums up to N = 10⁷
(relative error 2.7 × 10⁻⁹). Convergence shape is `O(1/N)` rather than
the generic `N^{-1/2}` (the `1/√N` coefficient fits to ~10⁻⁷). γ_{2I-spectral}
extracted via 4-parameter Richardson fit at N ∈ {10⁴, 10⁵, 10⁶, 10⁷}.

**Reference reconciliation:** CinC initially gave 0.7615 as a Mathar-style
reference for γ₀(ℚ(√5)). Methods A and B both produce 0.6047. After
computing L(1, χ₁₂) = log(2+√3)/√3 = **0.7603**, the misfile became clear:
0.7615 was a misrecollection of the L-value for ℚ(√3) (different field,
different constant). CinC owned the error cleanly; Task 2 proceeded with
0.6047 as the working value.

Full data in [task_2_cf_data.json](task_2_cf_data.json); table in
[task_2_cf_table.md](task_2_cf_table.md).

### Observable CF patterns (formal evaluation in Task 3)

- Row (b) 1/√3 is exactly periodic `[0; 1, 1, 2, 1, 2, 1, 2, ...]` (period
  `[1, 2]`) — Lagrange's theorem on quadratic irrationals. Structural to
  ℚ(√3), not to the brief's question.
- Rows (a) and (b) share 7-term CF prefix `[0; 1, 1, 2, 1, 2, 1]`. This is
  a **numerical-closeness artifact**: |γ − 1/√3| ≈ 1.3 × 10⁻⁴.
- Rows (d) and (f) share 5-term CF prefix `[0; 1, 1, 1, 1]`. Also a
  closeness artifact: |π/(3√3) − γ₀(ℚ(√5))| ≈ **8 × 10⁻⁵**. An L-value
  and a Stieltjes constant from different fields landing 8 × 10⁻⁵ apart;
  under uniform-prior closest-pair across 28 pairs, p ≈ 0.08 — marginal.
- Extreme partial quotients (a_k > 100) occur in rows (d), (e), (f), (g),
  consistent with K-L heavy-tail expectations once multiple testing is
  controlled.

---

## 4. Task 3 — Relationship analysis

All four sub-tasks evaluated against Task 5's empirical null
(see §6). Bonferroni corrections applied within sub-tasks where
appropriate. None of the four sub-tasks yields an empirical p < 0.05
after multiple-testing correction.

### 3a Pairwise partial-quotient correlations (28 pairs)

K-L analytic null: nothing survives Bonferroni (0.05/28 ≈ 0.0018).
Two uncorrected hits: (a, c) Spearman p = 0.049, (c, h) Pearson p = 0.039
(n = 5, very low power).

Empirical-null re-evaluation: 1/21 pairs at p < 0.05 ((b, g) with p = 0.048),
matching the 1.05 expected false-positive rate at α = 0.05.

**Verdict: clean null.**

### 3b Within-row block-pattern outliers

Many z > 3 outliers under K-L analytic null per row (a) through (g),
including a striking-looking z = 16.05 on `(13, 5, 1, 1)` in row (a) γ.

Row (b) 1/√3 has z up to 457 on `(1, 2, 1, 2, 1, 2)` — Lagrange's theorem
on quadratic irrationals, not structural to the question; excluded from
primary evaluation.

Empirical null reveals K-L iid CFs of length 300 routinely produce
max-z = 28 (p = 0.05), max-z = 49 (p = 0.01) for block length 4, due
to multiple testing over ~10⁵ candidate block instances per CF. Row
(a)'s z = 16.05 has **empirical p = 0.154** — unremarkable.

**Verdict: 0 / 7 rows at empirical p < 0.05.**

### 3c Cross-CF block matching (length 4–10)

Most cross-CF pairs show *fewer* common length-4 blocks than the K-L iid
random pair median (typical: 50–80 obs vs 66 K-L-median, with 196 from
the K-L analytic formula that overestimates). Candidate CFs are
internally somewhat less "K-L-typical" than random.

Two excess outliers, both **numerical-closeness artifacts**:
- (a, b): 296 obs L = 4 common blocks (vs median 66), driven by the 7-term
  shared prefix from |γ − 1/√3| ≈ 1.3 × 10⁻⁴.
- (d, f): 52 obs L = 4 common blocks. Median, not excess — but their shared
  5-term CF prefix from |π/(3√3) − γ₀(ℚ(√5))| ≈ 8 × 10⁻⁵ is itself an
  artifact, even though the L = 4 block count isn't dramatic.

Both excluded from primary evaluation.

**Verdict: 0 / 26 non-closeness pairs at empirical p < 0.05.**

### 3d Convergent comparison

Universal small-q Fibonacci/Lucas hits at `q ∈ {1, 1, 2, 3, 5, 7}` are
trivial (Fibonacci is dense in small integers). Two single-observation
hits at "n = 10" positions worth flagging:
- Row (a) γ: q₇ = **123** = L₁₀ (Lucas).
- Row (h) γ_{2I-spectral}: q₂ = **55** = F₁₀ (Fibonacci).

Empirical null: random K-L CFs average 1.18 Fibonacci hits and 0.92 Lucas
hits in first 10 convergents, with p = 0.05 thresholds at 3 hits each.
Observed candidate Fibonacci/Lucas hit counts are all in the null bulk
(highest p = 0.097 for rows (d) and (f) at 3 Fibonacci hits each).

Icosahedral integer set `{12, 20, 30, 60}` produces zero hits across all
8 candidate rows, consistent with the empirical null mean of 0.06
icosahedral hits per random K-L CF.

**Verdict: 0 / 8 rows at empirical p < 0.05 for any of Fib / Lucas / Icos.**

---

## 5. Task 4 — Euler–Maclaurin prefix prediction test

**FAILS the brief §3 necessary condition** at the aggregate level
(p = 0.086 vs the required p < 0.05).

### Construction

`γ_K(n) = H_n − ln n − 1/(2n) + Σ_{k=1}^{K} B_{2k}/(2k·n^{2k})`,
truncation residual `R(n, K) ≈ |B_{2(K+1)}/(2(K+1)·n^{2(K+1)})|`.

For (n, K) ∈ {100, 10³, 10⁴, 10⁵, 10⁶} × {0, 1, 2, 5, 10} (25 cells),
compute γ_K(n) at 500 dps, CF it, record divergence position from CF(γ).
Null cohort: 500 K-L random g* per cell, perturbed by noise of magnitude
R(n, K), divergence positions recorded.

### Result

Per-cell true percentiles (γ's divergence position within null distribution)
average 59.8% (vs 50% expected under H0). Combined across 25 cells via
Stouffer's method:

- **Stouffer's z = 1.3659**, **one-sided p = 0.086**.
- **0 / 25 cells** at individual p < 0.05.
- **0 / 25 cells** at Bonferroni (α = 0.002).
- Mean true percentile = 59.8% (slight rightward bias, not significant).

There is a real slight bias — γ's CF agrees with itself a bit longer than
a random Stieltjes-like constant under EM-style noise — but the brief's
explicit necessary condition (p < 0.05 vs null) is not met. The "sufficient
condition" (p < 0.001, clear structural pattern, or two-or-more constants
showing related structure) is not met by any margin.

Re-analysis details and proper Stouffer methodology in
[task_4_reanalysis.md](task_4_reanalysis.md).

---

## 6. Task 5 — Null distribution

1000 K-L random CFs of length 300 (seed 20260520) + 1000 random pairs.
Empirical thresholds tabulated for each statistic:

| statistic | median | p = 0.05 | p = 0.01 | p = 0.001 |
|---|---|---|---|---|
| max-z within-CF (L = 4) | 8.23 | 28.25 | 48.54 | 253.30 |
| max-z within-CF (L = 5) | 8.21 | 34.51 | 62.95 | 393.85 |
| \|Pearson r\| pair | 0.016 | 0.079 | 0.432 | 0.924 |
| \|Spearman r\| pair | 0.039 | 0.116 | 0.153 | 0.209 |
| common L = 4 blocks per pair | 66 | 85 | 95 | 102 |
| common L = 6 blocks per pair | 6 | 15 | 21 | 24 |
| shared CF prefix length | 0 | 2 | 2 | 4 |
| Fibonacci hits in first 10 convergents | 1 | 3 | 5 | 9 |
| Lucas hits in first 10 convergents | 1 | 3 | 5 | 7 |
| Icosahedral {12, 20, 30, 60} hits | 0 | 1 | 1 | 1 |

These thresholds drive the empirical re-evaluation of Task 3 sub-tasks
and the null cohort for Task 4. Full distributions and re-evaluations in
[task_5_null_distribution.md](task_5_null_distribution.md) and
[task_5_null_distribution.json](task_5_null_distribution.json).

---

## 7. Stop condition outcome — **HALT-NULL**

Per brief §3 necessary condition for H1 survival:

> at least one of —
> - A statistically significant (p < 0.05 vs null distribution) CF-level
>   relationship between γ and one of {γ_{ℚ(√5)}, γ_{ℚ(√−3)}, γ_{2I-spectral}}, OR
> - Bernoulli-derived prefix prediction for γ-CF matching observed γ-CF
>   prefix at p < 0.05 vs null distribution of random Stieltjes-like constants

**Neither half of the necessary condition is met.**

- Task 3a–3d empirical null: clean across all four sub-tasks (0 / 21,
  0 / 7, 0 / 26, 0 / 8 at p < 0.05 with closeness artifacts and the
  ℚ(√3)-Lagrange-periodicity row excluded).
- Task 4 Stouffer combined p = 0.086. Above 0.05.

Per brief §5: HALT-NULL.

### Brief §5 prescription

> "**HALT-NULL:** all tests at p ≥ 0.05 against null. The broader
> trace-shadow framework does not specifically predict γ structure beyond
> chance. Paper 15 v3.3 §4 stands at OBSERVED structural permanent ceiling;
> Paper 190 remains v1.1 (φ paper); recommend a §5 sharpening in v1.1 to
> acknowledge the negative result without overclaiming the φ content."

### What this finding *does not* establish

- Paper 184's trace-shadow framework as a corpus-wide structural reading
  remains intact; what this Phase 1 tested was whether γ *specifically*
  carries Bernoulli-derived CF structure beyond K-L random.
- γ-irrationality and γ-CF leading-term predictability remain 300-year-old
  open problems. This test does not refute them; it merely shows the
  Bootstrap framework's machinery (Euler–Maclaurin + 2I-native ζ
  candidates) does not predict observable CF structure at the precisions
  and lengths examined here.

### What the close-pair (d, f) observation suggests for follow-up

The closest non-trivial close-pair finding is:

> |L(1, χ₋₃) − γ₀(ℚ(√5))| = |π/(3√3) − γ₀(ℚ(√5))| ≈ 8 × 10⁻⁵.

This is at p ≈ 0.08 under a uniform-prior closest-pair-across-28 model.
Not significant, but the two constants are structurally distinct (an
L-value of one field and a Stieltjes constant of a *different* field).
If a future Paper 190-or-derivative wanted a candidate to investigate
on independent grounds, this pair is the natural target. **Out of scope
for this Phase 1.**

---

## 8. Files

| file | content |
|---|---|
| [Mr_Code_Brief_Paper_190_Phase_1_v0_2.md](Mr_Code_Brief_Paper_190_Phase_1_v0_2.md) | the v0.2 brief, unchanged at pre-reg commit `971938b` |
| [task_1_bernoulli_baseline.py](task_1_bernoulli_baseline.py) + [task_1_baseline.md](task_1_baseline.md) | Task 1 |
| [cf_tools.py](cf_tools.py), [dedekind_stieltjes.py](dedekind_stieltjes.py), [spectral_2I.py](spectral_2I.py), [smoke_test_modules.py](smoke_test_modules.py) | shared modules + 4-field validation |
| [task_2_cf_table.py](task_2_cf_table.py) + [task_2_cf_table.md](task_2_cf_table.md) + [task_2_cf_data.json](task_2_cf_data.json) | Task 2 |
| [task_3_relationships.py](task_3_relationships.py) + [task_3_relationships.md](task_3_relationships.md) + [task_3_relationships.json](task_3_relationships.json) | Task 3 |
| [task_4_em_prefix_prediction.py](task_4_em_prefix_prediction.py) + [task_4_em_prefix_prediction.md](task_4_em_prefix_prediction.md) + [task_4_em_prefix_prediction.json](task_4_em_prefix_prediction.json) | Task 4 (original; methodology notes in commit) |
| [task_4_reanalysis.py](task_4_reanalysis.py) + [task_4_reanalysis.md](task_4_reanalysis.md) + [task_4_reanalysis.json](task_4_reanalysis.json) | Task 4 re-analysis (definitive) |
| [task_5_null_distribution.py](task_5_null_distribution.py) + [task_5_null_distribution.md](task_5_null_distribution.md) + [task_5_null_distribution.json](task_5_null_distribution.json) | Task 5 |

---

*Mr Code Phase 1 v0.2 reconnaissance, scope-limited per Pattern 97.
Pre-registration discipline maintained throughout; both halves of the
necessary condition for H1 survival failed; HALT-NULL.*

🐕☕⬡
