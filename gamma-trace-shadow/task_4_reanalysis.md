# Task 4 re-analysis

**Pre-reg commit:** `971938b`

Re-analysis of Task 4 using proper percentile computation against raw null arrays (N_NULL = 500 K-L random g* per cell), and Stouffer's combined z-test for an overall effect across the 25 (n, K) cells.

Corrected K-L prediction: $\mathrm{pos} \sim -\ln(R)/(2L)$ with Lévy's constant $L = \pi^2/(12 \ln 2) \approx 1.1866$ (natural-log convention).

## Per-cell results (true percentile vs raw null array)

| n | K | R | K-L pred | obs div | null median | null p90 | true %ile | one-sided p |
|---|---|---|----------|---------|-------------|----------|-----------|-------------|
| 100 | 0 | 8.33e-06 | 4.9 | 8 | 6.0 | 8.0 | 81.9% | 0.1806 |
| 100 | 1 | 8.33e-11 | 9.8 | 13 | 11.0 | 14.0 | 75.3% | 0.2465 |
| 100 | 2 | 3.97e-15 | 14.0 | 17 | 15.0 | 19.0 | 71.8% | 0.2824 |
| 100 | 5 | 2.11e-26 | 24.9 | 26 | 26.0 | 31.0 | 48.4% | 0.5160 |
| 100 | 10 | 2.81e-42 | 40.3 | 37 | 41.0 | 48.0 | 20.5% | 0.7954 |
| 1000 | 0 | 8.33e-08 | 6.9 | 9 | 8.0 | 11.0 | 64.4% | 0.3563 |
| 1000 | 1 | 8.33e-15 | 13.7 | 16 | 15.0 | 19.0 | 62.8% | 0.3723 |
| 1000 | 2 | 3.97e-21 | 19.8 | 22 | 20.0 | 25.0 | 61.2% | 0.3882 |
| 1000 | 5 | 2.11e-38 | 36.6 | 36 | 38.0 | 44.0 | 36.8% | 0.6317 |
| 1000 | 10 | 2.81e-64 | 61.7 | 65 | 62.0 | 72.0 | 60.2% | 0.3982 |
| 10000 | 0 | 8.33e-10 | 8.8 | 10 | 9.0 | 13.0 | 50.2% | 0.4980 |
| 10000 | 1 | 8.33e-19 | 17.5 | 21 | 18.0 | 23.0 | 71.2% | 0.2884 |
| 10000 | 2 | 3.97e-27 | 25.6 | 28 | 26.0 | 32.0 | 61.0% | 0.3902 |
| 10000 | 5 | 2.11e-50 | 48.2 | 45 | 49.0 | 57.0 | 23.7% | 0.7635 |
| 10000 | 10 | 2.81e-86 | 83.0 | 89 | 84.0 | 94.0 | 70.6% | 0.2944 |
| 100000 | 0 | 8.33e-12 | 10.7 | 13 | 11.0 | 15.0 | 64.8% | 0.3523 |
| 100000 | 1 | 8.33e-23 | 21.4 | 23 | 22.0 | 27.0 | 55.4% | 0.4461 |
| 100000 | 2 | 3.97e-33 | 31.4 | 31 | 32.0 | 39.0 | 35.4% | 0.6457 |
| 100000 | 5 | 2.11e-62 | 59.8 | 62 | 60.0 | 69.0 | 56.6% | 0.4341 |
| 100000 | 10 | 2.81e-108 | 104.4 | 116 | 106.0 | 116.0 | 89.7% | 0.1028 |
| 1000000 | 0 | 8.33e-14 | 12.7 | 16 | 13.0 | 18.0 | 75.7% | 0.2425 |
| 1000000 | 1 | 8.33e-27 | 25.3 | 28 | 26.0 | 31.1 | 62.8% | 0.3723 |
| 1000000 | 2 | 3.97e-39 | 37.3 | 36 | 38.0 | 45.0 | 31.6% | 0.6836 |
| 1000000 | 5 | 2.11e-74 | 71.5 | 78 | 72.0 | 81.0 | 76.5% | 0.2345 |
| 1000000 | 10 | 2.81e-130 | 125.7 | 138 | 126.0 | 139.0 | 86.7% | 0.1327 |

## Combined analysis

- **Cells**: 25
- **Stouffer's combined z (one-sided)**: 1.3659
- **Stouffer's combined p (one-sided)**: **0.0860**
- **Cells at p<0.05 individually**: 0/25 (expected under H0: ~1.25)
- **Cells at p<0.002 (Bonferroni)**: 0/25
- **Mean true percentile**: 59.80% (H0 expected: 50.0%)

## Interpretation

Combined Stouffer's p = 0.0860 is above 0.05; the brief's 'necessary condition' for H1 survival is not met. No individual cell reaches p < 0.05. No cell survives Bonferroni correction.

If Stouffer's combined p < 0.05 but no individual cell survives Bonferroni, this is a marginal aggregate effect — γ's CF tracks itself slightly longer than a random Stieltjes-like constant approximated at comparable precision, but the bias is small and the brief's 'sufficient condition' (clear structural pattern, two-or-more constants at p < 0.001, or clean discriminating prediction) is not met.
