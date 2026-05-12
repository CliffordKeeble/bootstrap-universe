# Task A — Null test for rail coverage

Paper 119 v1.2 computational task A. Pre-registration: `task_brief.md`.

- Seed: `119`
- Trials: `100,000`
- Realised data: `{1, 2, 3, 4, 5, 6, 7, 8, 13}`
- Rails (Fibonacci ∪ Lucas ∪ {6}): `{1, 2, 3, 4, 5, 6, 7, 8, 11, 13}` (size 10)
- Discriminating range [9, 13]: `{9, 10, 11, 12, 13}`

## Sub-table 1 — Coverage by sequence

| Sequence | Size | Cov [1,13] | Cov [1,8] | Cov [9,13] |
|---|---|---|---|---|
| Rails | 10 | 9/9 = 1.00 | 8/8 = 1.00 | 1/1 = 1.00 |
| Squares | 3 | 2/9 = 0.22 | 2/8 = 0.25 | 0/1 = 0.00 |
| Triangulars | 4 | 3/9 = 0.33 | 3/8 = 0.38 | 0/1 = 0.00 |
| Primes | 6 | 5/9 = 0.56 | 4/8 = 0.50 | 1/1 = 1.00 |
| Combined | 11 | 8/9 = 0.89 | 7/8 = 0.88 | 1/1 = 1.00 |

## Sub-table 2 — Precision in [9, 13]

| Sequence | Predictions in [9,13] | Hits | Precision |
|---|---|---|---|
| Rails | {11, 13} | {13} | 1/2 = 0.50 |
| Squares | {9} | {} | 0/1 = 0.00 |
| Triangulars | {10} | {} | 0/1 = 0.00 |
| Primes | {11, 13} | {13} | 1/2 = 0.50 |
| Combined | {9, 10, 11, 13} | {13} | 1/4 = 0.25 |

## Sub-table 3 — Random null distribution (100,000 trials, size-10 subsets of [1,13])

| Statistic | Coverage [1,13] | Coverage [9,13] |
|---|---|---|
| Rail value | 1.0000 | 1.0000 |
| Mean | 0.7690 | 0.7680 |
| Std dev | 0.0809 | 0.4221 |
| Rail percentile (≤ rail) | 100.00% | 100.00% |
| Empirical P(random ≥ rail) | 0.0138 | 0.7680 |

## Status flags

- Realised data and rail structure: DERIVED (Paper 119 §5 branching).
- Choice of comparison sequences: pre-registered (squares, triangulars, primes, random size-10).
- Coverage rates and p-values: OBSERVED.
- Interpretation: OBSERVED-with-empirical-support; see findings.md for the decision against pre-registered criterion.
