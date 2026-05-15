# Paper 200 v0.3 Task Z' - within-regime null test results

Pre-registered in `task_z2_brief.md`. Implementation in
`task_z2_within_regime_null.py`. JSON payload in `task_z2_results.json`.
Supersedes Task Z (commit `61aff82`, halted on wrong-convention anchor).

- seed = `200`, N (5-series primes) = `200`
- convention: Definition A (5-series = primes congruent to 5 mod 6;
  chi_5 Legendre mod 5; chi_5(5) = 0, ramified prime included)

## Anchor verification (row-by-row)

All 31 rows from CinC's hand-computation match the
locked-script output. `S_5(31) = 0` (expected 0) -> **OK**.

## Tie-position cross-check in n in [1, 33]

- Pre-registered expected set: [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33]
- Computed set              : [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33]
- Match: **True** (14 ties, matching Paper 199 sec.6.1)

## Regime fractions

| Regime | n range | zeros | length | fraction |
|---|---|---:|---:|---:|
| **tight** | [1, 33] | 14 | 33 | **0.4242** |
| post-regime | [34, 66] | 5 | 33 | 0.1515 |
| asymptotic | [67, 200] | 25 | 134 | 0.1866 |

## Decision against pre-registered rule

- LOW  : f_tight < 0.30
- MID  : 0.30 <= f_tight <= 0.55
- HIGH : f_tight > 0.55

f_tight = **0.4242** -> **Bin MID**

## Zero indices in [1, 200]

1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33, 41, 43, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 83, 85, 89, 101, 103, 105, 107, 109, 113, 115, 117, 121, 125, 131, 133, 137, 163, 165

## First 50 values of S_5(n) (5-series)

`*` marks rows where S_5(n) = 0.

| n | p_5n | p mod 5 | chi_5 | S_5(n) | tie |
|---:|---:|---:|---:|---:|:---:|
| 1 | 5 | 0 | +0 | +0 | * |
| 2 | 11 | 1 | +1 | +1 |  |
| 3 | 17 | 2 | -1 | +0 | * |
| 4 | 23 | 3 | -1 | -1 |  |
| 5 | 29 | 4 | +1 | +0 | * |
| 6 | 41 | 1 | +1 | +1 |  |
| 7 | 47 | 2 | -1 | +0 | * |
| 8 | 53 | 3 | -1 | -1 |  |
| 9 | 59 | 4 | +1 | +0 | * |
| 10 | 71 | 1 | +1 | +1 |  |
| 11 | 83 | 3 | -1 | +0 | * |
| 12 | 89 | 4 | +1 | +1 |  |
| 13 | 101 | 1 | +1 | +2 |  |
| 14 | 107 | 2 | -1 | +1 |  |
| 15 | 113 | 3 | -1 | +0 | * |
| 16 | 131 | 1 | +1 | +1 |  |
| 17 | 137 | 2 | -1 | +0 | * |
| 18 | 149 | 4 | +1 | +1 |  |
| 19 | 167 | 2 | -1 | +0 | * |
| 20 | 173 | 3 | -1 | -1 |  |
| 21 | 179 | 4 | +1 | +0 | * |
| 22 | 191 | 1 | +1 | +1 |  |
| 23 | 197 | 2 | -1 | +0 | * |
| 24 | 227 | 2 | -1 | -1 |  |
| 25 | 233 | 3 | -1 | -2 |  |
| 26 | 239 | 4 | +1 | -1 |  |
| 27 | 251 | 1 | +1 | +0 | * |
| 28 | 257 | 2 | -1 | -1 |  |
| 29 | 263 | 3 | -1 | -2 |  |
| 30 | 269 | 4 | +1 | -1 |  |
| 31 | 281 | 1 | +1 | +0 | * |
| 32 | 293 | 3 | -1 | -1 |  |
| 33 | 311 | 1 | +1 | +0 | * |
| 34 | 317 | 2 | -1 | -1 |  |
| 35 | 347 | 2 | -1 | -2 |  |
| 36 | 353 | 3 | -1 | -3 |  |
| 37 | 359 | 4 | +1 | -2 |  |
| 38 | 383 | 3 | -1 | -3 |  |
| 39 | 389 | 4 | +1 | -2 |  |
| 40 | 401 | 1 | +1 | -1 |  |
| 41 | 419 | 4 | +1 | +0 | * |
| 42 | 431 | 1 | +1 | +1 |  |
| 43 | 443 | 3 | -1 | +0 | * |
| 44 | 449 | 4 | +1 | +1 |  |
| 45 | 461 | 1 | +1 | +2 |  |
| 46 | 467 | 2 | -1 | +1 |  |
| 47 | 479 | 4 | +1 | +2 |  |
| 48 | 491 | 1 | +1 | +3 |  |
| 49 | 503 | 3 | -1 | +2 |  |
| 50 | 509 | 4 | +1 | +3 |  |
