# Paper 200 v0.5 Task AA - Def A/B/C verification results

Pre-registered in `task_aa_brief.md`. Implementation in
`task_aa_definitions_verification.py`. JSON payload in
`task_aa_results.json`.

Authoritative source for Def A/B/C: Paper 199 sec.3,
lines 122-128 and 54 of `Paper_199_RunningSumOf5_Series Primes_v1_2.md`.

## Anchor cross-checks against Task Z' (commit c2c5eeb)

| Check | Expected | Computed | Status |
|---|---|---|---|
| p_{5,1} | 5 | 5 | OK |
| p_{5,5} | 29 | 29 | OK |
| p_{5,31} | 281 | 281 | OK |
| p_{5,33} | 311 | 311 | OK |
| S_A(31) | 0 | 0 | OK |
| Def A tie set in [1, 33] | 14 positions | 14 positions | OK |

## Tie counts under each definition (Paper 199 sec.3)

| Definition | Window | Tie count | Paper 199 sec.3 figure |
|---|---|---:|---:|
| **Def A** (include n=1) | [1, 33] | 14 | 14 |
| **Def B** (exclude n=1; same sum, window restricted) | [2, 33] | 13 | 13 |
| **Def C** (omit p=5; 32-step sum) | [1, 32] | 13 | 13 |

## Running sum at the boundary-crossing prime p = 281

| Definition | Index at p = 281 | S_5 at that index |
|---|---:|---:|
| **Def A** | k = 31 | **+0** |
| **Def B** | k = 31 (same indexing as A) | **+0** |
| **Def C** | j = 30 (after skipping ramified) | **+0** |

All three definitions give the same value (0) at the boundary-crossing
prime, confirming Paper 199 sec.3's robustness claim.

## CinC literal-prediction scoring

CinC's prediction table (brief) listed values for specific
(definition, index) pairs. Scoring literally:

| Prediction | CinC value | Computed | Match |
|---|---:|---:|:---:|
| S_A(31) | 0 | 0 | YES |
| S_B(30) under Paper 199 Def B | 0 | -1 | NO |
| S_C at boundary p=281 (j=30 in Def C) | 0 | +0 | YES |
| Def A 14-tie set in [1, 33] | 14 | 14 | YES |
| Def B 13 ties under 'shift by 1' reading | 13 | 13 | YES |

## Tie-position sets under each definition

- **Def A** [1, 33], 14 ties: [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33]
- **Def B** [2, 33], 13 ties: [3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33]
- **Def C** [1, 32], 13 ties: [2, 4, 6, 8, 10, 14, 16, 18, 20, 22, 26, 30, 32]

Note Def A tie set = Def B tie set + {1}; Def C tie set is the
Def A tie set in [2, 33], reindexed by subtracting 1.

## Decision against pre-registered rule

**Bin MISMATCH-MINOR.**

Substantive observation holds under all three definitions, but Paper 199 sec.3's Def B does not shift indexing as CinC's working reading assumes (Paper 199 Def B uses the same running sum as Def A with a restricted counting window; Paper 199 Def C is the definition that shifts indexing). CinC's predicted value for S_B(30) is the value Paper 199 Def C gives at j=30, not what Def B gives at k=30.

## Side-by-side S_5 values for n in [1, 50]

Indexing convention: n indexes Def A / Def B (same sum); the
`S_C @ j=n-1` column shows Def C's running sum at its index j = n-1,
which corresponds to the same underlying prime p_{5,n} as Def A's index n.
(`-` for n=1 since the ramified prime has no Def C index.)

| n | p_5n | p mod 5 | chi_5 | S_A=S_B | S_C @ j=n-1 | A tie | B tie | C tie |
|---:|---:|---:|---:|---:|---:|:-:|:-:|:-:|
| 1 | 5 | 0 | +0 | +0 | - | * |  |  |
| 2 | 11 | 1 | +1 | +1 | +1 |  |  |  |
| 3 | 17 | 2 | -1 | +0 | +0 | * | * | * |
| 4 | 23 | 3 | -1 | -1 | -1 |  |  |  |
| 5 | 29 | 4 | +1 | +0 | +0 | * | * | * |
| 6 | 41 | 1 | +1 | +1 | +1 |  |  |  |
| 7 | 47 | 2 | -1 | +0 | +0 | * | * | * |
| 8 | 53 | 3 | -1 | -1 | -1 |  |  |  |
| 9 | 59 | 4 | +1 | +0 | +0 | * | * | * |
| 10 | 71 | 1 | +1 | +1 | +1 |  |  |  |
| 11 | 83 | 3 | -1 | +0 | +0 | * | * | * |
| 12 | 89 | 4 | +1 | +1 | +1 |  |  |  |
| 13 | 101 | 1 | +1 | +2 | +2 |  |  |  |
| 14 | 107 | 2 | -1 | +1 | +1 |  |  |  |
| 15 | 113 | 3 | -1 | +0 | +0 | * | * | * |
| 16 | 131 | 1 | +1 | +1 | +1 |  |  |  |
| 17 | 137 | 2 | -1 | +0 | +0 | * | * | * |
| 18 | 149 | 4 | +1 | +1 | +1 |  |  |  |
| 19 | 167 | 2 | -1 | +0 | +0 | * | * | * |
| 20 | 173 | 3 | -1 | -1 | -1 |  |  |  |
| 21 | 179 | 4 | +1 | +0 | +0 | * | * | * |
| 22 | 191 | 1 | +1 | +1 | +1 |  |  |  |
| 23 | 197 | 2 | -1 | +0 | +0 | * | * | * |
| 24 | 227 | 2 | -1 | -1 | -1 |  |  |  |
| 25 | 233 | 3 | -1 | -2 | -2 |  |  |  |
| 26 | 239 | 4 | +1 | -1 | -1 |  |  |  |
| 27 | 251 | 1 | +1 | +0 | +0 | * | * | * |
| 28 | 257 | 2 | -1 | -1 | -1 |  |  |  |
| 29 | 263 | 3 | -1 | -2 | -2 |  |  |  |
| 30 | 269 | 4 | +1 | -1 | -1 |  |  |  |
| 31 | 281 | 1 | +1 | +0 | +0 | * | * | * |
| 32 | 293 | 3 | -1 | -1 | -1 |  |  |  |
| 33 | 311 | 1 | +1 | +0 | +0 | * | * | * |
| 34 | 317 | 2 | -1 | -1 | -1 |  |  |  |
| 35 | 347 | 2 | -1 | -2 | -2 |  |  |  |
| 36 | 353 | 3 | -1 | -3 | -3 |  |  |  |
| 37 | 359 | 4 | +1 | -2 | -2 |  |  |  |
| 38 | 383 | 3 | -1 | -3 | -3 |  |  |  |
| 39 | 389 | 4 | +1 | -2 | -2 |  |  |  |
| 40 | 401 | 1 | +1 | -1 | -1 |  |  |  |
| 41 | 419 | 4 | +1 | +0 | +0 | * | * | * |
| 42 | 431 | 1 | +1 | +1 | +1 |  |  |  |
| 43 | 443 | 3 | -1 | +0 | +0 | * | * | * |
| 44 | 449 | 4 | +1 | +1 | +1 |  |  |  |
| 45 | 461 | 1 | +1 | +2 | +2 |  |  |  |
| 46 | 467 | 2 | -1 | +1 | +1 |  |  |  |
| 47 | 479 | 4 | +1 | +2 | +2 |  |  |  |
| 48 | 491 | 1 | +1 | +3 | +3 |  |  |  |
| 49 | 503 | 3 | -1 | +2 | +2 |  |  |  |
| 50 | 509 | 4 | +1 | +3 | +3 |  |  |  |
