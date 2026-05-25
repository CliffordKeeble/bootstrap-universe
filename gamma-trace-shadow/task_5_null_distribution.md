# Task 5 - Null distribution (K-L random Stieltjes-like constants)

**Pre-reg commit:** `971938b`  
**Seed:** 20260520  
**Samples:** 1000 CFs of length 300, 1000 random pairs.

## Empirical thresholds

### Max-z within-CF (vs K-L analytic) by block length

| L | median | p=0.05 | p=0.01 | p=0.001 | max |
|---|--------|--------|--------|---------|-----|
| 2 | 5.64 | 15.63 | 27.39 | 67.01 | 67.01 |
| 3 | 7.39 | 22.26 | 47.73 | 162.90 | 162.90 |
| 4 | 8.23 | 28.25 | 48.54 | 253.30 | 253.30 |
| 5 | 8.21 | 34.51 | 62.95 | 393.85 | 393.85 |
| 6 | 6.46 | 30.05 | 69.71 | 229.67 | 229.67 |
| 7 | 0.00 | 23.02 | 53.01 | 186.30 | 186.30 |
| 8 | 0.00 | 16.21 | 49.56 | 191.71 | 191.71 |
| 9 | 0.00 | 5.80 | 31.51 | 298.09 | 298.09 |
| 10 | 0.00 | 0.00 | 20.00 | 83.37 | 83.37 |

### |Pearson r| and |Spearman r| (pairwise, n=300)

| stat | median | p=0.05 | p=0.01 | p=0.001 | max |
|------|--------|--------|--------|---------|-----|
| |Pearson r| | 0.0159 | 0.0785 | 0.4320 | 0.9243 | 0.9243 |
| |Spearman r| | 0.0387 | 0.1158 | 0.1527 | 0.2090 | 0.2090 |

### Count of common blocks across a random pair, by length

| L | median | mean | p=0.05 | p=0.01 | p=0.001 | max |
|---|--------|------|--------|--------|---------|-----|
| 4 | 66 | 65.6 | 85 | 95 | 102 | 102 |
| 5 | 24 | 24.0 | 38 | 46 | 56 | 56 |
| 6 | 6 | 7.0 | 15 | 21 | 24 | 24 |
| 7 | 1 | 1.7 | 5 | 9 | 13 | 13 |
| 8 | 0 | 0.4 | 2 | 4 | 7 | 7 |
| 9 | 0 | 0.1 | 1 | 2 | 4 | 4 |
| 10 | 0 | 0.0 | 0 | 1 | 3 | 3 |

### Shared CF prefix length (pair, K-L iid)

| median | mean | p=0.05 | p=0.01 | p=0.001 | max |
|--------|------|--------|--------|---------|-----|
| 0 | 0.29 | 2 | 2 | 4 | 4 |

### Fibonacci / Lucas / Icos hits in first 10 convergents

| set | mean | p=0.05 | p=0.01 | p=0.001 | max |
|-----|------|--------|--------|---------|-----|
| Fibonacci | 1.18 | 3 | 5 | 9 | 9 |
| Lucas | 0.92 | 3 | 5 | 7 | 7 |
| Icosahedral {12,20,30,60} | 0.06 | 1 | 1 | 1 | 1 |

## Re-evaluation of Task 3 flagged findings vs empirical null

### Task 3a Pearson r (re-evaluated)

| pair | observed r | empirical p (|r|) | flag |
|------|------------|-------------------|------|
| a-b | +0.0116 | 0.6620 | ns |
| a-c | -0.0256 | 0.2560 | ns |
| a-d | -0.0177 | 0.4500 | ns |
| a-e | +0.0150 | 0.5370 | ns |
| a-f | -0.0394 | 0.1130 | ns |
| a-g | -0.0156 | 0.5170 | ns |
| b-c | -0.0457 | 0.0940 | ns |
| b-d | -0.0192 | 0.4050 | ns |
| b-e | +0.0153 | 0.5250 | ns |
| b-f | -0.0242 | 0.2880 | ns |
| b-g | +0.0809 | 0.0480 | **SIGNIF** (p<0.05) |
| c-d | -0.0200 | 0.3840 | ns |
| c-e | -0.0189 | 0.4120 | ns |
| c-f | -0.0231 | 0.3090 | ns |
| c-g | -0.0268 | 0.2410 | ns |
| d-e | -0.0219 | 0.3430 | ns |
| d-f | -0.0210 | 0.3530 | ns |
| d-g | -0.0139 | 0.5750 | ns |
| e-f | -0.0240 | 0.2920 | ns |
| e-g | +0.0252 | 0.2660 | ns |
| f-g | +0.0040 | 0.9230 | ns |

### Task 3b max z (length 4, per row) (re-evaluated)

| row | observed max-z (L=4) | empirical p | flag |
|-----|----------------------|-------------|------|
| a | 16.05 | 0.1540 | ns |
| b | 120.76 | n/a | excluded (1/sqrt 3 quadratic-irrational periodicity) |
| c | 4.06 | 0.9340 | ns |
| d | 5.76 | 0.7770 | ns |
| e | 7.07 | 0.6240 | ns |
| f | 5.38 | 0.8200 | ns |
| g | 5.35 | 0.8310 | ns |
| h | 0.00 | 1.0000 | ns |

### Task 3c common length-4 blocks per pair (re-evaluated)

Note: "close pairs" (a,b) and (d,f) are excluded from primary evaluation because their common blocks are driven by trivial numerical closeness; closeness artifact is *measured* by shared-prefix length in the next column.

| pair | obs common L=4 | empirical p (vs K-L iid pair) | flag |
|------|----------------|-------------------------------|------|
| a-b | 296 | 0.0000 | close-pair artifact |
| a-c | 67 | 0.4660 | ns |
| a-d | 38 | 0.9910 | ns |
| a-e | 75 | 0.2120 | ns |
| a-f | 56 | 0.8060 | ns |
| a-g | 55 | 0.8250 | ns |
| a-h | 0 | 1.0000 | ns |
| b-c | 5 | 1.0000 | ns |
| b-d | 3 | 1.0000 | ns |
| b-e | 15 | 1.0000 | ns |
| b-f | 5 | 1.0000 | ns |
| b-g | 7 | 1.0000 | ns |
| b-h | 0 | 1.0000 | ns |
| c-d | 45 | 0.9650 | ns |
| c-e | 69 | 0.3950 | ns |
| c-f | 60 | 0.6800 | ns |
| c-g | 83 | 0.0730 | ns |
| c-h | 2 | 1.0000 | ns |
| d-e | 52 | 0.8940 | ns |
| d-f | 52 | 0.8940 | close-pair artifact |
| d-g | 73 | 0.2740 | ns |
| d-h | 0 | 1.0000 | ns |
| e-f | 46 | 0.9590 | ns |
| e-g | 65 | 0.5340 | ns |
| e-h | 0 | 1.0000 | ns |
| f-g | 82 | 0.0890 | ns |
| f-h | 0 | 1.0000 | ns |
| g-h | 0 | 1.0000 | ns |

### Task 3d Fibonacci / Lucas / Icos hits per row

| row | Fib obs | Fib p | Lucas obs | Lucas p | Icos obs | Icos p |
|-----|---------|-------|-----------|---------|----------|--------|
| a | 2 | 0.2770 | 2 | 0.1740 | 0 | 1.0000 |
| b | 2 | 0.2770 | 1 | 0.6620 | 0 | 1.0000 |
| c | 1 | 0.7440 | 0 | 1.0000 | 0 | 1.0000 |
| d | 3 | 0.0970 | 1 | 0.6620 | 0 | 1.0000 |
| e | 1 | 0.7440 | 1 | 0.6620 | 0 | 1.0000 |
| f | 3 | 0.0970 | 1 | 0.6620 | 0 | 1.0000 |
| g | 2 | 0.2770 | 1 | 0.6620 | 0 | 1.0000 |
| h | 1 | 0.7440 | 0 | 1.0000 | 0 | 1.0000 |

## Verdict for Phase 1 stop-conditions

- Task 3a empirical-null SIGNIF (p<0.05): **1 / 21**
- Task 3b empirical-null SIGNIF (p<0.05, excl. row b periodicity): **0 / 7**
- Task 3c empirical-null SIGNIF (p<0.05, excl. close pairs): **0 / 26**
- Task 3d Fib/Lucas/Icos SIGNIF (any of 3 stats per row, p<0.05): **0 / 8**
