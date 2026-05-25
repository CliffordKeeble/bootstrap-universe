# Task 1 — Bernoulli baseline verification

**Pre-reg commit:** `971938b`  
**n =** 10^6 &nbsp; **K =** 10 &nbsp; **working dps =** 60

## Computed value vs reference

```
gamma (computed)   = 0.5772156649015328606065120900824024310421593359399235988
gamma (mpmath ref) = 0.5772156649015328606065120900824024310421593359399235988
signed error       = 3.88938e-61
correct digits     ~ 60.4
next-term residual ~ 2.8146e-130
```

Note: observed error ~4×10<sup>−61</sup> reflects working-precision rounding accumulation in the H_n sum (60 dps), not the EM truncation bound at K=10 (~3×10<sup>−130</sup>). Verification passes by a wide margin against the ≥50-digit target.

## First 10 asymptotic terms of H_n − ln n − γ

| order | B<sub>2k</sub> (exact) | signed coefficient −B<sub>2k</sub>/(2k) | numeric @ n=10<sup>6</sup> |
|-------|------------------------|------------------------------------------|----------------------------|
| n<sup>−1</sup> | (half-term) | +1/2 | 5.0e-7 |
| n<sup>−2</sup> | 1/6 | -1/12 | -8.333333333333333e-14 |
| n<sup>−4</sup> | -1/30 | 1/120 | 8.333333333333333e-27 |
| n<sup>−6</sup> | 1/42 | -1/252 | -3.968253968253968e-39 |
| n<sup>−8</sup> | -1/30 | 1/240 | 4.166666666666667e-51 |
| n<sup>−10</sup> | 5/66 | -1/132 | -7.575757575757576e-63 |
| n<sup>−12</sup> | -691/2730 | 691/32760 | 2.109279609279609e-74 |
| n<sup>−14</sup> | 7/6 | -1/12 | -8.333333333333333e-86 |
| n<sup>−16</sup> | -3617/510 | 3617/8160 | 4.432598039215686e-97 |
| n<sup>−18</sup> | 43867/798 | -43867/14364 | -3.05395433027012e-108 |
| n<sup>−20</sup> | -174611/330 | 174611/6600 | 2.645621212121212e-119 |

## Verdict

**PASS.** Computed γ matches reference to ~60 decimal digits (target: ≥ 50). Euler–Maclaurin / Bernoulli machinery verified. Cleared to proceed to Task 2.
