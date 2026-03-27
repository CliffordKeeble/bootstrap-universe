# Knuth Circle Product (Fibonacci Multiplication)

Implementation and analysis of the Knuth circle product ⊚ — a commutative, associative multiplication on positive integers defined via Zeckendorf (Fibonacci) representations.

**Reference:** D. E. Knuth, "Fibonacci Multiplication", *Applied Mathematics Letters*, Vol. 1, No. 1, pp. 57–60, 1988.

Part of the [Bootstrap Universe Programme](https://zenodo.org/communities/bootstrap-universe/) — Paper 168 computational appendix.

## The Product

Every positive integer has a unique Zeckendorf representation as a sum of non-consecutive Fibonacci numbers. The circle product treats these representations as polynomials in φ (the golden ratio) and multiplies them:

```
  2 = F₃        →  φ¹
  3 = F₄        →  φ²
  2 ⊚ 2 = 3     (φ¹ · φ¹ = φ² → F₄ = 3)
  2 ⊚ 3 = 5     (φ¹ · φ² = φ³ → F₅ = 5)
  3 ⊚ 3 = 8     (φ² · φ² = φ⁴ → F₆ = 8)
```

## The Implementation Bug (and Fix)

**Nine versions were required to achieve a correct implementation.** The difficulty is in normalisation.

When multiplying multi-digit numbers, the convolution can produce counts ≥ 2 at a single position. The carry rule is:

```
  2·φᵖ  →  φᵖ⁺¹ + φᵖ⁻²
```

For small p, this pushes digits to **negative** base-φ positions (p−2 < 0). These positions have real φ-values:

```
  φ⁻¹ = φ − 1 ≈ 0.618
  φ⁻² = 2 − φ ≈ 0.382
  φ⁻³ = 2φ − 3 ≈ 0.236
```

**The bug in versions 1–8:** Negative positions were dropped (because F₀ = 0) or merged (because F₁ = F₂ = 1). This corrupts the representation, changing subsequent products and breaking associativity (~12% failure rate on random triples).

**The fix in version 9:** Track ALL base-φ positions including negative ones throughout computation. Only convert to integers at final output. Associativity: 0 failures in 5,000 random triples.

## Results

### Product Table (2–10)

```
 ⊚  |     2     3     4     5     6     7     8     9    10
----+-------------------------------------------------------
  2 |     3     5     7     8    10    11    13    15    16
  3 |     5     8    11    13    16    18    21    24    26
  4 |     7    11    15    18    22    25    29    33    36
  5 |     8    13    18    21    26    29    34    39    42
  6 |    10    16    22    26    32    36    42    48    52
  7 |    11    18    25    29    36    40    47    54    58
  8 |    13    21    29    34    42    47    55    63    68
  9 |    15    24    33    39    48    54    63    72    78
 10 |    16    26    36    42    52    58    68    78    84
```

### Algebraic Properties

| Property | Status |
|---|---|
| Identity (1 ⊚ n = n) | ✓ Verified, n ∈ [1, 500] |
| Commutativity | ✓ Verified, all pairs in [2, 100] |
| Associativity | ✓ Verified, 5000 random triples in [2, 100] |
| Distributivity over + | ✗ Fails (~25% of triples) |

The circle product forms a **commutative monoid** on ℤ⁺, but NOT a ring (no distributivity).

### Circle Primes

A *circle prime* is an integer n > 1 that cannot be written as a ⊚ b for any 1 < a, b < n.

Up to 500:
- **95** ordinary primes
- **88** circle primes
- **24** integers that are BOTH ordinary and circle prime
- **71** ordinary primes that are ⊚-composite (e.g. 3 = 2⊚2, 5 = 2⊚3, 7 = 2⊚4)
- **64** ordinary composites that are ⊚-irreducible (e.g. 4, 6, 9, 12, 14, 27, 30, 35)

### Unique Factorisation

**Unique factorisation under ⊚ FAILS.** Testing integers 2–200:
- 150 have unique circle-prime decomposition
- 34 have non-unique decomposition (e.g. 15 = 4⊚4 = 2⊚9)
- 15 have no circle-prime decomposition at all

### Distribution

Circle primes distribute nearly uniformly mod 5 (15, 18, 18, 18, 19), unlike ordinary primes which avoid 0 mod 5. This reflects the ramification of 5 in ℤ[φ].

87 of 88 circle primes up to 500 contain F₂ in their Zeckendorf representation. The sole exception is 2 = F₃.

## Usage

```python
from knuth_circle_product import knuth_product, find_circle_primes

# Basic multiplication
print(knuth_product(3, 5))   # 13
print(knuth_product(8, 8))   # 55

# Find circle primes up to 100
cprimes = find_circle_primes(100)
print(cprimes)  # [2, 4, 6, 9, 12, 14, 17, 19, 27, 30, 35, 38, ...]

# Full analysis
python3 knuth_circle_product.py
```

### Requirements

Python 3.8+ (standard library only, no dependencies).

## Significance for Paper 168

The failure of unique factorisation under ⊚ demonstrates that the Fundamental Theorem of Arithmetic (FTA) and Zeckendorf's Theorem — dual projections of ℤ[φ] — genuinely differ in algebraic depth. The two notions of "prime" are almost disjoint: what ordinary arithmetic calls irreducible, Fibonacci arithmetic can decompose, and vice versa.

The duality is geometric (both structures project from the same ring) but not algebraic (they do not share the same factorisation theory).

## Author

Dr. Clifford Keeble  
ORCID: [0009-0003-6828-2155](https://orcid.org/0009-0003-6828-2155)  
Bootstrap Universe Programme, Woodbridge, UK

## Licence

MIT — see [LICENCE](LICENCE).
