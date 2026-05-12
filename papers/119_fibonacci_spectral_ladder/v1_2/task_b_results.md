# Task B — Spectral gap of S^3/2T and S^3/2O

Paper 119 v1.2 computational task B. Pre-registration: `task_brief.md`.

Frobenius reciprocity computation: for each binary polyhedral group Γ,
the smallest l > 0 such that the trivial Γ-representation appears in
D_{l/2}|_Γ. Then λ₁ = l(l+2). Class data verified against Coxeter /
Cisneros-Molina 2000 (class sizes sum to |Γ|; multiplicities are exact
non-negative integers throughout).

## Result table

| Group | Order | l (smallest l > 0 with trivial-rep) | λ₁ = l(l+2) | On Fibonacci rail? | On Lucas rail? | On face-closure rail? | Off all rails? |
|---|---|---|---|---|---|---|---|
| 2T | 24 | 6 | 48 | — | L_4^2 - 1 | — | no |
| 2O | 48 | 8 | 80 | — | — | — | yes |
| 2I (verification ✓) | 120 | 12 | 168 | F_7^2 - 1 | — | — | no |

## Multiplicities mult_trivial(l) for l = 0, …, 20

Useful audit: shows which D_{l/2} contain the trivial sub at all, not just the first one.

| l | dim D_{l/2} = l+1 | λ = l(l+2) | mult(2T) | mult(2O) | mult(2I) |
|---|---|---|---|---|---|
| 0 | 1 | 0 | 1 | 1 | 1 |
| 1 | 2 | 3 | 0 | 0 | 0 |
| 2 | 3 | 8 | 0 | 0 | 0 |
| 3 | 4 | 15 | 0 | 0 | 0 |
| 4 | 5 | 24 | 0 | 0 | 0 |
| 5 | 6 | 35 | 0 | 0 | 0 |
| 6 | 7 | 48 | 1 | 0 | 0 |
| 7 | 8 | 63 | 0 | 0 | 0 |
| 8 | 9 | 80 | 1 | 1 | 0 |
| 9 | 10 | 99 | 0 | 0 | 0 |
| 10 | 11 | 120 | 0 | 0 | 0 |
| 11 | 12 | 143 | 0 | 0 | 0 |
| 12 | 13 | 168 | 2 | 1 | 1 |
| 13 | 14 | 195 | 0 | 0 | 0 |
| 14 | 15 | 224 | 1 | 0 | 0 |
| 15 | 16 | 255 | 0 | 0 | 0 |
| 16 | 17 | 288 | 1 | 1 | 0 |
| 17 | 18 | 323 | 0 | 0 | 0 |
| 18 | 19 | 360 | 2 | 1 | 0 |
| 19 | 20 | 399 | 0 | 0 | 0 |
| 20 | 21 | 440 | 2 | 1 | 1 |

## Conjugacy class data used

### 2T (order 24)

| Class | Size | Half-angle α |
|---|---|---|
| {1} | 1 | 0 |
| {-1} | 1 | pi |
| C3 | 4 | 2*pi/3 |
| C3' | 4 | 4*pi/3 |
| C6 | 4 | pi/3 |
| C6' | 4 | 5*pi/3 |
| C4 | 6 | pi/2 |

### 2O (order 48)

| Class | Size | Half-angle α |
|---|---|---|
| {1} | 1 | 0 |
| {-1} | 1 | pi |
| C4 (face) | 6 | pi/2 |
| C3 | 8 | 2*pi/3 |
| C6 | 8 | pi/3 |
| C8 | 6 | pi/4 |
| C8' | 6 | 3*pi/4 |
| C2 (edge) | 12 | pi/2 |

### 2I (order 120)

| Class | Size | Half-angle α |
|---|---|---|
| {1} | 1 | 0 |
| {-1} | 1 | pi |
| C10a | 12 | pi/5 |
| C5a | 12 | 2*pi/5 |
| C10b | 12 | 3*pi/5 |
| C5b | 12 | 4*pi/5 |
| C6 | 20 | pi/3 |
| C3 | 20 | 2*pi/3 |
| C4 | 30 | pi/2 |

## Status flags

- Frobenius formula: DERIVED.
- Conjugacy class data: verified-against-reference (Coxeter, Cisneros-Molina 2000); class sizes sum to |Γ| in each case.
- Multiplicities exact (sympy symbolic), non-negative integers.
- λ₁ values: computed.
- Cross-reference: λ₁(2I) = 168 matches the verified Paper 119 §2.2 value (and Ikeda 1980).
- Rail-placement: OBSERVED.
- Structural interpretation: OBSERVED, pending v1.2 incorporation; not asserted in this artefact.
