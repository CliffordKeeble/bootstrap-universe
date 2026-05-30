# Pre-registration — Sub B (zeta-denominator chance-overlap null)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [../BRIEF.md](../BRIEF.md) (CinC, 30 May 2026, Sub B section)
**Binding**: this document is committed to git **before populations Z and I are finalised** and before any computation runs.

## Population definitions (frozen at pre-registration)

**Population Z (zeta-value denominators, 50 elements)**:
Z_k = denom(ζ(−2k+1)) for k = 1, …, 50, computed via:
  ζ(−n) = −B_{n+1} / (n + 1) for n ≥ 1
  denom(ζ(−n)) = denom(B_{n+1}/(n+1)) — *after* reducing the fraction.

This is the "value denominator" — the denominator of the *value*
ζ(−n) as a reduced fraction, NOT the Bernoulli denominator alone.
The brief is explicit about this convention; Mr A's v0.2 catch on
Paper 203 §7 already flagged "value denominators" vs "Bernoulli
denominators" — we report value-denominators per the brief.

Computation: `sympy.zeta(s).rewrite(...)` or directly via
`sympy.bernoulli(n+1)` then form −B_{n+1}/(n+1) and take .q (denominator).

**Population I (icosahedral / 2I / E₈ meaningful integers, 15 elements)**:
Per brief, exactly:
  I = {1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, 120, 168, 240}
This is the brief's pre-registered list. I freeze it now and do not
amend.

## Match test

The "matches" Paper 203 §7 highlighted:
- denom(ζ(−1)) = 12 (≟ |I|-element 12)
- denom(ζ(−3)) = 120 (≟ |I|-element 120)
- denom(ζ(−7)) = 240 (≟ |I|-element 240)

These are *value* denominators (per the brief). I will recompute them
to verify, then count all matches across Z[1..50] ∩ I.

## Null model

Per brief: under independence of Z and I, the expected number of
matches is approximately Poisson with mean μ where

  μ = |Z| × |I| / M     with     M = max(Z ∪ I).

Then z = (actual − μ) / √μ.

I commit to this null. If Z and I overlap at more than μ + 3√μ, the
matches are STRIKING; in (μ + 2√μ, μ + 3√μ), CONSISTENT WITH FRAMEWORK
ONLY; else CONSISTENT WITH CHANCE.

**Refinement (allowed before run)**: I may sharpen the null by computing
a tighter expected match count. Specifically, given that Z is structured
(value denominators have particular form per von Staudt–Clausen), the
"uniform-density-in-{1..M}" assumption is a *coarse upper bound* on
the chance match probability. A tighter null would compute the actual
density of Z and use that. I commit to computing both:
- **Coarse Poisson null**: μ = |Z| · |I| / M.
- **Tight empirical null**: count how many elements of {1, …, M} actually
  occur in Z (= |Z| up to duplicates), and use that for the conditional
  match expectation. Report both z-scores.

## Anti-circularity check (per brief)

Z is computed from von Staudt–Clausen exactly (Bernoulli numbers via
Sympy). I is the icosahedral / 2I / E₈ meaningful-integer set from
group theory. Neither population is constructed from the other; their
intersection is the empirical match count.

## What I will NOT do

- Add or remove members of I after seeing the match data.
- Switch between value-denominators and Bernoulli-denominators
  conditionally on which gives more matches.
- Adjust the Poisson threshold after seeing the data.

## Files

- `sub_b.py` — implementation of populations Z, I, match count, null.
- `findings_paper_203_sub_b.md` — long-form report.
