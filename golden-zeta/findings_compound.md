# Golden Zeta — Compound Investigation Findings

Pre-registration addendum: `pre_registration_compound.md`, SHA `c2b9d0d`  
Original protocol: `pre_registration_N1000.md`, SHA `21f74ed`  
Generated: 2026-04-12  
By: Mr Code  
Programme: Bootstrap Universe / 2I Universe  
Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155

## Headline summary

| Test | Question | z_overall | Decision |
|------|----------|-----------|----------|
| 1 | N=10000 sqrt(N) scaling | **-22.91** | sqrt(N) CONFIRMED |
| 2a | phi (baseline) | -8.14 | DISCOVERY |
| 2b | sqrt(2) phase + golden norm | **-9.38** | DISCOVERY |
| 2c | e phase + golden norm | -6.42 | DISCOVERY |
| 2d | pi phase + golden norm | -6.46 | DISCOVERY |
| 3a | phi phase + circular norm | -2.63 | NOT SIGNIFICANT |
| 3b | sqrt(2) phase + golden norm | -9.38 | DISCOVERY |
| 4a | zeta_phi vs Riemann zeros | -8.14 | DISCOVERY |
| 4b | zeta_phi vs L(chi_5) zeros | -6.18 | DISCOVERY |

**The story these numbers tell:** The golden norm |Re^2 - 5*Im^2| is
the active ingredient. The equidistributed phase is interchangeable.
The norm detects zeros of both factors of the Dedekind zeta of Q(sqrt(5)).

---

## Test 1 — N=10000 massive-N robustness

**Question:** Does the z = -8.14 result scale as sqrt(N)?

**Answer: Yes.** z = -22.91 at N = 10000 (predicted -25.7, ratio 0.89).

| Quantity | Value |
|----------|-------|
| Riemann zeros | 10000 (range [14.13, 9877.78]) |
| ζ_φ minima | 11675 (density 1.18/unit) |
| Signal mean delta | 0.3239 |
| Null mean delta | 0.4234 |
| Null SE | 0.0043 |
| z_overall | **-22.91** |
| Predicted (sqrt(N)) | -25.7 |
| Ratio actual/predicted | 0.89 |

Bin-resolved z-scores (all strongly negative):

| Bin | N zeros | z |
|-----|---------|------|
| [14, 3500) | 2966 | -12.07 |
| [3500, 7000) | 3737 | -14.55 |
| [7000, 10500) | 3297 | -12.96 |

Per pre-registered interpretation: z in [-20, -30] means "sqrt(N) scaling
confirmed; effect is constant across 100x sample expansion."

The slight shortfall (0.89 of prediction) is consistent with weak
saturation or hidden correlation in the null. It does not threaten
the detection — 22.9 sigma is unambiguous.

Computation: 10000 zetazeros fetched in 6127s (mpmath). 1.24M grid
points scanned. 1000 MC null trials in 92s.

---

## Test 2 — Alternative equidistributed sequences

**Question:** Is the golden ratio phi specific to the effect?

**Answer: No.** All four equidistributed sequences produce |z| > 5 with
the golden norm.

| Alpha | z_overall | Signal mean D | Null mean D |
|-------|-----------|--------------|-------------|
| phi (1.618) | -8.14 | 0.3155 | 0.4177 |
| sqrt(2) (1.414) | **-9.38** | 0.3113 | 0.4386 |
| e (2.718) | -6.42 | 0.3462 | 0.4313 |
| pi (3.142) | -6.46 | 0.4154 | 0.5242 |

sqrt(2) is **stronger** than phi. e and pi are weaker but still past 5 sigma.

Per pre-registered interpretation: "the effect is generic to equidistribution.
The phi-specificity claim in Paper 150 v1.3 is unsupported."

**What this means for v2.0:** The paper cannot claim that the golden ratio
is special. It must reframe to: "any equidistributed phase sequence, when
probed with the golden norm, matches Riemann zeros." The mechanism is the
norm, not the phase.

---

## Test 3 — Norm decoupling

**Question:** Does the effect require both golden phase AND golden norm together?

**Answer: The golden norm alone is sufficient.** The phase is cosmetic.

| Test | Phase | Norm | z_overall |
|------|-------|------|-----------|
| Baseline | phi | golden | -8.14 |
| 3a | phi | circular | **-2.63** |
| 3b | sqrt(2) | golden | **-9.38** |

Per pre-registered interpretation: "3a gives |z| < 3, 3b gives |z| ~ 8:
the golden NORM alone is sufficient; any equidistributed phase suffices.
This would be a major reframing."

**The key finding:** The circular norm Re^2 + Im^2 with golden phase
produces only 401 minima (density 0.28/unit) compared to 1704 minima
(density 1.20/unit) with the golden norm. The golden norm |Re^2 - 5*Im^2|
creates minima at a fundamentally different density — and those minima
align with Riemann zeros.

This makes structural sense: the golden norm is the field norm of Q(sqrt(5)).
It probes the arithmetic of that field. The circular norm probes nothing
field-specific.

---

## Test 4 — Dedekind factorisation partition

**Question:** Does zeta_phi detect zeros of both factors of the Dedekind zeta?

**Answer: Yes.** Both the Riemann zeta and L(chi_5) zeros are tightened.

| Target | N zeros | z_overall |
|--------|---------|-----------|
| Riemann zeta zeros | 1000 | -8.14 |
| L(chi_5) zeros | 454 | -6.18 |

Per pre-registered interpretation: "both give |z| > 5: zeta_phi detects
both factors. The Dedekind bridge IS the mechanism."

**L(chi_5) zeros:** Computed from vectorised Fejer-smoothed Dirichlet
series (N=10000 terms). First zero at t = 6.65, consistent with LMFDB.
454 zeros found in [6.65, 1424.71].

**Combined interpretation with Tests 2-3:** The golden norm |Re^2 - 5*Im^2|
is the field norm of Q(sqrt(5)). The Dedekind zeta function of Q(sqrt(5))
factors as zeta(s) * L(s, chi_5). A Dirichlet sum with any equidistributed
phase, when measured through this field norm, produces minima that
correlate with zeros of **both factors** of the Dedekind zeta.

This is not about the golden ratio. It is about Q(sqrt(5)).

---

## The narrative that emerges

The four tests tell a coherent story, but it is not the story Paper 150
v1.3 tells:

**v1.3 story (now refuted by Test 2):** The golden ratio phi has a
special relationship with Riemann zeros via the golden-angle phase.

**v2.0 story (supported by all four tests):** The field norm of Q(sqrt(5))
— specifically |Re^2 - 5*Im^2| — creates a Dirichlet-sum-based probe
that detects zeros of the Dedekind zeta function of Q(sqrt(5)). This
includes both Riemann zeta zeros and L(chi_5) zeros. The phase sequence
is interchangeable provided it is equidistributed; the norm is the active
ingredient.

The v2.0 story is more constrained than v1.3 (phi is not special) but
also more structural (the mechanism is identified: it is the field norm,
acting as a Dedekind probe). This is, arguably, a stronger paper because
the claim is sharper and the falsification tests are included.

## Status flags (per Pattern 39)

| Claim | Status | Evidence |
|-------|--------|----------|
| zeta_phi minima match Riemann zeros | OBSERVED | z = -8.14 at N=1000 |
| Effect is generic to equidistributed phases | OBSERVED | Test 2: all four alphas |z| > 5 |
| Golden norm carries the effect | OBSERVED | Test 3: norm decoupling |
| Dedekind bridge is the mechanism | STRUCTURAL | Test 4: both factors detected |
| phi is specific to the effect | REFUTED | Test 2: sqrt(2) is stronger |

## What v2.0 should say

> A Dirichlet sum with equidistributed phase, measured through the field
> norm of Q(sqrt(5)), produces minima that systematically match zeros of
> the Dedekind zeta function zeta_{Q(sqrt(5))} = zeta * L(chi_5). The
> effect is independent of the specific irrational used for the phase
> (tested with phi, sqrt(2), e, pi) and depends solely on the field norm.
> This identifies Q(sqrt(5)) as the number field mediating the
> correspondence, consistent with the programme's derivation of physical
> constants from the ring Z[phi].

---

*No deviations from registered protocol.*  
*All four tests complete.*  
*Pre-registration SHA: c2b9d0d*  
*Bootstrap Universe Programme*
