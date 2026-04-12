# Golden Zeta Investigation — Findings

Paper 150 §3.1.1 re-examination  
Generated: 2026-04-12  
By: Mr Code  
For: CinC and Cliff  
Programme: Bootstrap Universe / 2I Universe

---

## The question

Does the matching between ζ_φ minima and Riemann zeros degrade
to chance above some t? Cliff remembered: "we did 20 zeros,
came back for more, the zeros weren't right."

## The answer

**No. The matching does not degrade at high t. It is uniformly
modest across the entire range.**

The signal is real but weak. ζ_φ's minima are about 30% closer
to Riemann zeros than uniform random minima would be, and this
holds at t<80, at 80<t<160, and at 160<t<237 with no differential
degradation. The problem is not that the matching breaks down at
high t — the problem is that it was never as strong as the paper's
framing suggested.

## The numbers

### Baseline (N=5000, W=1.0)

| Range | N | Match rate | Mean Δ | Median Δ | Max Δ |
|-------|---|-----------|--------|----------|-------|
| [0, 80) | 21 | 100% | 0.303 | 0.277 | 0.814 |
| [80, 160) | 37 | 97.3% | 0.282 | 0.198 | 1.016 |
| [160, 237) | 42 | 97.6% | 0.326 | 0.172 | 1.259 |

Match rate is 97%+ everywhere. Mean Δ varies from 0.28 to 0.33.
**No collapse at high t.**

### Window sensitivity

As W tightens, match rate drops **uniformly** across all bins:

| W | [0,80) rate | [80,160) rate | [160,237) rate |
|---|-----------|-------------|--------------|
| 1.0 | 100% | 97.3% | 97.6% |
| 0.8 | 95.2% | 94.6% | 90.5% |
| 0.5 | 81.0% | 83.8% | 76.2% |
| 0.3 | 57.1% | 62.2% | 59.5% |
| 0.2 | 38.1% | 51.4% | 57.1% |

No bin collapses preferentially. The [160,237) bin actually
**improves** relative to [0,80) at W=0.2. The degradation Cliff
remembered is not in the window sensitivity — the concern was
about something else.

### Null hypothesis (1000 Monte Carlo trials)

| Range | ζ_φ mean Δ | Null mean Δ | z-score | Verdict |
|-------|-----------|-------------|---------|---------|
| [0, 80) | 0.303 | 0.441 | -1.36 | WEAK |
| [80, 160) | 0.282 | 0.444 | -2.23 | GOOD |
| [160, 237) | 0.326 | 0.448 | -1.68 | WEAK |

ζ_φ beats random at z = -1.4 to -2.2 depending on bin. The best
performance is actually in the **middle bin** [80,160), not in
[0,80). No bin reaches z < -3 (the threshold for "many sigma").

**The matching is better than chance. It is not dramatically
better than chance.**

### Per-zero forensic (delta/null ratio)

| Range | Mean ratio | Median ratio | Zeros closer than null |
|-------|-----------|-------------|----------------------|
| [0, 80) | 0.695 | 0.610 | 15/21 (71%) |
| [80, 160) | 0.636 | 0.412 | 29/37 (78%) |
| [160, 237) | 0.740 | 0.412 | 30/42 (71%) |

Ratio < 1 means ζ_φ is closer than random. About 70-78% of zeros
in every bin are closer than null expectation. Uniform across
range — no high-t collapse.

### Scaling (N = 5000 → 20000)

| N | Mean Δ overall | Mean Δ [160,237) |
|---|---------------|----------------|
| 5000 | 0.305 | 0.326 |
| 10000 | 0.287 | 0.299 |
| 20000 | 0.282 | 0.292 |

Modest improvement. The gap is mostly structural, not a finite-N
artefact. Increasing N by 4x reduces mean Δ by about 10%.

## What this means for Paper 150

### What the paper got right

- ζ_φ's minima do correlate with Riemann zeros. The correlation
  is real and present across the full range [14, 237].
- The golden-angle phase construction produces a function whose
  minima are systematically closer to Riemann zeros than random
  points would be.
- There is no high-t degradation. The concern that triggered this
  investigation is not supported by the data.

### What the paper overstated

- The z-scores are -1.4 to -2.2, not many-sigma. The signal would
  not survive a strict multiple-testing correction across bins.
- The match rate at W=1.0 is high (97%+) but this is partly an
  artefact of the window being comparable to the mean spacing
  between ζ_φ minima (~0.9). A window of W=1.0 when minima occur
  every ~0.9 units of t means most zeros will be within W of
  _some_ minimum by accident.
- The null test reveals this: random minima at the same density
  achieve 89% match rate at W=1.0. ζ_φ achieves 97-100%. The
  improvement over chance is about 8 percentage points in match
  rate, or 30% in mean Δ.

### Recommended correction for §3.1.1

The section should:

1. Add the null comparison. State the z-scores honestly.
2. Frame the result as "ζ_φ minima are systematically closer to
   Riemann zeros than random points at the same density" rather
   than implying a near-perfect correspondence.
3. Remove or qualify any language suggesting the matching is
   dramatically better at t<80 than at t>80 — the data does not
   support differential performance.
4. Note that the signal strength (~2 sigma) is suggestive but
   not conclusive. It warrants further investigation, not a
   strong claim.

## The honest answer to Cliff's question

Cliff asked: does the matching degrade to chance above some t?

**No, it doesn't degrade. But it was never far from chance.**

The matching is uniformly ~30% better than random across the whole
range. That's real. It's interesting. It's worth reporting. But
it's a 1-2 sigma signal, not a 5-sigma discovery. §3.1.1 should
say this clearly.

The first-20 result (mean Δ=0.27) and the extended-to-100 result
(mean Δ=0.30) are essentially the same. The "degradation" Cliff
remembered was probably a misreading of the max-Δ column (which
does grow with range, as expected for any matching process) rather
than the mean-Δ column (which stays flat).

---

*Generated by zeta_phi.py, null_test.py, scaling_test.py, diagnostic.py*  
*Bootstrap Universe Programme*
