# Mr Code brief — icosahedral lag simulation v0.1 (space from time)

**Date**: 28 May 2026
**Brief author**: CinC
**Relationship to prior work**: This is a NEW investigation track, not a continuation of the 2I-network static line (v0.1 commit 06cbd9b, v0.2 commit 32dc678). Those tested whether a *static* graph carries S³/2I spectral structure; both returned clean FAILs (trees, then diverging refinement, then null-indistinguishability). The diagnosis was that they bolted space on by hand. This track tests the opposite: whether spatial geometry *emerges* from temporal coupling of icosahedral units.
**Pre-registration document**: commit to git BEFORE any computation.
**Estimated runtime**: exploratory; aim for ≤30 min compute. Likely well under.

---

## The conjecture

**S³ geometry emerges from the temporal coupling (lag) of parallel icosahedral units. Distance is encoded as delay; the metric is a consequence of how icosahedral dynamics synchronise across lag, not a pre-imposed structure.**

This is the "closure in time" insight made computational. The static line asked "does this fixed spatial graph have S³ structure?" This asks "does space emerge from time?" The unit of dynamics carries icosahedral (A₅ / 2I) symmetry; the coupling carries delay; the claim is that the emergent correlation geometry is S³-like, and that the icosahedral symmetry is load-bearing (non-icosahedral units won't do it).

## The critical anti-circularity requirement

**The coupling topology must NOT pre-impose the geometry we are trying to find.** This is the single most important constraint in the brief. If lag τ_ij is set proportional to a distance in some pre-chosen space, the simulation is circular — it will "find" the geometry it was handed. To avoid this:

- Coupling is **all-to-all** (complete graph) or **random** (Erdős–Rényi), NOT a spatial lattice.
- Lag is a **single global value τ**, or randomly assigned per edge, NOT distance-encoding.
- Geometry must be **reconstructed from the steady-state dynamics** (correlations), never read off the coupling structure.

If Mr Code cannot see how to satisfy this for a given construction choice, STOP and flag — a circular construction is worse than no result.

## Construction

**Dynamical unit (primary):** each of N units holds a state z_i on the Riemann sphere (ℂ ∪ {∞}). Each unit iterates an **A₅-equivariant rational map** R(z) — Klein's icosahedral dynamics. Use a documented icosahedrally-equivariant map (e.g. the Doyle–McMullen iteration map, or one built from the icosahedral invariants f₁₂, H₂₀, T₃₀); the exact map matters less than that it commutes with the A₅ action. Document the choice in the pre-registration.

**Dynamical unit (alternative, if the Riemann-sphere map proves awkward):** each unit holds a discrete 2I element g_i ∈ 2I (120 states). Dynamics updates g_i by left-composition with a 2I element selected from lagged neighbour states. Cleaner state space, discrete. Mr Code may use this if the rational-map version is numerically fragile; document which was used and why.

**Coupling with lag (discrete time):** at each integer step t, each unit updates from its own current state and the **lagged** states of its coupling partners:

```
z_i(t+1) = R(z_i(t)) ⊕ ε · Σ_j K_ij · coupling(z_j(t − τ), z_i(t))
```

where ⊕ is an appropriate combination on the Riemann sphere (e.g. Möbius nudge toward the coupled value, kept on-sphere), ε is coupling strength, K_ij is the (all-to-all or random) coupling matrix, and τ is the global lag. The coupling function pulls a unit toward (or repels from) its partners' past states. Mr Code chooses a sensible on-sphere coupling and documents it.

**Run to steady state / attractor.** Discard a transient, then record states over a window.

**Reconstruct effective geometry:**
1. Compute pairwise correlation C_ij between units over the recording window (e.g. time-averaged proximity on the sphere, or phase-locking value).
2. Convert to effective distance: d_ij = arccos(C_ij) or d_ij = −log(C_ij) (document choice).
3. From {d_ij}, compute the emergent-geometry observables below. (Embedding/MDS optional for visualisation, not required for the verdict.)

## Pre-registered observables and thresholds

**Observable 1 — spectral dimension d_s.** Build a weighted graph from the effective distances (or a proximity graph), run a diffusion / random walk, and extract the spectral dimension from the return-probability scaling p(t) ~ t^(−d_s/2).
- Target: **d_s ∈ [2.5, 3.5]** (S³ is 3-dimensional).
- A tree gives d_s near 1–2 and unbounded; flat ℝ³ gives 3 but non-compact; we want 3 AND compact (Observable 2).

**Observable 2 — compactness / volume saturation.** Count units within effective radius r of a typical unit, as a function of r.
- Target: growth ~ r³ at small r, then **saturation** (the count plateaus at N because the space closes). S³ is compact; flat space would not saturate.
- Pass: clear r³ regime followed by saturation. Fail: power-law growth with no saturation (non-compact), or no clean r³ regime.

**Observable 3 — curvature sign.** Use a coarse Ollivier–Ricci or triangle-comparison estimate on the effective-distance graph.
- Target: **positive** mean curvature (S³ is positively curved).
- Informational-to-pass: positive supports S³; flat/negative counts against.

**Observable 4 (secondary) — spectral gap.** If and only if Observables 1–3 pass, compute the Laplacian spectral gap of the emergent structure and compare (volume-normalised) to the S³/2I target ~50.4 (= 168·(2π²/120)^(2/3), the v0.2 corrected target).
- This is a SECONDARY check, deliberately not the headline — it's exactly what the static line couldn't reach, so it must not gate the verdict. Report it; don't let it drive PASS/FAIL.

## Pre-registered null hypotheses

Run the identical pipeline for each null and compare:

1. **Non-icosahedral units.** Replace the A₅-equivariant map with a generic (non-symmetric) rational map of comparable degree, or generic oscillators. If the nulls also produce d_s ≈ 3 + saturation, the icosahedral structure is NOT load-bearing and the result is generic coupled-oscillator behaviour.
2. **No delay (τ = 0).** Same icosahedral units, instantaneous coupling. If S³ structure appears without lag, then "space from time" is not what's producing it.
3. **Shuffled coupling.** Randomise K_ij while preserving its statistics. Tests whether any specific coupling accident is responsible.

**Distinguishability requirement:** the icosahedral-with-lag run must beat ALL THREE nulls on Observables 1 and 2 by a clear margin (pre-register: d_s closer to 3, and saturation present where nulls lack it). If the nulls match it, the conjecture is not supported.

## Stop-on-fail protocol

**Stop and report failure if:**
- No steady state / attractor is reached (the dynamics blow up or wander without settling) after reasonable transient — report as FAIL–no-attractor.
- Observable 1 (spectral dimension) falls outside [2.5, 3.5] at the natural N — report as FAIL–wrong-dimension.
- Observable 2 shows no saturation (non-compact) — report as FAIL–non-compact.
- Any null matches the icosahedral run on Observables 1–2 — report as FAIL–null-indistinguishable.

**Continue / PASS if:** Observables 1–3 pass AND all three nulls are beaten clearly. Then compute Observable 4 and report it as supplementary.

Do not tune ε, τ, or N after seeing the verdict-relevant data to rescue a fail. You may do a SMALL pre-data sweep of ε and τ to find a regime where the dynamics reaches a steady state at all (this is finding the operating point, not fitting the result) — document the sweep, fix the operating point, THEN run the verdict pipeline and nulls at the fixed point.

## Implementation notes

- Python; numpy/scipy. Riemann-sphere arithmetic via stereographic coordinates or unit-3-vector representation (avoid ∞ blow-ups — the 3-vector form on S² is numerically safest).
- N: start small (N=60 or 120, echoing |A₅|, |2I|) for the operating-point sweep; run the verdict at the largest N that's comfortable (aim N≈500–1000 if compute allows).
- Delay buffer: store the last τ steps of all states (a ring buffer); straightforward.
- Spectral dimension: standard return-probability method on the weighted graph; document the walk normalisation.
- Keep everything reproducible: fixed seed, log all parameters.
- This is exploratory — if something surprising but coherent emerges that isn't S³, report it; an unexpected coherent structure is itself a finding (Pattern 19).

## Deliverable

A markdown report (`findings_icosahedral_lag_v0_1.md`) committed alongside the pre-registration, structured as:

```
## Results

### Construction choices
- Which dynamical unit (Riemann-sphere map / discrete 2I) and why
- The A₅-equivariant map used (cite source or give the invariants)
- Coupling function, ε, τ, K_ij type, N — the fixed operating point
- The operating-point sweep that found it

### Observables table (icosahedral-with-lag vs the three nulls)
| Run | d_s | r³ regime? | saturation? | curvature sign | λ₁ (norm) |
|---|---|---|---|---|---|
| icosahedral + lag | ... | ... | ... | ... | ... |
| null: non-icosahedral | ... | ... | ... | ... | — |
| null: no delay | ... | ... | ... | ... | — |
| null: shuffled coupling | ... | ... | ... | ... | — |

### Verdict
One of:
- PASS — S³-like structure emerges, icosahedral symmetry load-bearing, beats all nulls
- PARTIAL — some observables pass, specifics reported
- FAIL–no-attractor / FAIL–wrong-dimension / FAIL–non-compact / FAIL–null-indistinguishable

### Honest limitations
- Circularity self-audit: confirm the coupling did NOT pre-impose geometry, and explain how you know
- What's not tested; where the construction might be doing the work instead of the physics

### Pattern flags
- Pattern 75 (null): how each null was satisfied
- Pattern 39 (DERIVED vs OBSERVED): the A₅ map is DERIVED; the emergent geometry is OBSERVED
- Pattern 19 (adversary): first attack a hostile reviewer mounts — almost certainly "your coupling secretly encodes the geometry." The anti-circularity section is the pre-emptive answer; make sure the report earns it.
```

## Hygiene

- The anti-circularity requirement is the binding constraint. A PASS from a circular construction is a FALSE pass and worse than a FAIL. Self-audit explicitly.
- Pre-registration thresholds (d_s ∈ [2.5, 3.5], saturation required, beat all three nulls) are fixed once committed.
- Do not read Paper 203, the conversation that generated this brief, or the 2I-network findings. The construction here is independent.
- Report negative results with the same emphasis as positive. Two clean FAILs precede this; a third honest FAIL is just as valuable as a PASS.
- If the dynamics simply won't settle, that itself may be programme-consistent ("S³/2I not stable for small networks") — report it as FAIL–no-attractor with that note, don't force it.

⌨️ over to you.

---

## Mr Code's report

*[To be filled in by Mr Code.]*
