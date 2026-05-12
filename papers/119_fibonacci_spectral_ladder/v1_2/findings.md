# Paper 119 v1.2 — findings (Tasks A and B)

Two pre-registered computational tasks (see `task_brief.md`) executed against
the criteria recorded before any computation was run. Decisions reported
against the pre-registered branches; no retro-fitting.

- Pre-registration commit: `72c0775`
- Audit gap: ≥ 30 min between pre-registration and results commit (computations
  actually executed in between)

---

## Task A — Null test for rail coverage

### Result

Detail in `task_a_results.md` and `task_a_results.json`. Headline numbers:

| Sequence | Cov [9,13] | Precision [9,13] |
|---|---|---|
| **Rails** | 1/1 = 1.00 | 1/2 = 0.50 |
| Squares | 0/1 = 0.00 | 0/1 = 0.00 |
| Triangulars | 0/1 = 0.00 | 0/1 = 0.00 |
| **Primes** | 1/1 = 1.00 | 1/2 = 0.50 |
| Combined | 1/1 = 1.00 | 1/4 = 0.25 |

Random null over 10⁵ trials of size-10 subsets of [1, 13]:

| | Coverage [1,13] | Coverage [9,13] |
|---|---|---|
| Rails | 1.0000 | 1.0000 |
| Random mean | 0.7690 | 0.7680 |
| P(random ≥ rails) | **0.0138** | **0.7680** |

### Decision against pre-registered criterion

The brief pre-registers three branches:

1. **Strong claim survives** if rails are *strictly* better than all comparison
   sequences in [9, 13] **and** in the top 5% of the random null.
2. **Strong claim fails** if rails are *matched or beaten* by any comparison
   sequence in [9, 13].
3. **Equivocal** if rails beat individual sequences but fall in the bulk of the null.

**Outcome: branch 2 — the strong claim fails by the pre-registered criterion.**

Rails are **tied** with primes on both coverage and precision in [9, 13]:
both predict `{11, 13}` against the realised hit `{13}`, both score 1/2.
Primes are not a constructed rail and were chosen pre-registration as a
control sequence. The criterion explicitly counts a tie as a failure ("matched
or beaten").

This is exactly the concern Mr A raised in §9.2: in [1, 8] the rails cover
trivially because *any* dense small-integer set covers [1, 8]; in [9, 13]
the rails make two predictions (`L_5 = 11`, `F_7 = 13`), only one of which
(`F_7 = 13`) is realised. The same call — `{11, 13}` predicted, `{13}` hit
— is made by the primes, with no need for the rail superstructure.

### Auxiliary observation (full range)

In [1, 13], rails do beat the null statistically: P(random ≥ rails) = 0.0138,
i.e. rails are in the top 1.38% of the null distribution. Rails also strictly
beat all individual comparison sequences in [1, 13] (rails 9/9, combined
8/9, primes 5/9, triangulars 3/9, squares 2/9). But the full range is the
range Mr A identified as *not* the discriminating test — it's the range where
rails win trivially. The discriminating range is the one where the criterion
was set, and the criterion fails there.

### Status flags

- Coverage rates and precisions: OBSERVED.
- Random-null distribution: OBSERVED (seed 119, 10⁵ trials).
- Decision (strong-claim fails): OBSERVED against pre-registered criterion;
  no methodology was altered between pre-registration and run.

---

## Task B — Spectral gap of S³/2T and S³/2O

### Result

Detail in `task_b_results.md` and `task_b_results.json`. Headline:

| Group | Order | l (smallest l > 0 with trivial sub) | λ₁ = l(l+2) | Rail placement |
|---|---|---|---|---|
| **2T** | 24 | **6** | **48** | **Lucas rail** (L₄² − 1 = 7² − 1 = 48) |
| **2O** | 48 | **8** | **80** | **OFF ALL RAILS** |
| 2I (verification) | 120 | 12 | 168 | F₇² − 1 = 13² − 1 = 168 ✓ |

The 2I verification matches Paper 119 §2.2 (and Ikeda 1980), confirming the
Frobenius-reciprocity pipeline is sound. Class-size totals (24, 48, 120)
check, and multiplicities are exact non-negative integers throughout l = 0…20
(symbolic computation in sympy).

### Decision against pre-registered criterion

The brief pre-registers three branches:

1. **Generalisation** if both λ₁(2T) and λ₁(2O) sit on at least one rail.
2. **Icosahedral specificity** if at least one of λ₁(2T), λ₁(2O) lands off all rails.
3. **Partial generalisation** if one of 2T, 2O lands on a rail and the other
   doesn't — "reported as such, with interpretation deferred to v1.2 discussion."

**Outcome: branch 3 — partial generalisation.**

2T lands cleanly on the **Lucas rail**: λ₁(2T) = 48 = 7² − 1 = L₄² − 1.
This is a non-trivial finding — the Lucas rail picks up the binary tetrahedral
spectral gap with the same `L_n² − 1` arithmetic that produced `L_5 = 11`
in the 2I context. (Note: the brief's pre-stated 2I rail-placement uses the
Fibonacci rail at λ₁ = 168 = F₇² − 1. So 2T → Lucas, 2I → Fibonacci.)

2O is genuinely off-rail: λ₁(2O) = 80, which is not of the form `F_n² − 1`,
not `L_n² − 1`, not 35.

### Interpretation note (status-flagged)

Per brief, structural interpretation is deferred to v1.2 paper discussion.
Two observations are worth flagging without asserting interpretation here:

- The rail family `{F_n² − 1, L_n² − 1, 35}` does not exhaust 2O. Either the
  rail family is icosahedral-plus-tetrahedral but not octahedral, or it is
  missing a generator that 2O would lie on (an "octahedral closure" analogous
  to face-closure {6}).
- The fact that 2T and 2I land on *different* rails (Lucas and Fibonacci
  respectively), but both Fibonacci-family in the broad sense, is itself a
  pattern that may matter for v1.2 framing.

These are observations, not claims.

### Status flags

- Frobenius reciprocity: DERIVED (standard rep theory).
- Conjugacy class data: verified-against-reference (Coxeter, Cisneros-Molina 2000).
  2T: 7 classes summing to 24. 2O: 8 classes summing to 48. 2I: 9 classes summing to 120.
- Multiplicities: exact (sympy symbolic), confirmed non-negative integers.
- 2I cross-reference (λ₁ = 168, F₇² − 1): matches Paper 119 §2.2.
- λ₁(2T) = 48, λ₁(2O) = 80: OBSERVED-and-cross-referenced (matches published
  values for binary-polyhedral spherical space form spectra, e.g. Ikeda 1980 /
  Bär).
- Rail placement decisions: OBSERVED against pre-registered criterion.
- Structural interpretation: deferred to v1.2 paper, not asserted here.

---

## Summary for v1.2 incorporation

- **Task A**: null test **fails the strong claim** (rails matched by primes
  in [9, 13]). Suggests v1.2 §9.2 must reframe the rail claim: rails carry
  predictive content over the full range [1, 13] (top 1.38% of null), but
  in the discriminating sub-range [9, 13] they are not distinguishable from
  the prime sequence. Mr A's concern is sustained by the computation.
- **Task B**: rail structure **partially generalises**. 2T sits on the Lucas
  rail (λ₁ = 48 = L₄² − 1); 2O is off all rails (λ₁ = 80). v1.2 §9.6 should
  promote 2T/2O from open question to derived comparison, and the framing
  should acknowledge that the rail family captures 2T and 2I but not 2O —
  a finer phenomenon than uniform icosahedral specificity, and finer than
  uniform binary-polyhedral generalisation.

Both results are reported as they came out of the pre-registered scripts.
No reframing has been applied between pre-registration and these findings.

🐕☕⬡
