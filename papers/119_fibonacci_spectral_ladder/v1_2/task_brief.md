# Mr Code Brief — Paper 119 v1.2 Computations

**Branch:** `paper-119-v1-2-computations`
**Folder:** `papers/119_fibonacci_spectral_ladder/v1_2/`
**Workflow:** pre-registration commit → run → results commit → PR → Cliff review → merge → CinC incorporates into v1.2
**Standing trust:** Green for both tasks (small computations, no external API calls, no side effects beyond local files)

---

## Context

Paper 119 (*The Fibonacci Spectral Ladder*) has been through two Mr Adversary v2.3 passes. The round-2 review awards 3 stars provisionally and explicitly asks for two computational checks before signing off on a 4-star paper:

1. **Mr A §9.2 — null test for rail coverage.** The combined rails Fibonacci ∪ Lucas ∪ {6} cover 10 of 13 integers in [1, 13]; the realised thresholds occupy 9 of those 10. Mr A's principal finding is that 8 of the 9 realised thresholds sit in [1, 8] where the rails trivially cover *any* dense small-integer sequence — so the framework's actual discriminating range is [9, 13], where it makes two non-trivial predictions (L₅ = 11, F₇ = 13) of which only one (F₇ = 13) is realised. Mr A: *"the framework's claim to discriminative content cannot be assessed without [the null test]. I would not let the paper leave the desk without it."*

2. **Mr A §9.6 — spectral gap of S³/2O and S³/2T.** Tests whether the Fibonacci–Lucas–face-closure rail structure is icosahedrally specific or generalises to other binary polyhedral space forms. Mr A: *"essential rather than open... cannot in good conscience let you postpone."*

Both tasks are small. Both feed the same v1.2 publication.

---

## Pre-registration discipline

**Pre-registration commit (before running either task):** this brief committed to the branch. The pre-registration captures: the methodology, the comparison sequences, the output format, and the criteria for "rail-coverage beats null" / "spectral gap lands on a rail."

After the pre-registration commit, no edits to this brief until results are in. Results commit follows in a separate commit. Any methodology change discovered necessary during execution gets a separate amendment commit *before* the results commit, with rationale.

**Audit gap:** the time between the pre-registration commit and the results commit is the audit moment. It must be long enough to demonstrate non-retrofitting — at least 30 minutes between commits, computations actually run in between.

---

## Task A — Null test for rail coverage

### Question

Does the rail structure (Fibonacci ∪ Lucas ∪ {6}) have non-trivial predictive power against comparable-density baseline sequences, particularly in the range [9, 13] where rails make discriminating predictions?

### Inputs (pre-registered, fixed)

**Realised data** (the nine first-appearance thresholds at which irreducible 2I representations first appear in the SU(2) → 2I branching, with l+1 reported):

```python
DATA = {1, 2, 3, 4, 5, 6, 7, 8, 13}
```

**Rail structure** (in [1, 13]):

```python
FIBONACCI = {1, 2, 3, 5, 8, 13}      # F_n for n ≥ 1, intersected with [1, 13]
LUCAS = {1, 2, 3, 4, 7, 11}           # L_n for n ≥ 0, intersected with [1, 13]
FACE = {6}                            # D! = |S_3| = 6
RAILS = FIBONACCI | LUCAS | FACE      # = {1, 2, 3, 4, 5, 6, 7, 8, 11, 13} = 10 values
```

**Range partition:**

```python
RANGE_FULL = set(range(1, 14))        # [1, 13], 13 values
RANGE_TRIVIAL = set(range(1, 9))      # [1, 8], 8 values — rails cover this trivially
RANGE_DISCRIMINATING = set(range(9, 14))  # [9, 13], 5 values — non-trivial test range
```

### Comparison sequences (pre-registered, Mr A's specific suggestions)

Three structured small-integer sequences, individually and combined:

```python
SQUARES = {n*n for n in range(1, 4)}                      # {1, 4, 9} — 3 values
TRIANGULARS = {n*(n+1)//2 for n in range(1, 5)}           # {1, 3, 6, 10} — 4 values
PRIMES = {2, 3, 5, 7, 11, 13}                             # 6 values
COMBINED_STRUCTURED = SQUARES | TRIANGULARS | PRIMES      # — should be 11 values
```

Plus a randomisation-based null:

```python
RANDOM_NULL: for N_TRIALS = 100_000, sample subsets of [1, 13]
             of size |RAILS| = 10, uniformly without replacement
```

### Methodology

For each candidate sequence S (rails, squares, triangulars, primes, combined, and each random sample):

1. Compute **coverage of DATA by S**:
   - In RANGE_FULL: |DATA ∩ S| / |DATA|
   - In RANGE_TRIVIAL: |DATA ∩ S ∩ RANGE_TRIVIAL| / |DATA ∩ RANGE_TRIVIAL|
   - In RANGE_DISCRIMINATING: |DATA ∩ S ∩ RANGE_DISCRIMINATING| / |DATA ∩ RANGE_DISCRIMINATING|

2. Compute **predictive precision** (S as a predictor of DATA):
   - In RANGE_FULL: |DATA ∩ S| / |S|
   - In RANGE_DISCRIMINATING: |DATA ∩ S ∩ RANGE_DISCRIMINATING| / |S ∩ RANGE_DISCRIMINATING| (if denominator > 0)

3. For the random null: compute distribution of full-range coverage rate and discriminating-range coverage rate. Report mean, std dev, and the empirical p-value for the rail-structure result against the random distribution.

### Pre-registered output format

A single results table `task_a_results.md` with three sub-tables:

**Sub-table 1: Coverage by sequence (recall: how much of DATA does each sequence cover?)**

| Sequence | Size | Cov [1,13] | Cov [1,8] | Cov [9,13] |
|---|---|---|---|---|
| Rails | 10 | 9/9 | 8/8 | 1/1 |
| Squares | 3 | ?/9 | ?/8 | ?/1 |
| Triangulars | 4 | ?/9 | ?/8 | ?/1 |
| Primes | 6 | ?/9 | ?/8 | ?/1 |
| Combined structured | 11 | ?/9 | ?/8 | ?/1 |

**Sub-table 2: Precision by sequence (recall: how predictive is each sequence?)**

| Sequence | Predictions in [9,13] | Hits in [9,13] | Precision in [9,13] |
|---|---|---|---|
| Rails | {11, 13} | {13} | 1/2 = 0.50 |
| Squares | {9} | {} | 0/1 = 0.00 |
| Triangulars | {10} | {} | 0/1 = 0.00 |
| Primes | {11, 13} | {13} | 1/2 = 0.50 |
| Combined | {9, 10, 11, 13} | {13} | 1/4 = 0.25 |

**Sub-table 3: Random null distribution (10⁵ trials of random size-10 subsets of [1, 13])**

| Statistic | Coverage [1,13] | Coverage [9,13] |
|---|---|---|
| Mean | ? | ? |
| Std dev | ? | ? |
| Rail percentile | ? | ? |
| Empirical p-value | ? | ? |

Plus a `task_a_results.json` with the raw distribution data for reproducibility.

### Status flags

- The realised data and rail structure are DERIVED (from the branching computation in Paper 119 §5).
- The choice of comparison sequences (squares, triangulars, primes, random size-10) is pre-registered before execution.
- The resulting coverage rates and p-values are OBSERVED.
- Any structural interpretation ("rails beat null", "rails fail null") is OBSERVED-with-empirical-support.

### Decision criterion (pre-registered)

**Strong claim survives** if rails achieve coverage rates and precision in RANGE_DISCRIMINATING strictly better than all comparison sequences AND in the top 5% of the random null distribution.

**Strong claim fails** if rails are matched or beaten by any comparison sequence in RANGE_DISCRIMINATING (the rails contribute no additional predictive content beyond what comparable-density sequences provide).

**Equivocal** if rails beat individual sequences but fall within the bulk of the random null distribution.

The decision is reported with the results, not retro-fitted to support the framing.

---

## Task B — Spectral gap of S³/2T and S³/2O

### Question

Where do the spectral gaps of the other two binary polyhedral space forms land? Specifically: does λ₁(S³/2T) and λ₁(S³/2O) sit on a Fibonacci rung (Fₙ² − 1), a Lucas rung (Lₙ² − 1), or the face-closure rung (D!² − 1 = 35)? Or neither?

### Background (DERIVED from standard rep theory)

For a binary polyhedral group Γ, the spectral gap of S³/Γ is l(l+2) at the smallest l > 0 such that the trivial Γ-representation appears in the SU(2) restriction D_{l/2}|_Γ. By Frobenius reciprocity:

```
mult(triv_Γ, D_{l/2}|_Γ) = (1/|Γ|) Σ_C |C| · χ_{l/2}(α_C)
```

where C runs over conjugacy classes of Γ, |C| is the class size, α_C is the half-angle of the SU(2) element in C, and χ_{l/2}(α) = sin((l+1)α) / sin(α) is the SU(2) character (with χ_{l/2}(0) = l+1 by limit).

### Group data (pre-registered, from standard references)

**Binary tetrahedral 2T, order 24, 7 conjugacy classes:**

| Class | Size | Half-angle α | Order |
|---|---|---|---|
| {1} | 1 | 0 | 1 |
| {−1} | 1 | π | 2 |
| C3 | 4 | 2π/3 | 3 |
| C3' | 4 | 4π/3 | 3 |
| C6 | 4 | π/3 | 6 |
| C6' | 4 | 5π/3 | 6 |
| C4 | 6 | π/2 | 4 |

(Total: 1+1+4+4+4+4+6 = 24 ✓)

**Binary octahedral 2O, order 48, 8 conjugacy classes:**

| Class | Size | Half-angle α | Order |
|---|---|---|---|
| {1} | 1 | 0 | 1 |
| {−1} | 1 | π | 2 |
| C4 (rot, full) | 6 | π/2 | 4 |
| C3 | 8 | 2π/3 | 3 |
| C6 | 8 | π/3 | 6 |
| C8 | 6 | π/4 | 8 |
| C8' | 6 | 3π/4 | 8 |
| C2 (refl-like) | 12 | π/2 | 4 |

(Total: 1+1+6+8+8+6+6+12 = 48 ✓ — verify class sizes against a reference; this is from memory and one of the order-4 classes may have a different half-angle)

**Implementation note:** Mr Code, please **verify the conjugacy class data against an authoritative source** (Coxeter or any standard binary polyhedral group reference) before computing. The 2I data in the published paper §2.2 is verified; the 2T/2O tables above are reconstructed and may have small errors — if anything looks off, fix from source and note the discrepancy.

### Methodology

For each Γ ∈ {2T, 2O}:

1. Verify the conjugacy class data against reference (Coxeter; Cisneros-Molina 2000; or equivalent).
2. For l = 1, 2, 3, ..., 20: compute

```python
def chi(l, alpha):
    if alpha == 0:
        return l + 1
    return math.sin((l+1) * alpha) / math.sin(alpha)

def mult_trivial(l, classes):
    return sum(size * chi(l, alpha) for size, alpha in classes) / sum(size for size, alpha in classes)
```

3. Find the smallest l > 0 such that `mult_trivial(l, Γ) ≥ 1` (numerically: ≥ 1 − 10⁻⁶).
4. Report λ₁ = l(l+2).
5. Cross-reference with Ikeda 1980 [ref 1] or Cisneros-Molina 2000 published values if available — note any discrepancy.
6. Check whether λ₁ ∈ {Fibonacci² − 1} ∪ {Lucas² − 1} ∪ {35} (face-closure):
   - Fibonacci² − 1 = {0, 0, 3, 8, 24, 63, 168, 440, 1155, ...} for F_n with n = 1, 2, 3, ...
   - Lucas² − 1 = {3, 0, 8, 15, 48, 120, 323, 840, ...} for L_n with n = 0, 1, 2, ...
   - {35}

### Pre-registered output format

A results table `task_b_results.md`:

| Group | Order | l (smallest l > 0 with trivial-rep) | λ₁ = l(l+2) | On Fibonacci rail? | On Lucas rail? | On face-closure rail? | Off all rails? |
|---|---|---|---|---|---|---|---|
| 2T | 24 | ? | ? | ? | ? | ? | ? |
| 2O | 48 | ? | ? | ? | ? | ? | ? |
| 2I (verification) | 120 | 12 | 168 | F₇² − 1 ✓ | — | — | — |

Plus a `task_b_results.json` with the full multiplicity computation `mult_trivial(l, Γ)` for l = 0, ..., 20 for each group, so the gap-finding step is auditable.

### Status flags

- The Frobenius formula computation is DERIVED.
- The conjugacy class data is verified-against-reference.
- The values λ₁ for 2T and 2O are computed-and-cross-referenced.
- Rail-placement (whether λ₁ lands on Fibonacci/Lucas/face) is OBSERVED.
- Any structural interpretation ("rail structure generalises" / "rail structure is icosahedrally specific") is OBSERVED-and-pending-Paper-119-v1.2-incorporation, not asserted at the brief level.

### Decision criterion (pre-registered)

**Generalisation** if both λ₁(2T) and λ₁(2O) sit on at least one of the three rails. The rail structure then applies to all three binary polyhedral space forms, weakening icosahedral specificity but revealing a wider phenomenon. Paper 119 v1.2 would reframe the rail claim accordingly and a follow-up paper on the binary polyhedral rail structure becomes a natural candidate.

**Icosahedral specificity** if at least one of λ₁(2T), λ₁(2O) lands off all rails. The rail structure is then privileged to 2I, strengthening the icosahedral case structurally. Paper 119 v1.2 retains the rail claim but emphasises icosahedral exclusivity.

**Partial generalisation** if one of 2T, 2O lands on a rail and the other doesn't. Reported as such, with interpretation deferred to Paper 119 v1.2 discussion.

The decision is reported with the results, not retro-fitted.

---

## Workflow

1. **Pre-registration commit**: Mr Code commits this brief to `paper-119-v1-2-computations` branch under `papers/119_fibonacci_spectral_ladder/v1_2/task_brief.md`. Commit message: `pre-reg: Paper 119 v1.2 computational tasks A (null test) and B (2T/2O spectral gap)`.

2. **Implementation**: Mr Code writes the two Python scripts under the same folder:
   - `task_a_null_test.py`
   - `task_b_2T_2O_gap.py`

3. **Execution**: Mr Code runs both, produces `task_a_results.md`, `task_a_results.json`, `task_b_results.md`, `task_b_results.json`.

4. **Results commit**: separate commit. Commit message: `results: Paper 119 v1.2 tasks A and B — [one-line summary of each outcome]`.

5. **PR**: Mr Code opens PR against main with summary of findings and any methodology amendments.

6. **Review**: Cliff reviews (or CinC if Cliff delegates).

7. **Merge**: into main once reviewed.

8. **Downstream**: CinC incorporates results into Paper 119 v1.2 §9.2 (null test in body) and §9.6 (2O/2T comparison, promoted from open question to derived comparison). Mr Library publishes v1.2.

---

## Notes for Mr Code

- **Both tasks are small.** Task A is ~50 lines of Python. Task B is ~80 lines. Both should run in seconds.
- **No external API calls required.** Mr A's reference to Ikeda 1980 [ref 1] is for cross-reference; if Ikeda's published values for 2T/2O gaps aren't immediately available, the Frobenius computation IS the verification.
- **If you find a structural issue** with the brief during execution (e.g., the 2T/2O conjugacy class data has an error and the published table doesn't match Coxeter), commit an amendment to the brief explaining the issue *before* committing results. This preserves the pre-registration audit.
- **Output structure matters more than fast iteration.** Mr A will read the v1.2 paper that cites these results; please produce clean tables and reproducible JSON.
- **Status flags should be visible in the results documents**, not just in this brief.

---

🐕☕⬡
