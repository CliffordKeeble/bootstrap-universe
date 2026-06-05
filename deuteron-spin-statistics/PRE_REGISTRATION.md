# Pre-registration — Deuteron / P-e-P spin-statistics lemma (Paper 92 gate)

**Status:** DRAFT for blessing. Not committed, nothing computed. The git
commit of this file (with CinC + Cliff having seen it) is the lock; compute
begins only after.

**Provenance:** CinC brief `BRIEF_paper92_hold_spinstatistics_20260605.md`,
5 June 2026, as redirected by Mr Code (deuteron-primary) and confirmed by
CinC. Paper 92 is HELD pending this verdict.

---

## The single bit

Paper 92 models the neutron as a proton with its 108-bridge occupied by a
shared electron core (P-e-P). The spin-statistics objection that buried the
pre-1932 proton-electron nuclear model collapses, in this framework, to **one
bit**: the spin parity of the bridge-occupying electron core.

Three observed facts test the same bit:

| system | Bootstrap constituents | naive spin-½ count | observed | needs core to be |
|---|---|---|---|---|
| bare neutron | p + core | ½ ⊗ ½ = 0 ⊕ 1 (integer only) | **½** | **non-½** |
| deuteron | p + core + p | three ½ → half-integer | **1** (integer) | **integer (0)** |
| nitrogen-14 | 14 p + 7 core | 21 × ½ → half-integer (fermion) | **boson** (spin 1) | **integer (0)** |

All three are cured iff the core carries **integer spin (spin-0)**; all three
are refuted iff it carries **half-integer spin (½)**. The deuteron (spin-1,
unambiguous, inside Paper 92's own flagship system) is the **primary** target;
the bare neutron is the simpler check beneath it; N-14 is the historical
headline (the count that killed the 1932 model).

## The lemma, stated for the machine

The binary icosahedral group **2I** (order 120, double cover of the icosahedral
group I ≅ A₅, order 60) is exactly the object that distinguishes integer from
half-integer spin: its central element −1 (the 2π rotation) acts as +1 on
single-valued (integer-spin) irreps and as −1 on spinorial (half-integer)
irreps. 2I has 9 irreps:

- **Integer / single-valued** (factor through A₅): dims **1, 3, 3′, 4, 5** (Σ² = 60)
- **Half-integer / spinorial** (faithful, −1 acts as −1): dims **2, 2′, 4′, 6** (Σ² = 60)

The free electron is a Dirac spinor → the fundamental spinorial irrep **2**
(j = ½). **The lemma:** does the *bridge-occupied* electron core carry an
integer (parity-even) or half-integer (parity-odd) 2I representation?

## The selection rule (DERIVED backbone — no modelling freedom)

The central element −1 acts as a scalar ±1 on every irrep; on a composite it
acts as the **product** of the factors. Therefore the integer/spinorial bit of
any composite is exactly its **spinor-count parity**: −1 ↦ (−1)^(number of
spinor constituents). Two consequences, both proved in step 2, that govern the
whole lemma:

- **Construction-independence.** Induction from any subgroup H ∋ −1 conserves
  the central character exactly (−1 acts as the same scalar on `Ind^{2I}_H(ρ)`
  as on ρ): induction *cannot* mix integer and spinorial blocks. Clebsch–Gordan
  composition obeys the same product rule. So **induced-from-stabiliser and
  direct angular-momentum composition provably agree on the bit** — the
  construction choice does not and cannot decide the verdict.
- **The bit is set upstream, by one input:** the spinor parity of the bound
  core. Given it, the composites are forced:

  | composite | constituents | parity ⇒ spin | observed | ⇒ core parity must be |
  |---|---|---|---|---|
  | neutron | p + core | (odd + c) | ½ (odd) | **even (integer)** |
  | deuteron | p + core + p | (even + c) | 1 (even) | **even (integer)** |
  | N-14 | 14p + 7 core | (even + 7c) | boson (even) | **even (integer)** |

  All three give the **same** necessary condition — a clean DERIVED *iff*:
  the P-e-P model is spin-statistics-consistent **iff the bound core is
  integer-spin (parity-even).** This much is rigorous and construction-free.

## The hinge (the one physical claim 2I rep theory frames but cannot settle)

The free electron is the spinor **2** (parity-odd). The model survives **iff**
the bridge occupation renders the bound core parity-**even**. Rep theory does
not force this — it isolates it to a single sharp question:

> Does the 108-bridge standing-wave occupation delete the electron's intrinsic
> spinor parity, or does the core retain spin-½ as any bound electron does?

A single-valued *spatial* standing wave (a function on the icosahedral
geometry) is integer — but the electron's intrinsic spin-½ is a separate factor
that spatial confinement does not remove. So outcome 1 requires a **geometric
absorption/pairing mechanism** for the spinor parity, not merely a single-valued
spatial mode. That mechanism either exists in the 2I geometry or it does not.

## The named forcing mechanism to test (pre-committed, CinC 5 June)

The candidate mechanism is **neutrino ejection — recovered from Paper 92's own
formation story, not invented for this lemma.** Paper 92 §3.1 states
independently (for energetic reasons: the core contracts from Bohr to Compton
scale): *"The extended electron field — the neutrino — is ejected at
formation."* The hypothesis to derive:

> The electron core's spinor (parity-odd) factor is carried out by the ejected
> spin-½ neutrino, leaving a parity-even bound remnant. By fermion-number
> conservation, electron(odd) → neutrino(odd) ⊗ remnant ⟹ remnant is
> parity-even. Test whether this makes the bound count parity-even
> **consistently** across neutron, deuteron, and N-14 — count the fermions
> after ejection.

This is a *requirement-type* mechanism (the ejection is in the model for
independent energetic reasons + beta-decay phenomenology), which is the bar:
the mechanism is not tuned to dodge N-14.

## Hard rule: permission ≠ requirement (pre-committed, CinC 5 June)

The load-bearing anti-fiat clause, now pointed at the real hinge:

- *"The spatial standing-wave mode is single-valued"* only **permits**
  parity-even; it does **not force** it (a bound electron keeps intrinsic
  spin-½ however its spatial mode is confined).
- The mechanism counts as **forcing** only if the geometry/conservation
  **requires** the parity flip — the neutrino ejection actually removing one
  spin-½ factor from the bound count, **derived not asserted**.
- If the strongest defensible statement at the end is "single-valued spatial
  mode permits it," that is **OUTCOME 3 (conditional), not OUTCOME 1.** The
  convenient parity is never assigned to reach the exciting answer.

## Method — FIXED BEFORE THE ANSWER IS KNOWN

This is the anti-fiat commitment. No step may be adjusted after a verdict is
glimpsed.

1. **Build 2I** explicitly as 120 unit quaternions over ℤ[φ]. Reuse
   `golden-dirac/`; verify order 120 and class structure independently.
2. **Compute the full 2I character table** (9 irreps, exact over ℚ(√5)). Verify
   orthonormality. Classify each irrep by χ(−1)/χ(1) (= +1 integer, = −1
   spinorial); cross-check the integer block against the A₅ table in
   `paper-157-rigour/`. **Prove the two backbone facts above** (central-character
   conservation under induction; the parity product rule).
3. **Compute the composite bookkeeping under BOTH constructions (mandatory
   cross-check, named in advance):**
   - **(A) Induced-from-stabiliser.** Identify the stabiliser H ≤ 2I of the
     occupied bridge axis; form `Ind^{2I}_H` of the core's local state.
     *Defence on the record:* this is the standard build of a localised
     excitation's global transformation from its site symmetry — the natural
     reading of "a mode pinned to one bridge," chosen before the answer is
     known, not because it lands in any particular block (and by the
     conservation theorem it provably cannot bias the block).
   - **(B) Direct angular-momentum (Clebsch–Gordan) composition.** Compose the
     constituent spins under SU(2) ↓ 2I and decompose — the reading closest to
     how angular momenta physically combine and to the observed-spin check.
   By the backbone, (A) and (B) must agree on the parity bit; **computing both
   is the audit that the bookkeeping was done right**, and any disagreement is
   a red flag of an error to be found, not a result.
4. **Determine the bound-core parity — the hinge — via the named mechanism.**
   Test neutrino ejection: count the bound fermions before and after the
   formation ejection across all three composites. The question is whether the
   model's standing ejection commitment (one spin-½ neutrino per core) **forces**
   the remnant parity-even via fermion-number conservation, or only permits it.
   Apply the hard rule: forcing requires the spin-½ factor genuinely leaving the
   bound count, derived; "single-valued spatial mode permits it" is not enough.
   Flag the determination's status explicitly (DERIVED only if conservation
   forces it given the independently-motivated ejection; CONJECTURED if it rests
   on an unforced choice).
5. **Verdict rule (pre-committed):**
   - geometry **forces** the bound core parity-even, surviving the "where does
     the spin-½ go" question → **OUTCOME 1**: all three cured, strongest result.
   - bound core **necessarily** retains odd (spinor) parity → **OUTCOME 2**:
     P-e-P neutron refuted, logged clean (the 150 move).
   - rep theory proves the *iff* but the bound-core parity is **not forced**
     without an unjustified modelling choice → **OUTCOME 3**: the lemma stands
     as "92 is consistent iff the core is single-valued," that claim flagged
     CONJECTURED, 92 held as conditional-on-it. **The convenient parity is
     never assigned to reach outcome 1** — absent a forcing mechanism, the
     honest verdict is 3, not 1.

## What would make the result wrong (stated up front)

- Choosing the stabiliser / induction in step 3 *after* seeing which choice
  gives integer spin. (Forbidden; step 3 is fixed here.)
- Reporting "dominantly integer" when the decomposition is mixed, without
  reporting the spinorial weight. (Mixed → report the weights → likely
  OUTCOME 3.)
- Silently contradicting the published vertex↔proton / face↔electron
  correspondence (α paper, Paper 127). The derivation must be consistent with
  it; if it forces a contradiction, that is a flag event, reported, not buried.

## Outputs

- `compute_2I_irreps.py` — builds 2I, character table, classification (steps 1–2).
- `bridge_mode_decomposition.py` — step 3–4 induced rep + decomposition.
- `findings.md` — verdict (1/2/3), the decomposition, sensitivity to the
  step-3 construction, status flags (DERIVED for the group theory; the
  physical mode identification flagged honestly).
- Verdict scoped to exactly what the lemma earns — no claim about the neutron
  beyond the bit decided here.

## Deferred (separate task — Q4 answered)

The degeneration-units redo is **not** in this pre-reg; it is independent of
the spin-statistics bit. Q4 (5 June): no prior compute exists (pre-GitHub /
pre-Mr-Code), so it is **fresh work under W-104 dimensional discipline**, not a
repair. **Blocker:** the published paper contains no degeneration derivation —
formation and degeneration were uncommitted live-compute. "Fresh" therefore
needs a target scoped by CinC: what the degeneration derivation should
**produce** (input → output → expected units). Candidates: binding energy
2.224 MeV derived dimensionally (currently OBSERVED, only the ratio mπ/Ebind ≈
63 is derived); formation/neutrino-ejection energetics; or something not in the
paper. No pre-reg needed once scoped — it is a rigour redo, not a null test.

🐕☕⬡
