# Paper 157 — Arithmetic Rein-Back (audit finding)

**Independent of the block-weight test.** Brief v0.2 §5. This is a
Pattern 39 / Pattern 75 headline error in a live paper (v1.3) and is
actioned at 157's own rein-back by Cliff/CinC. Mr Code logs it; does not
edit the paper.

## The error

Paper 157 v1.3 states the value of its own formula as **137.06** with
residual **0.02%**. Verified at 40 dp (mpmath, `scoping_rep_theory.py`):

```
L(1,χ₅) = 2 ln φ / √5       = 0.43040894096400404   (paper's 0.43041 is correct)
α⁻¹ = 59√5/(2 ln φ)         = 137.07893676152580    (NOT 137.06)
CODATA 2022                 = 137.035999177
residual                    = 0.04293758  →  0.0313%  (NOT 0.02%)
```

The 0.02% figure is internally consistent with 137.06 having been used as
the base (|137.06−137.036|/137.036 ≈ 0.018%), so the slip originates in
evaluating 59√5/(2 log φ) and propagates into the residual. The sign is
unchanged (overshoot). Net effect: **the formula is ~1.6× further from
CODATA than v1.3 claims.** OBSERVED status is unchanged and honest — the
number simply must be right.

## Occurrences to correct (137.06 → 137.0789, 0.02% → 0.031%)

1. **Abstract** — "α⁻¹ = … = 137.06, matching … to 0.02%".
2. **§4, eq. (3)** — "= 59√5 / (2 log φ) = 137.06"; and "The discrepancy
   is 0.02%."
3. **Table 3** — row `59/L(1,χ₅) — this paper, eq. (3) | 137.06 | 0.02% |
   OBSERVED`.
4. **§4.2** — the "0.02% gap" heading and "discrepancy of 0.02% between
   137.06 and 137.036".
5. **§8 Conclusion** — "This gives α⁻¹ = 137.06, within 0.02% of the
   CODATA 2022 value."

Suggested corrected phrasings: value **137.0789** (or 137.079 at 3 dp),
residual **0.031%**. The §4.2 "consistent with higher corrections"
sentence may be retained, but note the v0.1/v0.2 tests found the tower-
residual and block-weight routes to that correction are **not forced**
(see `GATE_FINDINGS.md`), so it should read as an open question, not a
near-term mechanism.
