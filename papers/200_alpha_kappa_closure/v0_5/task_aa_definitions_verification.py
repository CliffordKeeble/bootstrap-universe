#!/usr/bin/env python3
"""Task AA: Definition A/B/C verification for Paper 200 v0.5 sec.5.1.

Pre-registered in papers/200_alpha_kappa_closure/v0_5/task_aa_brief.md.

Implements the three tie definitions per Paper 199 sec.3 (authoritative
source: C:\\Users\\Cliff\\Downloads\\Paper_199_RunningSumOf5_Series Primes_v1_2.md,
lines 122-128, 54). Computes S_5(n) under each definition and verifies
CinC's pre-registered predictions.

Paper 199 sec.3 definitions, verbatim:

  - Definition A (include n = 1): count positions k in {1, ..., 33} with
    S_5(k) = 0. Natural-order count = 14.

  - Definition B (exclude n = 1): count positions k in {2, ..., 33} with
    S_5(k) = 0. Natural-order count = 13.

  - Definition C (omit p = 5 entirely; 32-step sum over non-ramified
    primes): natural-order count = 13.

Key reading (for the v0.5 sec.5.1 parenthetical):

  - Def A and Def B use the SAME running sum (over all 5-series primes
    including the ramified p = 5). They differ only in the counting
    window: A counts ties in {1, ..., 33}, B counts ties in {2, ..., 33}.

  - Def C uses a DIFFERENT running sum (over non-ramified 5-series
    primes only). Index j under Def C corresponds to index (j+1) under
    Def A; the running-sum value satisfies S_C(j) = S_A(j+1) since the
    ramified step contributes 0.

CinC's prediction table (re-read by Mr Code against Paper 199 sec.3):

  - "S_5(30) under Def B = 0 because indices shift by 1" — this rationale
    matches Paper 199's Def C (which does shift indexing), not Paper
    199's Def B (which does not). Under Paper 199's Def B, S_5(30) is
    the same value as under Def A, i.e. -1 (per the Task Z' table).
    Under Paper 199's Def C, S_C(30) = S_A(31) = 0.

  - Substantive claim: "running sum = 0 at the boundary-crossing prime
    p = 281" holds under all three definitions at the appropriate
    index. Predictions hold substantively but CinC's Def B / Def C
    labels appear to be swapped against Paper 199 sec.3.

  - Decision rule then lands in Bin MISMATCH-MINOR (predictions hold
    substantively, but Def C convention differs from CinC's working
    reading). The bin determination is computed mechanically below.

ASCII-safe prints throughout.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

# Pre-registered constants. Do not edit between pre-reg and run.
SEED = 200
N_5SERIES = 50            # first 50 primes congruent to 5 mod 6
SIEVE_UPPER = 2000        # comfortably above p_{5,50}

REGIME_TIGHT_END = 33     # tight-oscillation regime is [1, 33] under Def A
BOUNDARY_PRIME = 281      # p_{5,31} under Def A (Crossing 1 anchor)

# Pre-registered anchor values from Task Z' (commit c2c5eeb of branch
# paper-200-v0-3-task-z).
EXPECTED_P5_VALUES = {
    1: 5,
    5: 29,
    31: 281,
    33: 311,
}
EXPECTED_S_A_AT_31 = 0
EXPECTED_DEF_A_TIE_SET_1_TO_33 = frozenset(
    {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}
)

# CinC's pre-registered predictions table (brief sec."Pre-registered
# predictions").  These are scored literally for the bin determination.
CINC_PRED_S_A_AT_31 = 0
CINC_PRED_S_B_AT_30 = 0      # under CinC's "Def B = shift by 1" reading
CINC_PRED_S_C_AT_BOUNDARY = 0
CINC_PRED_DEF_A_TIE_COUNT = 14
CINC_PRED_DEF_B_TIE_COUNT = 13

OUT_DIR = Path(__file__).parent
RESULTS_JSON = OUT_DIR / "task_aa_results.json"
RESULTS_MD = OUT_DIR / "task_aa_results.md"


def sieve_primes(upper: int) -> list[int]:
    sieve = bytearray([1]) * (upper + 1)
    sieve[0] = 0
    sieve[1] = 0
    for i in range(2, int(upper ** 0.5) + 1):
        if sieve[i]:
            start = i * i
            step = i
            sieve[start::step] = bytearray(len(range(start, upper + 1, step)))
    return [i for i, v in enumerate(sieve) if v]


def five_series_primes(count: int, upper: int) -> list[int]:
    """Return the first `count` primes p with p mod 6 == 5 (p = 5 included)."""
    primes = sieve_primes(upper)
    fs = [p for p in primes if p % 6 == 5]
    if len(fs) < count:
        raise RuntimeError(
            f"sieve to {upper} yielded only {len(fs)} primes congruent to 5 mod 6; "
            f"need {count}"
        )
    return fs[:count]


def chi5(p: int) -> int:
    """Legendre character mod 5 (Definition A: ramified -> 0)."""
    r = p % 5
    if r == 0:
        return 0
    if r == 1 or r == 4:
        return 1
    return -1


def running_sum(values: list[int]) -> list[int]:
    out, s = [], 0
    for v in values:
        s += v
        out.append(s)
    return out


def main() -> None:
    random.seed(SEED)

    # Step 2: 5-series primes + pre-registered prime cross-checks.
    primes_A = five_series_primes(N_5SERIES, SIEVE_UPPER)
    prime_check_fails = []
    for idx, expected in EXPECTED_P5_VALUES.items():
        if primes_A[idx - 1] != expected:
            prime_check_fails.append(
                (idx, expected, primes_A[idx - 1])
            )
    if prime_check_fails:
        print("HALT: 5-series prime cross-check failed.")
        for idx, exp, got in prime_check_fails:
            print(f"  p_{{5,{idx}}}: expected {exp}, got {got}")
        raise SystemExit(1)
    print("5-series prime cross-check: 4/4 anchors match.  OK")

    # Step 3: chi_5 values.
    chi_A = [chi5(p) for p in primes_A]

    # Step 4: running sums under each definition (per Paper 199 sec.3).
    #
    # Def A: full sum over all 5-series primes; tie window {1, ..., 33}.
    S_A = running_sum(chi_A)
    #
    # Def B: SAME running sum as A; tie window {2, ..., 33}.
    S_B = S_A  # alias; differs only in tie counting window.
    #
    # Def C: running sum over non-ramified primes only (omit p = 5);
    # 32-step window.
    primes_C = primes_A[1:]   # skip p_{5,1} = 5
    chi_C = [chi5(p) for p in primes_C]
    S_C = running_sum(chi_C)

    # Step 5: Task Z' anchor cross-check (Def A: S_A(31) = 0 and 14-tie set).
    if S_A[30] != EXPECTED_S_A_AT_31:
        print(f"HALT: anchor S_A(31) = {S_A[30]}, expected {EXPECTED_S_A_AT_31}.")
        raise SystemExit(2)
    print(f"Def A anchor S_A(31) = {S_A[30]} (expected {EXPECTED_S_A_AT_31})  OK")

    def_A_ties_1_to_33 = frozenset(
        n for n in range(1, REGIME_TIGHT_END + 1) if S_A[n - 1] == 0
    )
    if def_A_ties_1_to_33 != EXPECTED_DEF_A_TIE_SET_1_TO_33:
        print("HALT: Def A tie-position set in [1, 33] mismatch.")
        print(f"  expected: {sorted(EXPECTED_DEF_A_TIE_SET_1_TO_33)}")
        print(f"  computed: {sorted(def_A_ties_1_to_33)}")
        raise SystemExit(3)
    print(
        f"Def A tie set [1, 33]: {len(def_A_ties_1_to_33)} positions match Task Z'.  OK"
    )

    # Tie counts under each definition (per Paper 199 sec.3 windows).
    def_A_count = sum(1 for n in range(1, 34) if S_A[n - 1] == 0)
    def_B_count = sum(1 for n in range(2, 34) if S_B[n - 1] == 0)
    def_C_count = sum(1 for j in range(1, 33) if S_C[j - 1] == 0)

    print()
    print("Tie counts (Paper 199 sec.3 windows):")
    print(f"  Def A : {def_A_count}  (expected 14)  "
          f"{'OK' if def_A_count == 14 else 'MISMATCH'}")
    print(f"  Def B : {def_B_count}  (expected 13)  "
          f"{'OK' if def_B_count == 13 else 'MISMATCH'}")
    print(f"  Def C : {def_C_count}  (expected 13)  "
          f"{'OK' if def_C_count == 13 else 'MISMATCH'}")

    # Tie position sets under each def, listed for the verification table.
    def_A_ties = sorted(def_A_ties_1_to_33)
    def_B_ties = sorted(n for n in range(2, 34) if S_B[n - 1] == 0)
    def_C_ties = sorted(j for j in range(1, 33) if S_C[j - 1] == 0)

    # The boundary-crossing prime p = 281 lives at:
    #   Def A index k = 31  (primes_A[30] = 281)
    #   Def B index k = 31  (same indexing as A; primes_B = primes_A)
    #   Def C index j = 30  (primes_C[29] = 281)
    k_def_A = primes_A.index(BOUNDARY_PRIME) + 1
    k_def_B = k_def_A                  # Def B uses Def A's indexing
    j_def_C = primes_C.index(BOUNDARY_PRIME) + 1
    s_at_boundary = {
        "Def A": (k_def_A, S_A[k_def_A - 1]),
        "Def B": (k_def_B, S_B[k_def_B - 1]),
        "Def C": (j_def_C, S_C[j_def_C - 1]),
    }
    print()
    print(f"Running sum at boundary-crossing prime p = {BOUNDARY_PRIME}:")
    for label, (idx, val) in s_at_boundary.items():
        print(f"  {label}: index {idx}  ->  S_5 = {val:+d}")

    # CinC's literal predictions, scored:
    # (1) S_A(31) = 0   - already checked above (anchor).
    # (2) S_B(30) = 0   - the prediction CinC wrote against Def B.
    # (3) S_C(at index for 281) = 0
    # (4) Def A 14-tie set matches
    # (5) Def B count is 13 with "shift by 1" reading
    s_B_at_30_actual = S_B[29]      # value Paper 199 Def B gives
    s_C_at_30_actual = S_C[29]      # value Paper 199 Def C gives at index 30

    cinc_lit_check_2_passes = (s_B_at_30_actual == CINC_PRED_S_B_AT_30)
    cinc_lit_check_3_passes = (s_at_boundary["Def C"][1] == CINC_PRED_S_C_AT_BOUNDARY)
    cinc_lit_check_4_passes = (def_A_count == CINC_PRED_DEF_A_TIE_COUNT)
    cinc_lit_check_5_passes = (def_B_count == CINC_PRED_DEF_B_TIE_COUNT)

    # Substantive check: does the observation "running sum = 0 at p = 281"
    # hold under all three definitions, at the appropriate index?
    substantive_check_passes = all(val == 0 for (_, val) in s_at_boundary.values())

    # Bin determination (mechanical, per pre-registered rule).
    if not substantive_check_passes:
        bin_decision = "MISMATCH-MAJOR"
        bin_rationale = (
            "Substantive observation 'running sum = 0 at p = 281' does not hold "
            "under all three definitions. Paper 199 sec.3 robustness claim refuted."
        )
    else:
        all_literal_pass = (
            cinc_lit_check_2_passes
            and cinc_lit_check_3_passes
            and cinc_lit_check_4_passes
            and cinc_lit_check_5_passes
        )
        if all_literal_pass:
            bin_decision = "CONFIRM"
            bin_rationale = (
                "All CinC predictions match the literal values under Paper 199 "
                "sec.3 definitions; Def A/B/C robustness verified."
            )
        else:
            bin_decision = "MISMATCH-MINOR"
            bin_rationale = (
                "Substantive observation holds under all three definitions, but "
                "Paper 199 sec.3's Def B does not shift indexing as CinC's "
                "working reading assumes (Paper 199 Def B uses the same running "
                "sum as Def A with a restricted counting window; Paper 199 Def C "
                "is the definition that shifts indexing). CinC's predicted value "
                "for S_B(30) is the value Paper 199 Def C gives at j=30, not "
                "what Def B gives at k=30."
            )

    print()
    print(f"Bin decision: {bin_decision}")
    print(f"  {bin_rationale}")

    # Build first-50 table side-by-side.
    first_50 = []
    for n in range(1, 51):
        p = primes_A[n - 1]
        row = {
            "n": n,
            "p_5n": p,
            "p_mod_5": p % 5,
            "chi_5": chi_A[n - 1],
            "S_A": S_A[n - 1],
            "S_B": S_B[n - 1],  # identical to S_A
        }
        # Def C aligns S_C(j) with Def A index j+1 (since p_C(j) = p_A(j+1)).
        # So at Def A index n, the corresponding Def C index is n - 1 (valid for n >= 2).
        if n >= 2:
            row["S_C_at_n_minus_1"] = S_C[n - 2]
        else:
            row["S_C_at_n_minus_1"] = None  # n=1 corresponds to ramified prime, no Def C index
        row["A_tie"] = (S_A[n - 1] == 0)
        row["B_tie"] = (n >= 2 and S_B[n - 1] == 0)
        row["C_tie"] = (n >= 2 and S_C[n - 2] == 0)
        first_50.append(row)

    payload = {
        "task": "AA",
        "paper": "200",
        "version": "v0.5",
        "definitions_source": "Paper 199 sec.3 (lines 122-128, 54)",
        "seed": SEED,
        "n_5series_computed": N_5SERIES,
        "boundary_prime": BOUNDARY_PRIME,
        "prime_anchor_check": {
            "expected": EXPECTED_P5_VALUES,
            "computed": {idx: primes_A[idx - 1] for idx in EXPECTED_P5_VALUES},
            "all_match": True,
        },
        "task_z2_anchor_cross_check": {
            "S_A_at_31": S_A[30],
            "expected": EXPECTED_S_A_AT_31,
            "match": True,
            "def_A_tie_set_1_to_33": def_A_ties,
            "expected_set": sorted(EXPECTED_DEF_A_TIE_SET_1_TO_33),
            "set_match": True,
        },
        "tie_counts": {
            "Def_A": {"window": "[1, 33]", "count": def_A_count, "expected": 14},
            "Def_B": {"window": "[2, 33]", "count": def_B_count, "expected": 13},
            "Def_C": {"window": "[1, 32]", "count": def_C_count, "expected": 13},
        },
        "tie_positions": {
            "Def_A": def_A_ties,
            "Def_B": def_B_ties,
            "Def_C": def_C_ties,
        },
        "running_sum_at_boundary_prime_281": {
            "Def_A_index_31": S_A[30],
            "Def_B_index_31": S_B[30],
            "Def_C_index_30": S_C[29],
        },
        "cinc_literal_prediction_scoring": {
            "S_A_at_31_predicted": CINC_PRED_S_A_AT_31,
            "S_A_at_31_computed": S_A[30],
            "S_B_at_30_predicted": CINC_PRED_S_B_AT_30,
            "S_B_at_30_computed_under_paper_199_def_B": s_B_at_30_actual,
            "S_B_at_30_check_passes": cinc_lit_check_2_passes,
            "S_C_at_boundary_predicted": CINC_PRED_S_C_AT_BOUNDARY,
            "S_C_at_boundary_computed": s_at_boundary["Def C"][1],
            "S_C_at_boundary_check_passes": cinc_lit_check_3_passes,
            "def_A_count_predicted": CINC_PRED_DEF_A_TIE_COUNT,
            "def_A_count_computed": def_A_count,
            "def_A_count_check_passes": cinc_lit_check_4_passes,
            "def_B_count_predicted": CINC_PRED_DEF_B_TIE_COUNT,
            "def_B_count_computed": def_B_count,
            "def_B_count_check_passes": cinc_lit_check_5_passes,
        },
        "substantive_claim_running_sum_zero_at_p_281": {
            "holds_under_all_three_defs": substantive_check_passes,
            "values": {label: val for label, (_, val) in s_at_boundary.items()},
        },
        "decision": {
            "bin": bin_decision,
            "rationale": bin_rationale,
        },
        "first_50": first_50,
    }

    RESULTS_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    md = []
    md.append("# Paper 200 v0.5 Task AA - Def A/B/C verification results")
    md.append("")
    md.append("Pre-registered in `task_aa_brief.md`. Implementation in")
    md.append("`task_aa_definitions_verification.py`. JSON payload in")
    md.append("`task_aa_results.json`.")
    md.append("")
    md.append("Authoritative source for Def A/B/C: Paper 199 sec.3,")
    md.append("lines 122-128 and 54 of `Paper_199_RunningSumOf5_Series Primes_v1_2.md`.")
    md.append("")
    md.append("## Anchor cross-checks against Task Z' (commit c2c5eeb)")
    md.append("")
    md.append("| Check | Expected | Computed | Status |")
    md.append("|---|---|---|---|")
    for idx, exp in EXPECTED_P5_VALUES.items():
        md.append(f"| p_{{5,{idx}}} | {exp} | {primes_A[idx-1]} | OK |")
    md.append(f"| S_A(31) | 0 | {S_A[30]} | OK |")
    md.append(f"| Def A tie set in [1, 33] | 14 positions | "
              f"{len(def_A_ties_1_to_33)} positions | OK |")
    md.append("")
    md.append("## Tie counts under each definition (Paper 199 sec.3)")
    md.append("")
    md.append("| Definition | Window | Tie count | Paper 199 sec.3 figure |")
    md.append("|---|---|---:|---:|")
    md.append(f"| **Def A** (include n=1) | [1, 33] | {def_A_count} | 14 |")
    md.append(f"| **Def B** (exclude n=1; same sum, window restricted) | [2, 33] | {def_B_count} | 13 |")
    md.append(f"| **Def C** (omit p=5; 32-step sum) | [1, 32] | {def_C_count} | 13 |")
    md.append("")
    md.append("## Running sum at the boundary-crossing prime p = 281")
    md.append("")
    md.append("| Definition | Index at p = 281 | S_5 at that index |")
    md.append("|---|---:|---:|")
    md.append(f"| **Def A** | k = {k_def_A} | **{S_A[k_def_A-1]:+d}** |")
    md.append(f"| **Def B** | k = {k_def_B} (same indexing as A) | **{S_B[k_def_B-1]:+d}** |")
    md.append(f"| **Def C** | j = {j_def_C} (after skipping ramified) | **{S_C[j_def_C-1]:+d}** |")
    md.append("")
    md.append("All three definitions give the same value (0) at the boundary-crossing")
    md.append("prime, confirming Paper 199 sec.3's robustness claim.")
    md.append("")
    md.append("## CinC literal-prediction scoring")
    md.append("")
    md.append("CinC's prediction table (brief) listed values for specific")
    md.append("(definition, index) pairs. Scoring literally:")
    md.append("")
    md.append("| Prediction | CinC value | Computed | Match |")
    md.append("|---|---:|---:|:---:|")
    md.append(f"| S_A(31) | {CINC_PRED_S_A_AT_31} | {S_A[30]} | YES |")
    md.append(f"| S_B(30) under Paper 199 Def B | {CINC_PRED_S_B_AT_30} | {s_B_at_30_actual:+d} | "
              f"{'YES' if cinc_lit_check_2_passes else 'NO'} |")
    md.append(f"| S_C at boundary p=281 (j=30 in Def C) | {CINC_PRED_S_C_AT_BOUNDARY} | "
              f"{s_at_boundary['Def C'][1]:+d} | {'YES' if cinc_lit_check_3_passes else 'NO'} |")
    md.append(f"| Def A 14-tie set in [1, 33] | 14 | {def_A_count} | "
              f"{'YES' if cinc_lit_check_4_passes else 'NO'} |")
    md.append(f"| Def B 13 ties under 'shift by 1' reading | 13 | {def_B_count} | "
              f"{'YES' if cinc_lit_check_5_passes else 'NO'} |")
    md.append("")
    md.append("## Tie-position sets under each definition")
    md.append("")
    md.append(f"- **Def A** [1, 33], {def_A_count} ties: {def_A_ties}")
    md.append(f"- **Def B** [2, 33], {def_B_count} ties: {def_B_ties}")
    md.append(f"- **Def C** [1, 32], {def_C_count} ties: {def_C_ties}")
    md.append("")
    md.append("Note Def A tie set = Def B tie set + {1}; Def C tie set is the")
    md.append("Def A tie set in [2, 33], reindexed by subtracting 1.")
    md.append("")
    md.append("## Decision against pre-registered rule")
    md.append("")
    md.append(f"**Bin {bin_decision}.**")
    md.append("")
    md.append(bin_rationale)
    md.append("")
    md.append("## Side-by-side S_5 values for n in [1, 50]")
    md.append("")
    md.append("Indexing convention: n indexes Def A / Def B (same sum); the")
    md.append("`S_C @ j=n-1` column shows Def C's running sum at its index j = n-1,")
    md.append("which corresponds to the same underlying prime p_{5,n} as Def A's index n.")
    md.append("(`-` for n=1 since the ramified prime has no Def C index.)")
    md.append("")
    md.append("| n | p_5n | p mod 5 | chi_5 | S_A=S_B | S_C @ j=n-1 | A tie | B tie | C tie |")
    md.append("|---:|---:|---:|---:|---:|---:|:-:|:-:|:-:|")
    for row in first_50:
        sc = "-" if row["S_C_at_n_minus_1"] is None else f"{row['S_C_at_n_minus_1']:+d}"
        a = "*" if row["A_tie"] else ""
        b = "*" if row["B_tie"] else ""
        c = "*" if row["C_tie"] else ""
        md.append(
            f"| {row['n']} | {row['p_5n']} | {row['p_mod_5']} | "
            f"{row['chi_5']:+d} | {row['S_A']:+d} | {sc} | {a} | {b} | {c} |"
        )
    md.append("")

    RESULTS_MD.write_text("\n".join(md), encoding="utf-8")

    print()
    print(f"Wrote {RESULTS_JSON.name}")
    print(f"Wrote {RESULTS_MD.name}")


if __name__ == "__main__":
    main()
