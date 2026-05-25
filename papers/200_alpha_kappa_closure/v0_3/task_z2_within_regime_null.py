#!/usr/bin/env python3
"""Task Z': within-regime null test for Paper 200 v0.3 (corrected convention).

Pre-registered in papers/200_alpha_kappa_closure/v0_3/task_z2_brief.md.

The 5-series is the sequence of primes congruent to 5 (mod 6). For each
n >= 1, p_{5,n} is the n-th such prime. The running sum is

    S_5(n) = sum_{j=1}^{n} chi_5(p_{5,j})

under Definition A (Paper 199 sec.6.1):
  chi_5(p) = +1  if p mod 5 in {1, 4}
  chi_5(p) = -1  if p mod 5 in {2, 3}
  chi_5(p) =  0  if p mod 5 == 0  (only the ramified prime p = 5)

Note: 5 mod 6 = 5, so the ramified prime IS the first 5-series prime,
and S_5(1) = 0 (counted as a tie at n=1).

Anchor verification table (CinC hand-computation, brief sec."Anchor
verification") is hard-coded below and cross-checked row-by-row. If any
row differs, halt cleanly.

Pre-registered expected tie-position set in n in [1, 33]:
  {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}  (14 positions).
If the computed set differs, halt cleanly.

ASCII-safe prints throughout (Windows console encoding).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

# Pre-registered constants. Do not edit between pre-reg and run.
SEED = 200
N_5SERIES = 200       # number of 5-series primes to compute
SIEVE_UPPER = 5000    # comfortably above p_{5,200}

REGIME_TIGHT = (1, 33)
REGIME_POST = (34, 66)
REGIME_ASYMP = (67, 200)

BIN_LOW_MAX = 0.30
BIN_MID_MAX = 0.55

# Pre-registered anchor-verification table. Tuples: (n, p_5n, p_mod_5, chi5, S5, tie).
# Source: task_z2_brief.md sec."Anchor verification" (CinC hand-computation).
EXPECTED_ANCHOR_TABLE = [
    (1,    5, 0,  0,  0, True),
    (2,   11, 1, +1,  1, False),
    (3,   17, 2, -1,  0, True),
    (4,   23, 3, -1, -1, False),
    (5,   29, 4, +1,  0, True),
    (6,   41, 1, +1,  1, False),
    (7,   47, 2, -1,  0, True),
    (8,   53, 3, -1, -1, False),
    (9,   59, 4, +1,  0, True),
    (10,  71, 1, +1,  1, False),
    (11,  83, 3, -1,  0, True),
    (12,  89, 4, +1,  1, False),
    (13, 101, 1, +1,  2, False),
    (14, 107, 2, -1,  1, False),
    (15, 113, 3, -1,  0, True),
    (16, 131, 1, +1,  1, False),
    (17, 137, 2, -1,  0, True),
    (18, 149, 4, +1,  1, False),
    (19, 167, 2, -1,  0, True),
    (20, 173, 3, -1, -1, False),
    (21, 179, 4, +1,  0, True),
    (22, 191, 1, +1,  1, False),
    (23, 197, 2, -1,  0, True),
    (24, 227, 2, -1, -1, False),
    (25, 233, 3, -1, -2, False),
    (26, 239, 4, +1, -1, False),
    (27, 251, 1, +1,  0, True),
    (28, 257, 2, -1, -1, False),
    (29, 263, 3, -1, -2, False),
    (30, 269, 4, +1, -1, False),
    (31, 281, 1, +1,  0, True),  # ANCHOR
]

# Pre-registered expected tie-position set in [1, 33], from brief.
EXPECTED_TIE_POSITIONS_1_TO_33 = frozenset(
    {1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 27, 31, 33}
)

OUT_DIR = Path(__file__).parent
RESULTS_JSON = OUT_DIR / "task_z2_results.json"
RESULTS_MD = OUT_DIR / "task_z2_results.md"


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
    """Return the first `count` primes p with p mod 6 == 5.

    p = 5 IS included (5 mod 6 = 5). Excludes 2 and 3.
    """
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
    out = []
    s = 0
    for v in values:
        s += v
        out.append(s)
    return out


def fraction_zeros(s_values: list[int], lo: int, hi: int) -> tuple[int, int, float]:
    indices = range(lo, hi + 1)
    length = len(indices)
    zeros = sum(1 for n in indices if s_values[n - 1] == 0)
    return zeros, length, (zeros / length if length else 0.0)


def bin_decision(f: float) -> str:
    if f < BIN_LOW_MAX:
        return "LOW"
    if f <= BIN_MID_MAX:
        return "MID"
    return "HIGH"


def main() -> None:
    random.seed(SEED)

    primes_5 = five_series_primes(N_5SERIES, SIEVE_UPPER)
    chi_values = [chi5(p) for p in primes_5]
    s_values = running_sum(chi_values)

    # Step 2: row-by-row anchor verification against CinC's hand-computed table.
    mismatches = []
    for expected in EXPECTED_ANCHOR_TABLE:
        n, e_p, e_mod5, e_chi, e_s, e_tie = expected
        c_p = primes_5[n - 1]
        c_mod5 = c_p % 5
        c_chi = chi_values[n - 1]
        c_s = s_values[n - 1]
        c_tie = (c_s == 0)
        if (c_p, c_mod5, c_chi, c_s, c_tie) != (e_p, e_mod5, e_chi, e_s, e_tie):
            mismatches.append({
                "n": n,
                "expected": {"p": e_p, "p_mod_5": e_mod5, "chi_5": e_chi,
                             "S_5": e_s, "tie": e_tie},
                "computed": {"p": c_p, "p_mod_5": c_mod5, "chi_5": c_chi,
                             "S_5": c_s, "tie": c_tie},
            })

    if mismatches:
        print("HALT: anchor verification table mismatch.")
        for m in mismatches:
            print(f"  n={m['n']}: expected {m['expected']}; got {m['computed']}")
        raise SystemExit(1)

    print("Anchor table cross-check: 31 rows match CinC's hand-computation.  OK")

    # Step 5: tie-position set in [1, 33] against pre-registered set.
    computed_ties_1_to_33 = frozenset(
        n for n in range(1, 34) if s_values[n - 1] == 0
    )
    if computed_ties_1_to_33 != EXPECTED_TIE_POSITIONS_1_TO_33:
        missing = sorted(EXPECTED_TIE_POSITIONS_1_TO_33 - computed_ties_1_to_33)
        extra = sorted(computed_ties_1_to_33 - EXPECTED_TIE_POSITIONS_1_TO_33)
        print("HALT: tie-position set mismatch in [1, 33].")
        print(f"  missing from computed: {missing}")
        print(f"  extra in computed:     {extra}")
        raise SystemExit(2)

    print(
        f"Tie-position set [1, 33]: {sorted(computed_ties_1_to_33)} "
        f"({len(computed_ties_1_to_33)} ties)  OK"
    )

    # Step 6 & 7: regime counts, fractions, bin decision.
    zeros_tight, len_tight, f_tight = fraction_zeros(s_values, *REGIME_TIGHT)
    zeros_post, len_post, f_post = fraction_zeros(s_values, *REGIME_POST)
    zeros_asymp, len_asymp, f_asymp = fraction_zeros(s_values, *REGIME_ASYMP)
    decision = bin_decision(f_tight)

    print()
    print("Regime         range        zeros / length    fraction")
    print(
        f"tight       [{REGIME_TIGHT[0]:>3}, {REGIME_TIGHT[1]:>3}]    "
        f"{zeros_tight:>5} / {len_tight:<5}    {f_tight:.4f}"
    )
    print(
        f"post-regime [{REGIME_POST[0]:>3}, {REGIME_POST[1]:>3}]    "
        f"{zeros_post:>5} / {len_post:<5}    {f_post:.4f}"
    )
    print(
        f"asymptotic  [{REGIME_ASYMP[0]:>3}, {REGIME_ASYMP[1]:>3}]    "
        f"{zeros_asymp:>5} / {len_asymp:<5}    {f_asymp:.4f}"
    )
    print()
    print("Pre-registered bin boundaries:")
    print(f"  LOW : f_tight <  {BIN_LOW_MAX:.2f}")
    print(f"  MID : {BIN_LOW_MAX:.2f} <= f_tight <= {BIN_MID_MAX:.2f}")
    print(f"  HIGH: f_tight >  {BIN_MID_MAX:.2f}")
    print(f"Decision: f_tight = {f_tight:.4f}  ->  Bin {decision}")

    # First-50 table for results.md.
    first_50 = []
    for n in range(1, 51):
        p = primes_5[n - 1]
        first_50.append({
            "n": n,
            "p_5n": p,
            "p_mod_5": p % 5,
            "chi_5": chi_values[n - 1],
            "S_5": s_values[n - 1],
            "tie": s_values[n - 1] == 0,
        })

    zero_indices_full = [n for n in range(1, N_5SERIES + 1) if s_values[n - 1] == 0]
    zero_indices_post = [n for n in zero_indices_full if REGIME_POST[0] <= n <= REGIME_POST[1]]
    zero_indices_asymp = [n for n in zero_indices_full if REGIME_ASYMP[0] <= n <= REGIME_ASYMP[1]]

    payload = {
        "task": "Z'",
        "paper": "200",
        "version": "v0.3",
        "supersedes": "Task Z (commit 61aff82, halted on anchor under wrong-convention brief)",
        "definition": "A (5-series: primes p with p mod 6 == 5; chi_5 Legendre mod 5; chi_5(5) = 0)",
        "seed": SEED,
        "n_5series_computed": N_5SERIES,
        "anchor": {
            "row_count": len(EXPECTED_ANCHOR_TABLE),
            "all_rows_match": True,
            "S_5_at_anchor_index_31": s_values[30],
            "expected": 0,
        },
        "tie_positions_1_to_33": {
            "expected": sorted(EXPECTED_TIE_POSITIONS_1_TO_33),
            "computed": sorted(computed_ties_1_to_33),
            "match": True,
            "count": len(computed_ties_1_to_33),
        },
        "regimes": {
            "tight": {
                "range": list(REGIME_TIGHT),
                "zeros": zeros_tight,
                "length": len_tight,
                "fraction": f_tight,
            },
            "post": {
                "range": list(REGIME_POST),
                "zeros": zeros_post,
                "length": len_post,
                "fraction": f_post,
                "zero_positions": zero_indices_post,
            },
            "asymptotic": {
                "range": list(REGIME_ASYMP),
                "zeros": zeros_asymp,
                "length": len_asymp,
                "fraction": f_asymp,
                "zero_positions": zero_indices_asymp,
            },
        },
        "bin_boundaries": {
            "LOW_max_exclusive": BIN_LOW_MAX,
            "MID_max_inclusive": BIN_MID_MAX,
            "rule": (
                "LOW if f_tight < 0.30; "
                "MID if 0.30 <= f_tight <= 0.55; "
                "HIGH if f_tight > 0.55"
            ),
        },
        "decision": {
            "f_tight": f_tight,
            "bin": decision,
        },
        "zero_indices_in_1_to_200": zero_indices_full,
        "first_50": first_50,
    }

    RESULTS_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    md = []
    md.append("# Paper 200 v0.3 Task Z' - within-regime null test results")
    md.append("")
    md.append("Pre-registered in `task_z2_brief.md`. Implementation in")
    md.append("`task_z2_within_regime_null.py`. JSON payload in `task_z2_results.json`.")
    md.append("Supersedes Task Z (commit `61aff82`, halted on wrong-convention anchor).")
    md.append("")
    md.append(f"- seed = `{SEED}`, N (5-series primes) = `{N_5SERIES}`")
    md.append("- convention: Definition A (5-series = primes congruent to 5 mod 6;")
    md.append("  chi_5 Legendre mod 5; chi_5(5) = 0, ramified prime included)")
    md.append("")
    md.append("## Anchor verification (row-by-row)")
    md.append("")
    md.append(f"All {len(EXPECTED_ANCHOR_TABLE)} rows from CinC's hand-computation match the")
    md.append(f"locked-script output. `S_5(31) = {s_values[30]}` (expected 0) -> **OK**.")
    md.append("")
    md.append("## Tie-position cross-check in n in [1, 33]")
    md.append("")
    md.append(f"- Pre-registered expected set: {sorted(EXPECTED_TIE_POSITIONS_1_TO_33)}")
    md.append(f"- Computed set              : {sorted(computed_ties_1_to_33)}")
    md.append(f"- Match: **True** ({len(computed_ties_1_to_33)} ties, matching Paper 199 sec.6.1)")
    md.append("")
    md.append("## Regime fractions")
    md.append("")
    md.append("| Regime | n range | zeros | length | fraction |")
    md.append("|---|---|---:|---:|---:|")
    md.append(f"| **tight** | [{REGIME_TIGHT[0]}, {REGIME_TIGHT[1]}] | {zeros_tight} | {len_tight} | **{f_tight:.4f}** |")
    md.append(f"| post-regime | [{REGIME_POST[0]}, {REGIME_POST[1]}] | {zeros_post} | {len_post} | {f_post:.4f} |")
    md.append(f"| asymptotic | [{REGIME_ASYMP[0]}, {REGIME_ASYMP[1]}] | {zeros_asymp} | {len_asymp} | {f_asymp:.4f} |")
    md.append("")
    md.append("## Decision against pre-registered rule")
    md.append("")
    md.append(f"- LOW  : f_tight < {BIN_LOW_MAX:.2f}")
    md.append(f"- MID  : {BIN_LOW_MAX:.2f} <= f_tight <= {BIN_MID_MAX:.2f}")
    md.append(f"- HIGH : f_tight > {BIN_MID_MAX:.2f}")
    md.append("")
    md.append(f"f_tight = **{f_tight:.4f}** -> **Bin {decision}**")
    md.append("")
    md.append(f"## Zero indices in [1, {N_5SERIES}]")
    md.append("")
    md.append(", ".join(str(n) for n in zero_indices_full) if zero_indices_full else "(none)")
    md.append("")
    md.append("## First 50 values of S_5(n) (5-series)")
    md.append("")
    md.append("`*` marks rows where S_5(n) = 0.")
    md.append("")
    md.append("| n | p_5n | p mod 5 | chi_5 | S_5(n) | tie |")
    md.append("|---:|---:|---:|---:|---:|:---:|")
    for row in first_50:
        flag = "*" if row["tie"] else ""
        md.append(
            f"| {row['n']} | {row['p_5n']} | {row['p_mod_5']} | "
            f"{row['chi_5']:+d} | {row['S_5']:+d} | {flag} |"
        )
    md.append("")

    RESULTS_MD.write_text("\n".join(md), encoding="utf-8")

    print()
    print(f"Wrote {RESULTS_JSON.name}")
    print(f"Wrote {RESULTS_MD.name}")


if __name__ == "__main__":
    main()
