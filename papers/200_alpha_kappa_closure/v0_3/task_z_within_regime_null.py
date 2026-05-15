#!/usr/bin/env python3
"""Task Z: within-regime null test for Paper 200 v0.3.

Pre-registered in papers/200_alpha_kappa_closure/v0_3/task_brief.md.

Computes chi_5(p_n) for the first 200 primes under Definition A (Paper 199
convention: ramified prime p = 5 contributes chi_5(5) = 0 to the running
sum; running sum indexed by prime ordinal n with p_1 = 2).

Counts the fraction of S_5(n) = 0 events in three pre-registered regimes:

  - Tight        n in [1, 33]
  - Post-regime  n in [34, 66]
  - Asymptotic   n in [67, 200]

Applies the pre-registered decision rule to f_tight:

  - Bin LOW   if f_tight <  0.30
  - Bin MID   if 0.30 <= f_tight <= 0.55
  - Bin HIGH  if f_tight >  0.55

Verifies the structural anchor S_5(31) = 0 before reporting; halts if it
does not reproduce.

ASCII-safe prints throughout (Windows console encoding).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

# Pre-registered constants. Do not edit between pre-reg and run.
SEED = 200
N_PRIMES = 200
ANCHOR_INDEX = 31  # S_5(31) must equal 0 (Paper 199 / Paper 200 sec.5)

REGIME_TIGHT = (1, 33)
REGIME_POST = (34, 66)
REGIME_ASYMP = (67, 200)

BIN_LOW_MAX = 0.30   # f_tight <  0.30  -> LOW
BIN_MID_MAX = 0.55   # 0.30 <= f_tight <= 0.55 -> MID; > 0.55 -> HIGH

OUT_DIR = Path(__file__).parent
RESULTS_JSON = OUT_DIR / "task_z_results.json"
RESULTS_MD = OUT_DIR / "task_z_results.md"


def primes_up_to_count(count: int) -> list[int]:
    """Return the first `count` primes by simple sieve over a generous bound.

    For count = 200, the 200th prime is 1223; sieve to 2000 is plenty.
    """
    if count <= 0:
        return []
    # Bound chosen with a comfortable margin above prime(200) = 1223.
    upper = max(2000, count * 12)
    sieve = bytearray([1]) * (upper + 1)
    sieve[0] = 0
    sieve[1] = 0
    for i in range(2, int(upper ** 0.5) + 1):
        if sieve[i]:
            step = i
            start = i * i
            sieve[start::step] = bytearray(len(range(start, upper + 1, step)))
    primes = [i for i, v in enumerate(sieve) if v]
    if len(primes) < count:
        raise RuntimeError(
            f"sieve upper bound {upper} produced only {len(primes)} primes; "
            f"need {count}"
        )
    return primes[:count]


def chi5(p: int) -> int:
    """Dirichlet character chi_5 (Definition A: ramified p = 5 -> 0)."""
    r = p % 5
    if r == 0:
        return 0
    if r == 1 or r == 4:
        return 1
    return -1  # r in {2, 3}


def running_sum(values: list[int]) -> list[int]:
    out = []
    s = 0
    for v in values:
        s += v
        out.append(s)
    return out


def fraction_zeros(s_values: list[int], lo: int, hi: int) -> tuple[int, int, float]:
    """Return (zero_count, regime_length, fraction) for n in [lo, hi] inclusive.

    s_values is 1-indexed conceptually: s_values[k] corresponds to S_5(k+1).
    We translate to 1-based indexing here.
    """
    indices = range(lo, hi + 1)
    length = len(indices)
    zeros = sum(1 for n in indices if s_values[n - 1] == 0)
    frac = zeros / length if length else 0.0
    return zeros, length, frac


def bin_decision(f_tight: float) -> str:
    if f_tight < BIN_LOW_MAX:
        return "LOW"
    if f_tight <= BIN_MID_MAX:
        return "MID"
    return "HIGH"


def main() -> None:
    # Seed is documented in script header; this RNG is not used in the
    # computation (deterministic). Seeding here for paper-trail consistency.
    random.seed(SEED)

    primes = primes_up_to_count(N_PRIMES)

    # Brief-mandated cross-checks against a standard prime table.
    assert primes[0] == 2, f"p_1 should be 2, got {primes[0]}"
    assert primes[4] == 11, f"p_5 should be 11, got {primes[4]}"
    assert primes[30] == 127, f"p_31 should be 127, got {primes[30]}"

    chi_values = [chi5(p) for p in primes]
    s_values = running_sum(chi_values)

    # Anchor verification: S_5(31) must equal 0.
    s31 = s_values[ANCHOR_INDEX - 1]
    if s31 != 0:
        print(f"HALT: structural anchor failed. S_5({ANCHOR_INDEX}) = {s31}, expected 0.")
        raise SystemExit(1)

    zeros_tight, len_tight, f_tight = fraction_zeros(s_values, *REGIME_TIGHT)
    zeros_post, len_post, f_post = fraction_zeros(s_values, *REGIME_POST)
    zeros_asymp, len_asymp, f_asymp = fraction_zeros(s_values, *REGIME_ASYMP)

    decision = bin_decision(f_tight)

    # Cross-check vs Paper 199 sec.6.1 "14 ties in 33 primes" figure.
    # Brief reads ties_in_[1,33] = 14 -> ratio 14/33 approx 0.4242.
    p199_expected_zeros_in_1_33 = 14
    zeros_in_1_to_33 = sum(1 for n in range(1, 34) if s_values[n - 1] == 0)
    p199_match = (zeros_in_1_to_33 == p199_expected_zeros_in_1_33)

    # Console summary (ASCII safe).
    print("Task Z - Paper 200 v0.3 within-regime null test")
    print(f"seed = {SEED}, n_primes = {N_PRIMES}")
    print(f"anchor: S_5({ANCHOR_INDEX}) = {s31} (expected 0)  OK")
    print()
    print("Regime          range        zeros / length    fraction")
    print(f"tight       [{REGIME_TIGHT[0]:>3}, {REGIME_TIGHT[1]:>3}]    "
          f"{zeros_tight:>5} / {len_tight:<5}    {f_tight:.4f}")
    print(f"post-regime [{REGIME_POST[0]:>3}, {REGIME_POST[1]:>3}]    "
          f"{zeros_post:>5} / {len_post:<5}    {f_post:.4f}")
    print(f"asymptotic  [{REGIME_ASYMP[0]:>3}, {REGIME_ASYMP[1]:>3}]    "
          f"{zeros_asymp:>5} / {len_asymp:<5}    {f_asymp:.4f}")
    print()
    print(f"Pre-registered bin boundaries:")
    print(f"  LOW : f_tight <  {BIN_LOW_MAX:.2f}")
    print(f"  MID : {BIN_LOW_MAX:.2f} <= f_tight <= {BIN_MID_MAX:.2f}")
    print(f"  HIGH: f_tight >  {BIN_MID_MAX:.2f}")
    print(f"Decision: f_tight = {f_tight:.4f}  ->  Bin {decision}")
    print()
    print(f"Paper 199 sec.6.1 cross-check: zeros in [1, 33] = {zeros_in_1_to_33} "
          f"(Paper 199 figure: {p199_expected_zeros_in_1_33})  "
          f"{'MATCH' if p199_match else 'MISMATCH - flag for CinC'}")

    # Build first-50 table for results.md.
    table_rows = []
    for n in range(1, 51):
        p = primes[n - 1]
        chi = chi_values[n - 1]
        s = s_values[n - 1]
        is_zero = "*" if s == 0 else ""
        table_rows.append({
            "n": n,
            "p_n": p,
            "p_n_mod_5": p % 5,
            "chi_5": chi,
            "S_5": s,
            "is_zero": s == 0,
        })

    # Indices where S_5 = 0 anywhere in [1, 200].
    zero_indices = [n for n in range(1, N_PRIMES + 1) if s_values[n - 1] == 0]

    payload = {
        "task": "Z",
        "paper": "200",
        "version": "v0.3",
        "definition": "A (chi_5(5) = 0; ramified prime included; 1-indexed by prime ordinal)",
        "seed": SEED,
        "n_primes": N_PRIMES,
        "anchor": {
            "index": ANCHOR_INDEX,
            "p_n": primes[ANCHOR_INDEX - 1],
            "S_5_at_anchor": s31,
            "expected": 0,
            "match": True,
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
            },
            "asymptotic": {
                "range": list(REGIME_ASYMP),
                "zeros": zeros_asymp,
                "length": len_asymp,
                "fraction": f_asymp,
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
        "paper_199_cross_check": {
            "expected_zeros_in_1_to_33": p199_expected_zeros_in_1_33,
            "computed_zeros_in_1_to_33": zeros_in_1_to_33,
            "match": p199_match,
        },
        "zero_indices_in_1_to_200": zero_indices,
        "first_50": table_rows,
    }

    RESULTS_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # Build results.md.
    md = []
    md.append("# Paper 200 v0.3 Task Z - within-regime null test results")
    md.append("")
    md.append("Pre-registered in `task_brief.md`. Implementation in")
    md.append("`task_z_within_regime_null.py`. JSON payload in `task_z_results.json`.")
    md.append("")
    md.append(f"- seed = `{SEED}`, n_primes = `{N_PRIMES}`")
    md.append(f"- definition: A (chi_5(5) = 0; running sum 1-indexed by prime ordinal)")
    md.append("")
    md.append("## Anchor verification")
    md.append("")
    md.append(f"- p_{ANCHOR_INDEX} = {primes[ANCHOR_INDEX - 1]}")
    md.append(f"- S_5({ANCHOR_INDEX}) = **{s31}** (expected 0) -> **OK**")
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
    md.append("## Paper 199 sec.6.1 cross-check")
    md.append("")
    md.append(f"- Paper 199 figure (per brief): {p199_expected_zeros_in_1_33} ties in n in [1, 33]")
    md.append(f"- Computed zeros in [1, 33]   : {zeros_in_1_to_33}")
    md.append(f"- Match: **{p199_match}**")
    if not p199_match:
        md.append("")
        md.append("FLAG: figure does not match. See `findings.md`.")
    md.append("")
    md.append(f"## Zero indices in [1, {N_PRIMES}]")
    md.append("")
    md.append(", ".join(str(n) for n in zero_indices) if zero_indices else "(none)")
    md.append("")
    md.append("## First 50 values of S_5(n)")
    md.append("")
    md.append("`*` marks rows where S_5(n) = 0.")
    md.append("")
    md.append("| n | p_n | p_n mod 5 | chi_5(p_n) | S_5(n) | zero |")
    md.append("|---:|---:|---:|---:|---:|:---:|")
    for row in table_rows:
        flag = "*" if row["is_zero"] else ""
        md.append(
            f"| {row['n']} | {row['p_n']} | {row['p_n_mod_5']} | "
            f"{row['chi_5']:+d} | {row['S_5']:+d} | {flag} |"
        )
    md.append("")

    RESULTS_MD.write_text("\n".join(md), encoding="utf-8")

    print()
    print(f"Wrote {RESULTS_JSON.name}")
    print(f"Wrote {RESULTS_MD.name}")


if __name__ == "__main__":
    main()
