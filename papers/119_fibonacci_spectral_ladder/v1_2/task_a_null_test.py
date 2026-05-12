#!/usr/bin/env python3
"""Task A: null test for rail coverage. Paper 119 v1.2.

Pre-registered in papers/119_fibonacci_spectral_ladder/v1_2/task_brief.md.
Random seed and N_TRIALS are fixed below; do not edit between pre-reg and run.
"""

from __future__ import annotations

import json
import random
import statistics
from pathlib import Path

SEED = 119
N_TRIALS = 100_000

DATA = {1, 2, 3, 4, 5, 6, 7, 8, 13}
FIBONACCI = {1, 2, 3, 5, 8, 13}
LUCAS = {1, 2, 3, 4, 7, 11}
FACE = {6}
RAILS = FIBONACCI | LUCAS | FACE

RANGE_FULL = set(range(1, 14))
RANGE_TRIVIAL = set(range(1, 9))
RANGE_DISCRIMINATING = set(range(9, 14))

SQUARES = {n * n for n in range(1, 4)}
TRIANGULARS = {n * (n + 1) // 2 for n in range(1, 5)}
PRIMES = {2, 3, 5, 7, 11, 13}
COMBINED_STRUCTURED = SQUARES | TRIANGULARS | PRIMES

SEQUENCES = {
    "Rails": RAILS,
    "Squares": SQUARES,
    "Triangulars": TRIANGULARS,
    "Primes": PRIMES,
    "Combined": COMBINED_STRUCTURED,
}


def coverage(data: set[int], seq: set[int], rng: set[int]) -> float | None:
    denom = data & rng
    if not denom:
        return None
    return len(denom & seq) / len(denom)


def precision(data: set[int], seq: set[int], rng: set[int]) -> float | None:
    denom = seq & rng
    if not denom:
        return None
    return len(data & denom) / len(denom)


def fmt_set(s: set[int]) -> str:
    if not s:
        return "{}"
    return "{" + ", ".join(str(x) for x in sorted(s)) + "}"


def fmt_frac(num: int, den: int) -> str:
    if den == 0:
        return "n/a"
    return f"{num}/{den} = {num/den:.2f}"


def main() -> None:
    out_dir = Path(__file__).parent

    # Sub-table 1: coverage by sequence
    table1_rows = []
    for name, seq in SEQUENCES.items():
        d_full = DATA & RANGE_FULL
        d_triv = DATA & RANGE_TRIVIAL
        d_disc = DATA & RANGE_DISCRIMINATING
        h_full = DATA & seq & RANGE_FULL
        h_triv = DATA & seq & RANGE_TRIVIAL
        h_disc = DATA & seq & RANGE_DISCRIMINATING
        table1_rows.append(
            {
                "name": name,
                "size": len(seq),
                "members": sorted(seq),
                "cov_full_num": len(h_full),
                "cov_full_den": len(d_full),
                "cov_trivial_num": len(h_triv),
                "cov_trivial_den": len(d_triv),
                "cov_discrim_num": len(h_disc),
                "cov_discrim_den": len(d_disc),
            }
        )

    # Sub-table 2: precision in [9, 13]
    table2_rows = []
    for name, seq in SEQUENCES.items():
        preds = sorted(seq & RANGE_DISCRIMINATING)
        hits = sorted(DATA & seq & RANGE_DISCRIMINATING)
        table2_rows.append(
            {
                "name": name,
                "predictions": preds,
                "hits": hits,
                "precision_num": len(hits),
                "precision_den": len(preds),
            }
        )

    # Sub-table 3: random null
    random.seed(SEED)
    pool = sorted(RANGE_FULL)
    sample_size = len(RAILS)

    cov_full_dist = []
    cov_disc_dist = []
    for _ in range(N_TRIALS):
        sample = set(random.sample(pool, sample_size))
        cov_full_dist.append(coverage(DATA, sample, RANGE_FULL))
        cov_disc_dist.append(coverage(DATA, sample, RANGE_DISCRIMINATING))

    rail_cov_full = coverage(DATA, RAILS, RANGE_FULL)
    rail_cov_disc = coverage(DATA, RAILS, RANGE_DISCRIMINATING)

    mean_full = statistics.mean(cov_full_dist)
    std_full = statistics.stdev(cov_full_dist)
    mean_disc = statistics.mean(cov_disc_dist)
    std_disc = statistics.stdev(cov_disc_dist)

    p_full_ge = sum(1 for c in cov_full_dist if c >= rail_cov_full) / N_TRIALS
    p_disc_ge = sum(1 for c in cov_disc_dist if c >= rail_cov_disc) / N_TRIALS
    pct_full = sum(1 for c in cov_full_dist if c <= rail_cov_full) / N_TRIALS * 100
    pct_disc = sum(1 for c in cov_disc_dist if c <= rail_cov_disc) / N_TRIALS * 100

    results = {
        "seed": SEED,
        "n_trials": N_TRIALS,
        "data": sorted(DATA),
        "rails": sorted(RAILS),
        "rail_components": {
            "fibonacci": sorted(FIBONACCI),
            "lucas": sorted(LUCAS),
            "face": sorted(FACE),
        },
        "ranges": {
            "full": sorted(RANGE_FULL),
            "trivial": sorted(RANGE_TRIVIAL),
            "discriminating": sorted(RANGE_DISCRIMINATING),
        },
        "sub_table_1_coverage": table1_rows,
        "sub_table_2_precision_discrim": table2_rows,
        "sub_table_3_random_null": {
            "rail_coverage_full": rail_cov_full,
            "rail_coverage_discrim": rail_cov_disc,
            "mean_full": mean_full,
            "std_full": std_full,
            "mean_discrim": mean_disc,
            "std_discrim": std_disc,
            "rail_percentile_full": pct_full,
            "rail_percentile_discrim": pct_disc,
            "p_value_full_ge_rail": p_full_ge,
            "p_value_discrim_ge_rail": p_disc_ge,
        },
    }

    (out_dir / "task_a_results.json").write_text(
        json.dumps(results, indent=2, sort_keys=False), encoding="utf-8"
    )

    # Markdown
    md = []
    md.append("# Task A — Null test for rail coverage")
    md.append("")
    md.append("Paper 119 v1.2 computational task A. Pre-registration: `task_brief.md`.")
    md.append("")
    md.append(f"- Seed: `{SEED}`")
    md.append(f"- Trials: `{N_TRIALS:,}`")
    md.append(f"- Realised data: `{fmt_set(DATA)}`")
    md.append(f"- Rails (Fibonacci ∪ Lucas ∪ {{6}}): `{fmt_set(RAILS)}` (size {len(RAILS)})")
    md.append(f"- Discriminating range [9, 13]: `{fmt_set(RANGE_DISCRIMINATING)}`")
    md.append("")
    md.append("## Sub-table 1 — Coverage by sequence")
    md.append("")
    md.append("| Sequence | Size | Cov [1,13] | Cov [1,8] | Cov [9,13] |")
    md.append("|---|---|---|---|---|")
    for r in table1_rows:
        md.append(
            "| {n} | {s} | {a} | {b} | {c} |".format(
                n=r["name"],
                s=r["size"],
                a=fmt_frac(r["cov_full_num"], r["cov_full_den"]),
                b=fmt_frac(r["cov_trivial_num"], r["cov_trivial_den"]),
                c=fmt_frac(r["cov_discrim_num"], r["cov_discrim_den"]),
            )
        )
    md.append("")
    md.append("## Sub-table 2 — Precision in [9, 13]")
    md.append("")
    md.append("| Sequence | Predictions in [9,13] | Hits | Precision |")
    md.append("|---|---|---|---|")
    for r in table2_rows:
        preds_str = fmt_set(set(r["predictions"]))
        hits_str = fmt_set(set(r["hits"]))
        prec_str = fmt_frac(r["precision_num"], r["precision_den"]) if r["precision_den"] else "n/a (no predictions)"
        md.append(f"| {r['name']} | {preds_str} | {hits_str} | {prec_str} |")
    md.append("")
    md.append(f"## Sub-table 3 — Random null distribution ({N_TRIALS:,} trials, size-10 subsets of [1,13])")
    md.append("")
    md.append("| Statistic | Coverage [1,13] | Coverage [9,13] |")
    md.append("|---|---|---|")
    md.append(f"| Rail value | {rail_cov_full:.4f} | {rail_cov_disc:.4f} |")
    md.append(f"| Mean | {mean_full:.4f} | {mean_disc:.4f} |")
    md.append(f"| Std dev | {std_full:.4f} | {std_disc:.4f} |")
    md.append(f"| Rail percentile (≤ rail) | {pct_full:.2f}% | {pct_disc:.2f}% |")
    md.append(f"| Empirical P(random ≥ rail) | {p_full_ge:.4f} | {p_disc_ge:.4f} |")
    md.append("")
    md.append("## Status flags")
    md.append("")
    md.append("- Realised data and rail structure: DERIVED (Paper 119 §5 branching).")
    md.append("- Choice of comparison sequences: pre-registered (squares, triangulars, primes, random size-10).")
    md.append("- Coverage rates and p-values: OBSERVED.")
    md.append("- Interpretation: OBSERVED-with-empirical-support; see findings.md for the decision against pre-registered criterion.")
    md.append("")
    (out_dir / "task_a_results.md").write_text("\n".join(md), encoding="utf-8")

    print(f"Task A complete. Rail cov [1,13]={rail_cov_full:.4f} ({pct_full:.2f}%-ile), "
          f"Rail cov [9,13]={rail_cov_disc:.4f} ({pct_disc:.2f}%-ile)")
    print(f"P(random >= rail) full = {p_full_ge:.4f}, discrim = {p_disc_ge:.4f}")


if __name__ == "__main__":
    main()
