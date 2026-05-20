"""
Task 3 - Side-by-side relationship analysis for Paper 190 Phase 1 v0.2.

Pre-reg: Mr_Code_Brief_Paper_190_Phase_1_v0_2.md (commit 971938b).

Sub-tasks (per brief sec 4 Task 3):
 3a Pearson correlation of partial-quotient sequences (Bonferroni-corrected).
    Spearman rank correlation included as robust sanity check.
 3b Block-pattern detection: scan each CF for repeated blocks of length
    2..10; compare to Khinchin-Levy expected frequency.
 3c Cross-CF block matching: blocks of length >= 4 in two different CFs.
 3d Convergent comparison: shared denominators / Fibonacci / Lucas /
    icosahedral integers (12, 20, 30, 60).

Analytical Khinchin-Levy null is used here for closed-form expected
counts. Empirical null distribution (Task 5) replaces the analytical
form where it is inadequate. Findings flagged at p < 0.05 (uncorrected)
and p < Bonferroni for significance.
"""

import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
from mpmath import mp, mpf, log as mlog, nstr

ROOT = Path(__file__).parent
CF_DATA = ROOT / "task_2_cf_data.json"
OUT_MD = ROOT / "task_3_relationships.md"
OUT_JSON = ROOT / "task_3_relationships.json"

# Khinchin-Levy distribution: P(a_k = m) = log2(1 + 1/(m*(m+2))).
# K-L applies to a_k for k >= 1 (a_0 is the integer part, unconstrained).


def kl_prob(m):
    """Khinchin-Levy probability of a single partial quotient being m (m >= 1)."""
    return math.log2(1.0 + 1.0 / (m * (m + 2)))


def kl_table(m_max=10000):
    """Cumulative table of K-L probabilities for m = 1..m_max."""
    return np.array([kl_prob(m) for m in range(1, m_max + 1)])


KL_TABLE = kl_table(10000)
KL_TAIL = 1.0 - KL_TABLE.sum()   # probability mass at m > m_max (tiny)
KL_q = float((KL_TABLE ** 2).sum())   # collision probability under K-L


# ---------------------------------------------------------------------------
# Fibonacci / Lucas / icosahedral integer sets (for 3d)
# ---------------------------------------------------------------------------

def fibonacci_up_to(N):
    fibs = [1, 1]
    while fibs[-1] + fibs[-2] <= N:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


def lucas_up_to(N):
    lucs = [1, 3]
    while lucs[-1] + lucs[-2] <= N:
        lucs.append(lucs[-1] + lucs[-2])
    return lucs


def icosahedral_integers():
    return {12, 20, 30, 60}


# ---------------------------------------------------------------------------
# Load Task 2 CF data
# ---------------------------------------------------------------------------

with open(CF_DATA, "r", encoding="utf-8") as f:
    DATA = json.load(f)

ROWS = ["a", "b", "c", "d", "e", "f", "g", "h"]
NAMES = {k: DATA["candidates"][k]["name"] for k in ROWS}
CFS = {k: DATA["candidates"][k]["cf"] for k in ROWS}
CONVS = {k: DATA["candidates"][k]["convergents"] for k in ROWS}

# Partial quotients a_k for k >= 1 (drop a_0 for K-L analysis)
PQS = {k: CFS[k][1:] for k in ROWS}
LENGTHS = {k: len(PQS[k]) for k in ROWS}


# ---------------------------------------------------------------------------
# 3a Pearson + Spearman correlation matrix
# ---------------------------------------------------------------------------

def task_3a():
    print("\n=== Task 3a: pairwise correlations ===")
    pairs = []
    for i, ri in enumerate(ROWS):
        for j, rj in enumerate(ROWS):
            if i >= j:
                continue
            xi = PQS[ri]
            xj = PQS[rj]
            n = min(len(xi), len(xj))
            if n < 3:
                pairs.append((ri, rj, n, None, None, None, None))
                continue
            xi_t = np.array(xi[:n], dtype=float)
            xj_t = np.array(xj[:n], dtype=float)
            r_p, p_p = stats.pearsonr(xi_t, xj_t)
            r_s, p_s = stats.spearmanr(xi_t, xj_t)
            pairs.append((ri, rj, n, r_p, p_p, r_s, p_s))

    n_pairs = len(pairs)
    bonf = 0.05 / n_pairs
    print(f"Number of pairs: {n_pairs}; Bonferroni threshold = 0.05/{n_pairs} = {bonf:.5g}")

    rows_out = []
    flags_uncorr = []
    flags_bonf = []
    for (ri, rj, n, r_p, p_p, r_s, p_s) in pairs:
        flag_u = (p_p is not None and p_p < 0.05) or (p_s is not None and p_s < 0.05)
        flag_b = (p_p is not None and p_p < bonf) or (p_s is not None and p_s < bonf)
        if flag_u:
            flags_uncorr.append((ri, rj))
        if flag_b:
            flags_bonf.append((ri, rj))
        rows_out.append({
            "row_i": ri, "row_j": rj, "n_paired": n,
            "pearson_r": r_p, "pearson_p": p_p,
            "spearman_r": r_s, "spearman_p": p_s,
            "flag_uncorrected": flag_u,
            "flag_bonferroni": flag_b,
        })
        if r_p is not None:
            print(f"  ({ri},{rj})  n={n:4d}  Pearson r={r_p:+.4f} p={p_p:.4g}  "
                  f"Spearman r={r_s:+.4f} p={p_s:.4g}  "
                  f"{'BONF' if flag_b else ('uncorr' if flag_u else '')}")
        else:
            print(f"  ({ri},{rj})  n={n:4d}  (insufficient data)")
    return {
        "pairs": rows_out,
        "n_pairs": n_pairs,
        "bonferroni_threshold": bonf,
        "flags_uncorrected": flags_uncorr,
        "flags_bonferroni": flags_bonf,
    }


# ---------------------------------------------------------------------------
# 3b Block-pattern detection (per-row)
# ---------------------------------------------------------------------------

def block_kl_expected_count(block, sequence_length):
    """Expected number of occurrences of a specific block in a K-L iid sequence."""
    prob = 1.0
    for m in block:
        if m < 1 or m > len(KL_TABLE):
            return None
        prob *= float(KL_TABLE[m - 1])
    n_slots = sequence_length - len(block) + 1
    if n_slots <= 0:
        return 0.0
    return n_slots * prob


def task_3b():
    print("\n=== Task 3b: within-row block-pattern detection ===")
    results = {}
    for k in ROWS:
        pqs = PQS[k]
        if len(pqs) < 4:
            results[k] = {"length": len(pqs), "blocks": {}, "outliers": []}
            print(f"  ({k}) length={len(pqs)} - too short, skipping")
            continue
        row_blocks = {}
        outliers = []
        for L in range(2, 11):
            counter = Counter(tuple(pqs[i:i + L]) for i in range(len(pqs) - L + 1))
            top = counter.most_common(20)
            block_data = []
            for block, count in top:
                if count < 2:
                    break
                expected = block_kl_expected_count(block, len(pqs))
                if expected is None or expected == 0:
                    z = None
                else:
                    # Poisson-approx z-score (count - lambda) / sqrt(lambda)
                    z = (count - expected) / math.sqrt(expected)
                block_data.append({
                    "block": list(block),
                    "count": count,
                    "expected_KL": expected,
                    "z_poisson": z,
                })
                if z is not None and z > 3.0:
                    outliers.append({"row": k, "block": list(block),
                                      "count": count, "expected": expected, "z": z})
            row_blocks[L] = block_data
        results[k] = {"length": len(pqs), "blocks": row_blocks, "outliers": outliers}
        if outliers:
            print(f"  ({k}) length={len(pqs)} -- outliers (z>3 vs K-L):")
            for o in outliers[:10]:
                bs = ",".join(str(x) for x in o['block'])
                print(f"      block=({bs}) count={o['count']} expected={o['expected']:.3f} "
                      f"z={o['z']:.2f}")
        else:
            print(f"  ({k}) length={len(pqs)} -- no z>3 outliers vs K-L")
    return results


# ---------------------------------------------------------------------------
# 3c Cross-CF block matching (length >= 4)
# ---------------------------------------------------------------------------

def task_3c():
    print("\n=== Task 3c: cross-CF block matching (lengths 4..10) ===")
    cross = {}
    for i, ri in enumerate(ROWS):
        for j, rj in enumerate(ROWS):
            if i >= j:
                continue
            xi = PQS[ri]; xj = PQS[rj]
            if len(xi) < 4 or len(xj) < 4:
                continue
            len_i = len(xi); len_j = len(xj)
            pair_key = f"{ri}-{rj}"
            cross[pair_key] = {}
            for L in range(4, 11):
                if len_i < L or len_j < L:
                    continue
                blocks_i = set()
                blocks_pos_i = defaultdict(list)
                for k in range(len_i - L + 1):
                    b = tuple(xi[k:k + L])
                    blocks_i.add(b)
                    blocks_pos_i[b].append(k)
                matches = []
                for k in range(len_j - L + 1):
                    b = tuple(xj[k:k + L])
                    if b in blocks_i:
                        matches.append({"block": list(b),
                                         "pos_i": blocks_pos_i[b],
                                         "pos_j": k,
                                         "kl_prob": math.prod(
                                             float(KL_TABLE[m - 1])
                                             if 1 <= m <= len(KL_TABLE) else 1e-30
                                             for m in b),
                                         })
                expected_under_kl = (len_i - L + 1) * (len_j - L + 1) * (KL_q ** L)
                cross[pair_key][L] = {
                    "matches": matches,
                    "n_observed": len(matches),
                    "n_expected_KL": expected_under_kl,
                }
            # Summary line
            summary = ", ".join(f"L={L}: obs={cross[pair_key][L]['n_observed']} "
                                 f"exp={cross[pair_key][L]['n_expected_KL']:.2f}"
                                 for L in sorted(cross[pair_key].keys()))
            print(f"  ({ri},{rj}): {summary}")
    # Top surprise matches (lowest individual K-L probability)
    surprise = []
    for pair_key, by_L in cross.items():
        for L, dat in by_L.items():
            for m in dat["matches"]:
                surprise.append((m["kl_prob"], pair_key, L, m["block"],
                                  m["pos_i"], m["pos_j"]))
    surprise.sort()
    print("\n  Top 20 rarest cross-CF matches (lowest K-L probability):")
    for (p, pk, L, blk, pi, pj) in surprise[:20]:
        print(f"    {pk}  L={L}  block={blk}  KL-prob={p:.3g}  pos_i={pi[:3]} pos_j={pj}")
    return {"per_pair": cross, "top_surprise": [
        {"kl_prob": p, "pair": pk, "L": L, "block": blk, "pos_i": pi, "pos_j": pj}
        for (p, pk, L, blk, pi, pj) in surprise[:50]
    ]}


# ---------------------------------------------------------------------------
# 3d Convergent comparison
# ---------------------------------------------------------------------------

def task_3d():
    print("\n=== Task 3d: convergent comparison ===")
    # Max denominator across all rows to size the Fib/Lucas tables
    max_q = max(c[1] for row in ROWS for c in CONVS[row])
    fibs = set(fibonacci_up_to(max_q))
    lucs = set(lucas_up_to(max_q))
    icos = icosahedral_integers()

    per_row = {}
    for k in ROWS:
        rows_data = []
        for idx, (p, q) in enumerate(CONVS[k]):
            tags = []
            if q in fibs:
                tags.append("Fib")
            if q in lucs:
                tags.append("Lucas")
            if q in icos:
                tags.append("Icos")
            rows_data.append({"k": idx, "p": p, "q": q, "tags": tags})
        per_row[k] = rows_data

    # Shared denominators across rows
    q_to_rows = defaultdict(list)
    for k in ROWS:
        for idx, (p, q) in enumerate(CONVS[k]):
            q_to_rows[q].append((k, idx))
    shared = {q: pos for q, pos in q_to_rows.items() if len(pos) > 1 and q > 1}

    # Common factors among denominators of last reported convergent per row
    last_q = {k: CONVS[k][-1][1] for k in ROWS}

    # Print per-row tags
    for k in ROWS:
        tags_line = []
        for d in per_row[k]:
            label = f"q_{d['k']}={d['q']}"
            if d['tags']:
                label += f"[{','.join(d['tags'])}]"
            tags_line.append(label)
        print(f"  ({k}) {NAMES[k]}: {', '.join(tags_line)}")
    print(f"\n  Shared denominators (q appears in 2+ rows, q > 1):")
    for q in sorted(shared.keys()):
        positions = ", ".join(f"({r},k={i})" for r, i in shared[q])
        print(f"    q={q}: {positions}")

    return {"per_row": per_row, "shared_denominators": {
        str(q): [{"row": r, "k": i} for (r, i) in pos]
        for q, pos in shared.items()
    }}


# ---------------------------------------------------------------------------
# Aggregate output
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Task 3 - Side-by-side relationship analysis")
    print("=" * 70)
    print(f"K-L collision probability q = sum P(m)^2 = {KL_q:.6f}")
    print(f"K-L tail mass at m > {len(KL_TABLE)}: {KL_TAIL:.3e}")
    for k in ROWS:
        print(f"  Row ({k}) {NAMES[k]}: CF length = {len(CFS[k])}, "
              f"K-L analysable PQ count = {LENGTHS[k]}")

    res_a = task_3a()
    res_b = task_3b()
    res_c = task_3c()
    res_d = task_3d()

    # Save JSON
    output = {
        "metadata": {
            "pre_reg_commit": "971938b",
            "task": "Task 3 of Mr_Code_Brief_Paper_190_Phase_1_v0_2.md",
            "kl_collision_probability_q": KL_q,
        },
        "task_3a_correlations": res_a,
        "task_3b_block_patterns": res_b,
        "task_3c_cross_block_matching": res_c,
        "task_3d_convergents": res_d,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT_JSON}")

    # Write markdown summary
    md = []
    md.append("# Task 3 - Side-by-side relationship analysis")
    md.append("")
    md.append(f"**Pre-reg commit:** `971938b`")
    md.append("")
    md.append("## 3a Pairwise correlations")
    md.append("")
    md.append(f"Pairs tested: **{res_a['n_pairs']}** "
              f"(Bonferroni threshold = 0.05/{res_a['n_pairs']} = "
              f"{res_a['bonferroni_threshold']:.4g}).")
    md.append("")
    md.append("| pair | n | Pearson r | Pearson p | Spearman r | Spearman p | flag |")
    md.append("|------|---|-----------|-----------|------------|------------|------|")
    for p in res_a["pairs"]:
        if p["pearson_r"] is None:
            md.append(f"| ({p['row_i']},{p['row_j']}) | {p['n_paired']} | - | - | - | - | n/a |")
            continue
        flag = ("**BONF**" if p["flag_bonferroni"]
                else ("uncorr" if p["flag_uncorrected"] else ""))
        md.append(f"| ({p['row_i']},{p['row_j']}) | {p['n_paired']} | "
                  f"{p['pearson_r']:+.4f} | {p['pearson_p']:.4g} | "
                  f"{p['spearman_r']:+.4f} | {p['spearman_p']:.4g} | {flag} |")
    md.append("")

    md.append("## 3b Within-row block-pattern outliers (z > 3 vs K-L)")
    md.append("")
    md.append("| row | block length | block | count | expected (K-L) | z |")
    md.append("|-----|--------------|-------|-------|----------------|---|")
    any_outliers = False
    for k in ROWS:
        for o in res_b[k]["outliers"]:
            any_outliers = True
            bs = ",".join(str(x) for x in o["block"])
            md.append(f"| {k} | {len(o['block'])} | ({bs}) | {o['count']} | "
                      f"{o['expected']:.3f} | {o['z']:.2f} |")
    if not any_outliers:
        md.append("| - | - | - | - | - | - |")
        md.append("")
        md.append("*No within-row block z > 3 vs K-L found.*")
    md.append("")

    md.append("## 3c Cross-CF block matching summary (length 4 only shown)")
    md.append("")
    md.append("| pair | observed | expected (K-L) | obs/exp |")
    md.append("|------|----------|----------------|---------|")
    for pair_key, by_L in res_c["per_pair"].items():
        if 4 in by_L:
            d = by_L[4]
            ratio = d["n_observed"] / d["n_expected_KL"] if d["n_expected_KL"] > 0 else float('inf')
            md.append(f"| {pair_key} | {d['n_observed']} | {d['n_expected_KL']:.2f} | {ratio:.2f} |")
    md.append("")
    md.append("### Top 20 rarest cross-CF matches (lowest K-L probability)")
    md.append("")
    md.append("| pair | L | block | K-L prob | pos_i (first) | pos_j |")
    md.append("|------|---|-------|----------|---------------|-------|")
    for s in res_c["top_surprise"][:20]:
        blk = ",".join(str(x) for x in s["block"])
        pos_i = s["pos_i"][0] if s["pos_i"] else "-"
        md.append(f"| {s['pair']} | {s['L']} | ({blk}) | {s['kl_prob']:.3g} | "
                  f"{pos_i} | {s['pos_j']} |")
    md.append("")

    md.append("## 3d Convergent comparison")
    md.append("")
    md.append("Tags: **Fib** = in Fibonacci sequence, **Lucas** = in Lucas sequence, "
              "**Icos** = in icosahedral integers {12, 20, 30, 60}.")
    md.append("")
    for k in ROWS:
        md.append(f"### Row ({k}) {NAMES[k]}")
        md.append("")
        md.append("| k | p_k | q_k | tags |")
        md.append("|---|-----|-----|------|")
        for d in res_d["per_row"][k]:
            tags = ",".join(d["tags"]) if d["tags"] else ""
            md.append(f"| {d['k']} | {d['p']} | {d['q']} | {tags} |")
        md.append("")

    md.append("### Shared denominators across rows (q > 1)")
    md.append("")
    md.append("| q | rows |")
    md.append("|---|------|")
    for q, positions in sorted(((int(q), pos) for q, pos in
                                  res_d["shared_denominators"].items())):
        ps = ", ".join(f"({p['row']},k={p['k']})" for p in positions)
        md.append(f"| {q} | {ps} |")
    md.append("")

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
