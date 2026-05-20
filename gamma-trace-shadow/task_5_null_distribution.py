"""
Task 5 - Null distribution (Khinchin-Levy random Stieltjes-like constants).

Pre-reg: Mr_Code_Brief_Paper_190_Phase_1_v0_2.md (commit 971938b).

Generates K=1000 random length-L=300 CFs with partial quotients sampled
from the Khinchin-Levy distribution (the natural prior for "almost all"
irrationals). For each sample, computes the same statistics as Task 3
and Task 4, building an empirical null distribution per statistic.

Outputs empirical p=0.05 and p=0.001 thresholds plus an evaluation of
each Task 3 flagged finding against the empirical null.

Random sampling
- Inverse-CDF: F_KL(m) = 1 - log2((m+2)/(m+1)),
  so given U ~ Uniform(0,1), the partial quotient is
  ceil((2 - 2^(1-U)) / (2^(1-U) - 1)).
- Heavy tail captured up to ~1e10 in float64 (U within 1e-10 of 1).
- Seed: 20260520 (today's date as int) for reproducibility.
"""

import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).parent
TASK3 = ROOT / "task_3_relationships.json"
TASK2 = ROOT / "task_2_cf_data.json"
OUT_JSON = ROOT / "task_5_null_distribution.json"
OUT_MD = ROOT / "task_5_null_distribution.md"

SEED = 20260520
N_SAMPLES = 1000
LEN_CF = 300
M_PAIRS = 1000   # number of random pairs sampled for pairwise null

# Khinchin-Levy table (analytic, for comparison)
def kl_prob(m):
    return math.log2(1.0 + 1.0 / (m * (m + 2)))


KL_TABLE = np.array([kl_prob(m) for m in range(1, 10001)])
KL_q = float((KL_TABLE ** 2).sum())


def sample_kl(rng, n):
    """Sample n partial quotients iid from K-L distribution."""
    U = rng.uniform(0.0, 1.0, size=n)
    # Avoid U exactly at 0 or 1 (rng.uniform excludes 1 already).
    U = np.clip(U, 1e-15, 1.0 - 1e-15)
    x = np.power(2.0, 1.0 - U)
    m_cont = (2.0 - x) / (x - 1.0)
    return np.ceil(m_cont).astype(np.int64)


# Fibonacci / Lucas
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


def convergent_denominators(cf, m=10):
    """Compute first m convergent denominators from a CF list."""
    q_prev, q_curr = 1, 0
    qs = []
    for a in cf[:m]:
        q_next = a * q_curr + q_prev
        qs.append(q_next)
        q_prev, q_curr = q_curr, q_next
    return qs


def block_kl_expected_count(block, sequence_length):
    prob = 1.0
    for m in block:
        if m < 1 or m > len(KL_TABLE):
            return None
        prob *= float(KL_TABLE[m - 1])
    n_slots = sequence_length - len(block) + 1
    if n_slots <= 0:
        return 0.0
    return n_slots * prob


def max_z_in_cf(cf_pqs, block_length):
    """Maximum z-score (count vs K-L expected) across all blocks of given length."""
    counter = Counter(tuple(cf_pqs[i:i + block_length])
                      for i in range(len(cf_pqs) - block_length + 1))
    max_z = -math.inf
    max_block = None
    max_count = 0
    for block, count in counter.items():
        if count < 2:
            continue
        expected = block_kl_expected_count(block, len(cf_pqs))
        if expected is None or expected == 0:
            continue
        z = (count - expected) / math.sqrt(expected)
        if z > max_z:
            max_z = z
            max_block = block
            max_count = count
    return max_z, max_block, max_count


def count_common_blocks(cf_a, cf_b, block_length):
    blocks_a = set()
    for k in range(len(cf_a) - block_length + 1):
        blocks_a.add(tuple(cf_a[k:k + block_length]))
    count = 0
    for k in range(len(cf_b) - block_length + 1):
        if tuple(cf_b[k:k + block_length]) in blocks_a:
            count += 1
    return count


def shared_prefix_length(cf_a, cf_b):
    n = 0
    for x, y in zip(cf_a, cf_b):
        if x == y:
            n += 1
        else:
            break
    return n


# ---------------------------------------------------------------------------
# Main sampling
# ---------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(SEED)
    print(f"Generating {N_SAMPLES} random K-L CFs of length {LEN_CF} (seed={SEED})...")
    cfs = [sample_kl(rng, LEN_CF).tolist() for _ in range(N_SAMPLES)]

    # ---- Per-CF statistics ----
    print(f"\nComputing per-CF statistics over {N_SAMPLES} samples...")
    max_z_by_L = {L: [] for L in range(2, 11)}
    for ci, cf in enumerate(cfs):
        if ci % 200 == 0:
            print(f"  [{ci}/{N_SAMPLES}] per-CF stats...")
        for L in range(2, 11):
            z, _, _ = max_z_in_cf(cf, L)
            max_z_by_L[L].append(z if z != -math.inf else 0.0)

    # ---- Pairwise statistics ----
    print(f"\nComputing pairwise statistics over {M_PAIRS} random pairs...")
    pearson_r = []
    spearman_r = []
    common_blk_by_L = {L: [] for L in range(4, 11)}
    shared_prefix = []
    for pi in range(M_PAIRS):
        if pi % 200 == 0:
            print(f"  [{pi}/{M_PAIRS}] pairwise stats...")
        i, j = rng.choice(N_SAMPLES, size=2, replace=False)
        a = np.array(cfs[i], dtype=float)
        b = np.array(cfs[j], dtype=float)
        rp, _ = stats.pearsonr(a, b)
        rs, _ = stats.spearmanr(a, b)
        pearson_r.append(rp)
        spearman_r.append(rs)
        for L in range(4, 11):
            common_blk_by_L[L].append(count_common_blocks(cfs[i], cfs[j], L))
        shared_prefix.append(shared_prefix_length(cfs[i], cfs[j]))

    # ---- Convergent / Fib / Lucas / Icos ----
    print(f"\nConvergent + Fib/Lucas/Icos analysis over {N_SAMPLES} samples...")
    fib_set = set(fibonacci_up_to(10**8))
    luc_set = set(lucas_up_to(10**8))
    icos_set = {12, 20, 30, 60}

    fib_hits_first10 = []
    luc_hits_first10 = []
    icos_hits_first10 = []
    L_10_position_hit_F = 0
    L_10_position_hit_L = 0
    for cf in cfs:
        qs = convergent_denominators(cf, m=10)
        nf = sum(1 for q in qs if q in fib_set and q > 1)
        nl = sum(1 for q in qs if q in luc_set and q > 1)
        ni = sum(1 for q in qs if q in icos_set)
        fib_hits_first10.append(nf)
        luc_hits_first10.append(nl)
        icos_hits_first10.append(ni)
        if len(qs) >= 8 and qs[7] in luc_set and qs[7] > 1:
            L_10_position_hit_L += 1
        if len(qs) >= 3 and qs[2] in fib_set and qs[2] > 1:
            L_10_position_hit_F += 1

    # ---- Empirical thresholds ----
    def thresholds(arr):
        arr = np.sort(arr)
        return {
            "min": float(arr[0]),
            "p_05": float(arr[int(0.95 * len(arr))]),       # one-sided upper
            "p_01": float(arr[int(0.99 * len(arr))]),
            "p_001": float(arr[int(0.999 * len(arr))]),
            "max": float(arr[-1]),
            "median": float(arr[len(arr) // 2]),
            "mean": float(np.mean(arr)),
        }

    thresh_maxz = {L: thresholds(max_z_by_L[L]) for L in max_z_by_L}
    thresh_pearson = thresholds(np.abs(np.array(pearson_r)))
    thresh_spearman = thresholds(np.abs(np.array(spearman_r)))
    thresh_common = {L: thresholds(common_blk_by_L[L]) for L in common_blk_by_L}
    thresh_prefix = thresholds(shared_prefix)
    thresh_fib_hits = thresholds(fib_hits_first10)
    thresh_luc_hits = thresholds(luc_hits_first10)
    thresh_icos_hits = thresholds(icos_hits_first10)

    # ---- Evaluate Task 3 observed against empirical null ----
    with open(TASK3, "r", encoding="utf-8") as f:
        task3 = json.load(f)
    with open(TASK2, "r", encoding="utf-8") as f:
        task2 = json.load(f)

    # 3a Pearson against empirical null
    obs_pearson_r = [
        (p["row_i"], p["row_j"], p["pearson_r"])
        for p in task3["task_3a_correlations"]["pairs"]
        if p["pearson_r"] is not None and p["n_paired"] >= 100
    ]
    # 3b max z-score per row (length 4)
    obs_3b_maxz_L4 = []
    for k, row in task3["task_3b_block_patterns"].items():
        if row["length"] < 4:
            continue
        max_z_row = max((o["z"] for o in row["outliers"] if len(o["block"]) == 4),
                         default=0.0)
        obs_3b_maxz_L4.append((k, max_z_row))

    # 3c common length-4 blocks per pair (excluding closeness artifacts)
    obs_3c_L4 = {}
    for pair_key, by_L in task3["task_3c_cross_block_matching"]["per_pair"].items():
        if "4" in by_L or 4 in by_L:
            key = 4 if 4 in by_L else "4"
            obs_3c_L4[pair_key] = by_L[key]["n_observed"]

    # 3d Fib/Lucas hits per row in first 10 convergents
    obs_3d = {}
    for k, conv_data in task3["task_3d_convergents"]["per_row"].items():
        fib_hits = sum(1 for d in conv_data
                       if "Fib" in d["tags"] and d["q"] > 1)
        luc_hits = sum(1 for d in conv_data
                       if "Lucas" in d["tags"] and d["q"] > 1)
        icos_hits = sum(1 for d in conv_data if "Icos" in d["tags"])
        obs_3d[k] = {
            "fib_hits_first10": fib_hits,
            "luc_hits_first10": luc_hits,
            "icos_hits_first10": icos_hits,
        }

    # Empirical p-value helpers
    def emp_p_upper(arr, x):
        return float((np.array(arr) >= x).mean())

    def emp_p_upper_abs(arr, x):
        return float((np.abs(np.array(arr)) >= abs(x)).mean())

    # Evaluate
    eval_3a = []
    pearson_arr = np.array(pearson_r)
    for (ri, rj, r) in obs_pearson_r:
        p_emp = emp_p_upper_abs(pearson_arr, r)
        eval_3a.append({"pair": f"{ri}-{rj}", "r": r, "p_empirical": p_emp})

    eval_3b = []
    for (k, max_z) in obs_3b_maxz_L4:
        if k == "b":
            # Skip the quadratic-irrational periodicity row
            eval_3b.append({"row": k, "max_z_L4": max_z,
                             "p_empirical": "N/A (quadratic irrational, periodic)",
                             "note": "row b is 1/sqrt 3; Lagrange periodicity is structural to quadratic irrationals, not to the question"})
            continue
        p_emp = emp_p_upper(max_z_by_L[4], max_z)
        eval_3b.append({"row": k, "max_z_L4": max_z, "p_empirical": p_emp})

    eval_3c = []
    common_L4_arr = np.array(common_blk_by_L[4])
    for pair_key, obs_count in obs_3c_L4.items():
        p_emp = emp_p_upper(common_L4_arr, obs_count)
        eval_3c.append({"pair": pair_key, "obs_common_L4": obs_count, "p_empirical": p_emp})

    eval_3d = {}
    fib_arr = np.array(fib_hits_first10)
    luc_arr = np.array(luc_hits_first10)
    icos_arr = np.array(icos_hits_first10)
    for k, d in obs_3d.items():
        eval_3d[k] = {
            "fib_hits_first10": d["fib_hits_first10"],
            "fib_p_empirical": emp_p_upper(fib_arr, d["fib_hits_first10"]),
            "luc_hits_first10": d["luc_hits_first10"],
            "luc_p_empirical": emp_p_upper(luc_arr, d["luc_hits_first10"]),
            "icos_hits_first10": d["icos_hits_first10"],
            "icos_p_empirical": emp_p_upper(icos_arr, d["icos_hits_first10"]),
        }

    # ---- Save JSON ----
    output = {
        "metadata": {
            "pre_reg_commit": "971938b",
            "task": "Task 5 of Mr_Code_Brief_Paper_190_Phase_1_v0_2.md",
            "seed": SEED,
            "n_samples": N_SAMPLES,
            "len_cf": LEN_CF,
            "m_pairs": M_PAIRS,
            "kl_collision_probability_q": KL_q,
        },
        "thresholds": {
            "max_z_in_cf_by_block_length": thresh_maxz,
            "pearson_r_abs": thresh_pearson,
            "spearman_r_abs": thresh_spearman,
            "common_blocks_by_length": thresh_common,
            "shared_prefix_length": thresh_prefix,
            "fib_hits_first10": thresh_fib_hits,
            "lucas_hits_first10": thresh_luc_hits,
            "icos_hits_first10": thresh_icos_hits,
        },
        "evaluation_task_3a": eval_3a,
        "evaluation_task_3b": eval_3b,
        "evaluation_task_3c": eval_3c,
        "evaluation_task_3d": eval_3d,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT_JSON}")

    # ---- Markdown summary ----
    md = []
    md.append("# Task 5 - Null distribution (K-L random Stieltjes-like constants)")
    md.append("")
    md.append(f"**Pre-reg commit:** `971938b`  ")
    md.append(f"**Seed:** {SEED}  ")
    md.append(f"**Samples:** {N_SAMPLES} CFs of length {LEN_CF}, "
              f"{M_PAIRS} random pairs.")
    md.append("")
    md.append("## Empirical thresholds")
    md.append("")
    md.append("### Max-z within-CF (vs K-L analytic) by block length")
    md.append("")
    md.append("| L | median | p=0.05 | p=0.01 | p=0.001 | max |")
    md.append("|---|--------|--------|--------|---------|-----|")
    for L in range(2, 11):
        t = thresh_maxz[L]
        md.append(f"| {L} | {t['median']:.2f} | {t['p_05']:.2f} | {t['p_01']:.2f} | "
                  f"{t['p_001']:.2f} | {t['max']:.2f} |")
    md.append("")

    md.append("### |Pearson r| and |Spearman r| (pairwise, n=300)")
    md.append("")
    md.append("| stat | median | p=0.05 | p=0.01 | p=0.001 | max |")
    md.append("|------|--------|--------|--------|---------|-----|")
    for label, t in [("|Pearson r|", thresh_pearson), ("|Spearman r|", thresh_spearman)]:
        md.append(f"| {label} | {t['median']:.4f} | {t['p_05']:.4f} | "
                  f"{t['p_01']:.4f} | {t['p_001']:.4f} | {t['max']:.4f} |")
    md.append("")

    md.append("### Count of common blocks across a random pair, by length")
    md.append("")
    md.append("| L | median | mean | p=0.05 | p=0.01 | p=0.001 | max |")
    md.append("|---|--------|------|--------|--------|---------|-----|")
    for L in range(4, 11):
        t = thresh_common[L]
        md.append(f"| {L} | {t['median']:.0f} | {t['mean']:.1f} | "
                  f"{t['p_05']:.0f} | {t['p_01']:.0f} | "
                  f"{t['p_001']:.0f} | {t['max']:.0f} |")
    md.append("")

    md.append("### Shared CF prefix length (pair, K-L iid)")
    md.append("")
    t = thresh_prefix
    md.append(f"| median | mean | p=0.05 | p=0.01 | p=0.001 | max |")
    md.append(f"|--------|------|--------|--------|---------|-----|")
    md.append(f"| {t['median']:.0f} | {t['mean']:.2f} | "
              f"{t['p_05']:.0f} | {t['p_01']:.0f} | "
              f"{t['p_001']:.0f} | {t['max']:.0f} |")
    md.append("")

    md.append("### Fibonacci / Lucas / Icos hits in first 10 convergents")
    md.append("")
    md.append("| set | mean | p=0.05 | p=0.01 | p=0.001 | max |")
    md.append("|-----|------|--------|--------|---------|-----|")
    for label, t in [("Fibonacci", thresh_fib_hits),
                       ("Lucas", thresh_luc_hits),
                       ("Icosahedral {12,20,30,60}", thresh_icos_hits)]:
        md.append(f"| {label} | {t['mean']:.2f} | {t['p_05']:.0f} | "
                  f"{t['p_01']:.0f} | {t['p_001']:.0f} | {t['max']:.0f} |")
    md.append("")

    md.append("## Re-evaluation of Task 3 flagged findings vs empirical null")
    md.append("")
    md.append("### Task 3a Pearson r (re-evaluated)")
    md.append("")
    md.append("| pair | observed r | empirical p (|r|) | flag |")
    md.append("|------|------------|-------------------|------|")
    for e in eval_3a:
        flag = "**SIGNIF** (p<0.05)" if e["p_empirical"] < 0.05 else "ns"
        md.append(f"| {e['pair']} | {e['r']:+.4f} | {e['p_empirical']:.4f} | {flag} |")
    md.append("")

    md.append("### Task 3b max z (length 4, per row) (re-evaluated)")
    md.append("")
    md.append("| row | observed max-z (L=4) | empirical p | flag |")
    md.append("|-----|----------------------|-------------|------|")
    for e in eval_3b:
        if e["p_empirical"] == "N/A (quadratic irrational, periodic)":
            md.append(f"| {e['row']} | {e['max_z_L4']:.2f} | n/a | excluded "
                      f"(1/sqrt 3 quadratic-irrational periodicity) |")
        else:
            flag = "**SIGNIF**" if e["p_empirical"] < 0.05 else "ns"
            md.append(f"| {e['row']} | {e['max_z_L4']:.2f} | {e['p_empirical']:.4f} | {flag} |")
    md.append("")

    md.append("### Task 3c common length-4 blocks per pair (re-evaluated)")
    md.append("")
    md.append("Note: \"close pairs\" (a,b) and (d,f) are excluded from primary "
              "evaluation because their common blocks are driven by trivial "
              "numerical closeness; closeness artifact is *measured* by "
              "shared-prefix length in the next column.")
    md.append("")
    md.append("| pair | obs common L=4 | empirical p (vs K-L iid pair) | flag |")
    md.append("|------|----------------|-------------------------------|------|")
    for e in eval_3c:
        if e["pair"] in ("a-b", "d-f"):
            md.append(f"| {e['pair']} | {e['obs_common_L4']} | "
                      f"{e['p_empirical']:.4f} | close-pair artifact |")
        else:
            flag = "**SIGNIF**" if e["p_empirical"] < 0.05 else "ns"
            md.append(f"| {e['pair']} | {e['obs_common_L4']} | "
                      f"{e['p_empirical']:.4f} | {flag} |")
    md.append("")

    md.append("### Task 3d Fibonacci / Lucas / Icos hits per row")
    md.append("")
    md.append("| row | Fib obs | Fib p | Lucas obs | Lucas p | Icos obs | Icos p |")
    md.append("|-----|---------|-------|-----------|---------|----------|--------|")
    for k in ["a", "b", "c", "d", "e", "f", "g", "h"]:
        e = eval_3d[k]
        md.append(f"| {k} | {e['fib_hits_first10']} | {e['fib_p_empirical']:.4f} "
                  f"| {e['luc_hits_first10']} | {e['luc_p_empirical']:.4f} "
                  f"| {e['icos_hits_first10']} | {e['icos_p_empirical']:.4f} |")
    md.append("")

    md.append("## Verdict for Phase 1 stop-conditions")
    md.append("")
    # Count flagged-significant items
    n_3a_signif = sum(1 for e in eval_3a if e["p_empirical"] < 0.05)
    n_3b_signif = sum(1 for e in eval_3b
                     if e["p_empirical"] != "N/A (quadratic irrational, periodic)"
                     and e["p_empirical"] < 0.05)
    n_3c_signif = sum(1 for e in eval_3c
                     if e["pair"] not in ("a-b", "d-f") and e["p_empirical"] < 0.05)
    n_3d_signif = sum(1 for k in ["a","b","c","d","e","f","g","h"]
                     if eval_3d[k]["fib_p_empirical"] < 0.05
                     or eval_3d[k]["luc_p_empirical"] < 0.05
                     or eval_3d[k]["icos_p_empirical"] < 0.05)
    md.append(f"- Task 3a empirical-null SIGNIF (p<0.05): **{n_3a_signif} / {len(eval_3a)}**")
    md.append(f"- Task 3b empirical-null SIGNIF (p<0.05, excl. row b periodicity): **{n_3b_signif} / {len(eval_3b)-1}**")
    md.append(f"- Task 3c empirical-null SIGNIF (p<0.05, excl. close pairs): **{n_3c_signif} / {len(eval_3c)-2}**")
    md.append(f"- Task 3d Fib/Lucas/Icos SIGNIF (any of 3 stats per row, p<0.05): **{n_3d_signif} / 8**")
    md.append("")

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
