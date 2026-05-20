"""
Task 4 re-analysis: proper percentile computation + Stouffer combined p.

Loads cached gamma_K(n) values from task_4_em_prefix_prediction.json (saves
the slow H_n re-computation), generates fresh null arrays per (n, K) cell
WITHOUT the linear-interpolation shortcut, computes the actual percentile
of gamma's divergence position in the null array, and combines per-cell
one-sided p-values via Stouffer's z.

Also fixes the analytical prediction: K-L precision consumption rate gives
position ~ -ln(R)/(2L) where L = pi^2/(12 ln 2). The original task 4 code
printed -ln(R)/L (off by factor 2) as the prediction.
"""

import json
import math
from pathlib import Path

import numpy as np
from mpmath import mp, mpf
from scipy.stats import norm

from cf_tools import continued_fraction

ROOT = Path(__file__).parent
TASK4_JSON = ROOT / "task_4_em_prefix_prediction.json"
OUT_JSON = ROOT / "task_4_reanalysis.json"
OUT_MD = ROOT / "task_4_reanalysis.md"

DPS = 500
CF_MAX = 250
N_NULL = 500   # bigger null cohort for proper percentile resolution
SEED = 20260521

LEVY = math.pi ** 2 / (12.0 * math.log(2))  # natural-log Lévy constant


def sample_kl_cf(rng, length):
    U = rng.uniform(0.0, 1.0, size=length)
    U = np.clip(U, 1e-15, 1.0 - 1e-15)
    x = np.power(2.0, 1.0 - U)
    return np.ceil((2.0 - x) / (x - 1.0)).astype(np.int64).tolist()


def cf_to_mpf(cf):
    if not cf:
        return mpf(0)
    val = mpf(cf[-1])
    for k in range(len(cf) - 2, -1, -1):
        val = mpf(cf[k]) + mpf(1) / val
    return val


def divergence_position(cf_a, cf_b):
    for i, (x, y) in enumerate(zip(cf_a, cf_b)):
        if x != y:
            return i
    return min(len(cf_a), len(cf_b))


def main():
    with open(TASK4_JSON, "r", encoding="utf-8") as f:
        cached = json.load(f)
    grid = cached["grid_gamma_K_n"]

    # Generate null cohort once
    rng = np.random.default_rng(SEED)
    mp.dps = DPS
    print(f"Generating {N_NULL} K-L random g* in (0,1)...")
    g_stars = []
    for _ in range(N_NULL):
        cf = [0] + sample_kl_cf(rng, CF_MAX + 50)
        g = cf_to_mpf(cf)
        g_stars.append(g)

    print("\nRe-analysing each (n, K) cell:")
    print(f"  {'n':>8} {'K':>3} {'R':>11} {'KL_pred':>8} {'obs':>5} "
          f"{'null_med':>9} {'null_p90':>9} {'true_pct':>9} {'p_1sided':>9}")

    results = []
    for g_cell in grid:
        n = g_cell["n"]
        K = g_cell["K"]
        R = float(g_cell["R_residual"])
        obs = g_cell["observed_div_pos"]

        # Corrected K-L prediction: pos ~ -ln(R)/(2L)
        kl_pred = -math.log(R) / (2 * LEVY)

        # Build raw null array for this cell
        R_mpf = mpf(R)
        null_divs = []
        # Reset RNG per cell? Use a derived seed to keep this deterministic.
        # We re-use the same g_stars across cells to reduce variance.
        cell_rng = np.random.default_rng(SEED + n * 100 + K)
        for g in g_stars:
            noise = R_mpf * mpf(float(cell_rng.uniform(-1.0, 1.0)))
            g_pert = g + noise
            cf_p = continued_fraction(g_pert, N=CF_MAX)
            cf_c = continued_fraction(g, N=CF_MAX)
            null_divs.append(divergence_position(cf_p, cf_c))
        null_divs = np.array(null_divs)

        # True percentile: P(null >= obs) one-sided upper
        # (γ "more agreement" -> right tail of null)
        n_geq = int((null_divs >= obs).sum())
        n_strict_gt = int((null_divs > obs).sum())
        # One-sided p: probability null >= obs (allowing ties on the conservative side)
        p_one_sided = (n_geq + 0.5) / (len(null_divs) + 1)  # smoothed
        true_percentile = 100.0 * (1.0 - p_one_sided)

        results.append({
            "n": n, "K": K, "R": R,
            "kl_predicted_pos_corrected": kl_pred,
            "observed_div_pos": obs,
            "null_median": float(np.median(null_divs)),
            "null_p90": float(np.percentile(null_divs, 90)),
            "null_p95": float(np.percentile(null_divs, 95)),
            "null_mean": float(np.mean(null_divs)),
            "true_percentile": true_percentile,
            "p_one_sided": p_one_sided,
            "n_null": len(null_divs),
            "n_null_ge_obs": n_geq,
        })

        print(f"  {n:>8} {K:>3} {R:>11.3e} {kl_pred:>8.1f} {obs:>5} "
              f"{np.median(null_divs):>9.1f} {np.percentile(null_divs, 90):>9.1f} "
              f"{true_percentile:>8.1f}% {p_one_sided:>9.4f}")

    # Stouffer's combined z across all 25 cells
    p_vals = np.array([r["p_one_sided"] for r in results])
    z_per_cell = norm.ppf(1.0 - p_vals)   # one-sided z from p-value
    stouffer_z = float(np.sum(z_per_cell) / math.sqrt(len(z_per_cell)))
    stouffer_p = float(1.0 - norm.cdf(stouffer_z))

    # Also compute: count of cells at p < 0.05 individually
    n_signif_individual = int((p_vals < 0.05).sum())
    # Bonferroni for 25 cells
    bonf_thresh = 0.05 / len(p_vals)
    n_signif_bonf = int((p_vals < bonf_thresh).sum())

    print(f"\n=== Combined analysis ===")
    print(f"Cells: {len(results)}")
    print(f"Stouffer's combined z (one-sided): {stouffer_z:.4f}")
    print(f"Stouffer's combined p (one-sided): {stouffer_p:.4f}")
    print(f"Cells at p<0.05 individually: {n_signif_individual}/25")
    print(f"Cells at p<{bonf_thresh:.5g} (Bonferroni): {n_signif_bonf}/25")
    print(f"Mean true percentile: {np.mean([r['true_percentile'] for r in results]):.2f}% "
          f"(H0 expected: 50%)")

    # Save JSON
    output = {
        "metadata": {
            "pre_reg_commit": "971938b",
            "task": "Task 4 re-analysis (proper percentile + Stouffer)",
            "seed": SEED,
            "n_null_per_cell": N_NULL,
            "dps": DPS,
            "levy_constant_natural_log": LEVY,
            "kl_prediction_formula": "pos ~ -ln(R)/(2L), L = pi^2/(12 ln 2)",
        },
        "per_cell": results,
        "combined": {
            "stouffer_z_one_sided": stouffer_z,
            "stouffer_p_one_sided": stouffer_p,
            "n_signif_individual_p05": n_signif_individual,
            "n_signif_bonferroni": n_signif_bonf,
            "bonferroni_threshold": bonf_thresh,
            "mean_true_percentile": float(np.mean([r["true_percentile"] for r in results])),
        },
    }
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT_JSON}")

    # Markdown
    md = []
    md.append("# Task 4 re-analysis")
    md.append("")
    md.append(f"**Pre-reg commit:** `971938b`")
    md.append("")
    md.append("Re-analysis of Task 4 using proper percentile computation against "
              f"raw null arrays (N_NULL = {N_NULL} K-L random g* per cell), and "
              "Stouffer's combined z-test for an overall effect across the 25 "
              "(n, K) cells.")
    md.append("")
    md.append("Corrected K-L prediction: $\\mathrm{pos} \\sim -\\ln(R)/(2L)$ "
              f"with Lévy's constant $L = \\pi^2/(12 \\ln 2) \\approx {LEVY:.4f}$ "
              "(natural-log convention).")
    md.append("")

    md.append("## Per-cell results (true percentile vs raw null array)")
    md.append("")
    md.append("| n | K | R | K-L pred | obs div | null median | null p90 | true %ile | one-sided p |")
    md.append("|---|---|---|----------|---------|-------------|----------|-----------|-------------|")
    for r in results:
        md.append(f"| {r['n']} | {r['K']} | {r['R']:.2e} | "
                  f"{r['kl_predicted_pos_corrected']:.1f} | "
                  f"{r['observed_div_pos']} | "
                  f"{r['null_median']:.1f} | "
                  f"{r['null_p90']:.1f} | "
                  f"{r['true_percentile']:.1f}% | "
                  f"{r['p_one_sided']:.4f} |")
    md.append("")

    md.append("## Combined analysis")
    md.append("")
    md.append(f"- **Cells**: 25")
    md.append(f"- **Stouffer's combined z (one-sided)**: {stouffer_z:.4f}")
    md.append(f"- **Stouffer's combined p (one-sided)**: **{stouffer_p:.4f}**")
    md.append(f"- **Cells at p<0.05 individually**: {n_signif_individual}/25 "
              f"(expected under H0: ~1.25)")
    md.append(f"- **Cells at p<{bonf_thresh:.5g} (Bonferroni)**: {n_signif_bonf}/25")
    md.append(f"- **Mean true percentile**: "
              f"{np.mean([r['true_percentile'] for r in results]):.2f}% "
              "(H0 expected: 50.0%)")
    md.append("")

    md.append("## Interpretation")
    md.append("")
    interp = []
    if stouffer_p < 0.05:
        interp.append(f"Combined Stouffer's p = {stouffer_p:.4f} is below 0.05, "
                       f"crossing the brief's 'necessary condition' threshold for "
                       f"H1 survival (p < 0.05 vs null).")
    else:
        interp.append(f"Combined Stouffer's p = {stouffer_p:.4f} is above 0.05; "
                       f"the brief's 'necessary condition' for H1 survival is not met.")

    if n_signif_individual == 0:
        interp.append("No individual cell reaches p < 0.05.")
    else:
        interp.append(f"{n_signif_individual} cell(s) reach individual p < 0.05.")

    if n_signif_bonf == 0:
        interp.append("No cell survives Bonferroni correction.")

    md.append(" ".join(interp))
    md.append("")

    md.append("If Stouffer's combined p < 0.05 but no individual cell survives "
              "Bonferroni, this is a marginal aggregate effect — γ's CF tracks "
              "itself slightly longer than a random Stieltjes-like constant "
              "approximated at comparable precision, but the bias is small and "
              "the brief's 'sufficient condition' (clear structural pattern, "
              "two-or-more constants at p < 0.001, or clean discriminating "
              "prediction) is not met.")
    md.append("")

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
