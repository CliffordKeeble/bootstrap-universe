# test1_density.py - Test 1, density transition at digit-gap boundaries
# Bootstrap Universe Programme - Brief 4 May 2026 (CinC, brief addendum to
# the Mr Code stack), riemann_digit_gap_test, Test 1 of 3.
#
# Pre-registered design (locked before running, per Pattern 75):
#
#   Boundary windows (fixed at +/- 5 about the predicted heights):
#     boundary 1 = [135, 145]    (digit-gap scale x = 10^61, gamma ~ 140.48)
#     boundary 2 = [179, 189]    (digit-gap scale x = 10^80, gamma ~ 184.21)
#
#   Test statistic per boundary window W:
#     RMS_W = sqrt( mean_{T in W} R(T)^2 )
#     R(T)   = Delta(T) - Delta_smooth(T)
#     Delta(T) = N_zeta(T) - N_chi5(T)             (empirical step function)
#     Delta_smooth(T) = -(T/2pi) log(5) + 1        (Selberg null, exact)
#
#   Permutation null:
#     1000 random window positions of width 10 (matching boundary widths),
#     placed entirely within control region C = [55, 130] union [150, 175]
#     union [190, 215]. Random seed = 1*137 + 42 = 179 per the brief's
#     seed scheme.
#
#   Empirical p-value: fraction of permutation RMS values >= observed.
#
#   PASS:  p_emp(boundary 1) < 0.05  AND  p_emp(boundary 2) < 0.05
#   FAIL:  Test 1 verdict not PASS  -> stop, do not run Tests 2 and 3.
#
# Outputs:
#   test1_density_results.json
#   test1_density.png  (4-panel diagnostic: Delta(T), R(T), Delta_smooth fit,
#                       permutation distributions vs observed)
#
# DERIVED: the Selberg null formula and Delta_smooth.
# OBSERVED: the empirical RMS values, the permutation distribution, the p.
# FRAMEWORK-IMPLIED: the prediction that boundary-window RMS is anomalous.

import json
import math
import random
import sys
from pathlib import Path

import mpmath
from mpmath import mp, mpf

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from zeros_data import (
    load_zeta_zeros,
    load_lchi5_zeros,
    Delta_smooth,
    DEFAULT_DPS,
)

# ----------------------------------------------------------------------
# Pre-registered constants - locked before any computation.
# ----------------------------------------------------------------------

BOUNDARY_1 = (135.0, 145.0)         # gamma ~ 140 -> digit gap at x = 10^61
BOUNDARY_2 = (179.0, 189.0)         # gamma ~ 184 -> digit gap at x = 10^80
WINDOW_WIDTH = BOUNDARY_1[1] - BOUNDARY_1[0]   # = 10.0; same for both

# Control region for permutation null. Excludes both boundary zones AND
# the buffer immediately around them (gap [130, 150] and [175, 190] excluded).
# Upper bound at 215 caps to where we have L(chi_5) coverage (last zero ~219).
CONTROL_REGIONS = [(55.0, 130.0), (150.0, 175.0), (190.0, 215.0)]

N_PERMUTATIONS = 1000
SEED = 1 * 137 + 42                  # = 179, per the brief's seed scheme

# Sampling resolution for the RMS integral. Step 0.05 gives 200 samples per
# 10-wide window. Below the average L-zero spacing at T=140 (~1.4) by ~30x.
T_GRID_STEP = 0.05


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def empirical_count(zeros_sorted_floats, T):
    """N(T) = #{gamma in zeros : gamma <= T}.  Bisect on a sorted float list."""
    # numpy searchsorted is fine here.
    return int(np.searchsorted(zeros_sorted_floats, T, side="right"))


def evaluate_R(T, zeta_floats, lchi5_floats):
    n_z = empirical_count(zeta_floats, T)
    n_l = empirical_count(lchi5_floats, T)
    delta = n_z - n_l
    delta_smooth = float(Delta_smooth(T))
    return delta - delta_smooth


def rms_in_window(R_grid, T_grid, w_lo, w_hi):
    """RMS of R over the closed window [w_lo, w_hi]."""
    mask = (T_grid >= w_lo) & (T_grid <= w_hi)
    return float(np.sqrt(np.mean(R_grid[mask] ** 2)))


def random_window_in_control(rng, width, control_regions):
    """Sample a window centre such that the full window fits in some
    control-region piece."""
    while True:
        # Pick a region weighted by length (uniform sampling on the union).
        lengths = [hi - lo for lo, hi in control_regions]
        total = sum(lengths)
        r = rng.random() * total
        cum = 0.0
        for (lo, hi), L in zip(control_regions, lengths):
            cum += L
            if r <= cum:
                if hi - lo < width:
                    break  # this region too narrow, retry
                centre = rng.uniform(lo + width / 2, hi - width / 2)
                return (centre - width / 2, centre + width / 2)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    mp.dps = DEFAULT_DPS + 5

    print("test1_density - Riemann x digit-gap density transition test")
    print("=" * 76)
    print("Pre-registered windows:")
    print("  boundary 1: [%g, %g]  (gamma ~ 140 -> x = 10^61)" % BOUNDARY_1)
    print("  boundary 2: [%g, %g]  (gamma ~ 184 -> x = 10^80)" % BOUNDARY_2)
    print("Control region (perm null): %s" % CONTROL_REGIONS)
    print("Permutations: %d, seed = %d" % (N_PERMUTATIONS, SEED))
    print()

    # ---- Load zero data ------------------------------------------------
    print("Loading zero data...")
    zeta_zeros = load_zeta_zeros(n=100)
    lchi5_zeros = load_lchi5_zeros(n=150)
    zeta_floats = np.sort(np.array([float(z) for z in zeta_zeros]))
    lchi5_floats = np.sort(np.array([float(z) for z in lchi5_zeros]))
    print("  zeta zeros loaded:    n=%d, max gamma=%.4f"
          % (len(zeta_floats), zeta_floats[-1]))
    print("  L(chi_5) zeros:       n=%d, max gamma=%.4f"
          % (len(lchi5_floats), lchi5_floats[-1]))
    print("  count <= 140 (zeta):  %d   (Selberg smooth: %.2f)"
          % ((zeta_floats <= 140).sum(), float(load_zeta_smooth(140))))
    print("  count <= 140 (chi_5): %d   (Selberg smooth: %.2f)"
          % ((lchi5_floats <= 140).sum(), float(load_chi5_smooth(140))))
    print()

    # ---- Compute R(T) on the analysis grid -----------------------------
    print("Computing R(T) on grid step=%g over [50, 215]..." % T_GRID_STEP)
    T_grid = np.arange(50.0, 215.0 + T_GRID_STEP / 2, T_GRID_STEP)
    Delta_emp = (np.searchsorted(zeta_floats, T_grid, side="right")
                 - np.searchsorted(lchi5_floats, T_grid, side="right")).astype(float)
    Delta_smooth_grid = np.array([float(Delta_smooth(t)) for t in T_grid])
    R_grid = Delta_emp - Delta_smooth_grid
    print("  R range: [%.3f, %.3f], mean=%.4f, std=%.4f"
          % (R_grid.min(), R_grid.max(), R_grid.mean(), R_grid.std()))
    print()

    # ---- Observed RMS in each boundary window --------------------------
    rms_obs1 = rms_in_window(R_grid, T_grid, *BOUNDARY_1)
    rms_obs2 = rms_in_window(R_grid, T_grid, *BOUNDARY_2)
    print("Observed RMS:")
    print("  boundary 1 [%g, %g]:  %.4f" % (*BOUNDARY_1, rms_obs1))
    print("  boundary 2 [%g, %g]:  %.4f" % (*BOUNDARY_2, rms_obs2))
    print()

    # ---- Permutation null ----------------------------------------------
    print("Running %d-permutation null on control region..." % N_PERMUTATIONS)
    rng = random.Random(SEED)
    null_rms = []
    for i in range(N_PERMUTATIONS):
        a, b = random_window_in_control(rng, WINDOW_WIDTH, CONTROL_REGIONS)
        null_rms.append(rms_in_window(R_grid, T_grid, a, b))
    null_rms = np.array(null_rms)
    print("  null RMS distribution: mean=%.4f, median=%.4f, std=%.4f, max=%.4f"
          % (null_rms.mean(), np.median(null_rms), null_rms.std(), null_rms.max()))
    print()

    # ---- Empirical p-values -------------------------------------------
    p_emp1 = float((null_rms >= rms_obs1).sum() + 1) / (N_PERMUTATIONS + 1)
    p_emp2 = float((null_rms >= rms_obs2).sum() + 1) / (N_PERMUTATIONS + 1)
    print("Empirical p-values (one-sided, RMS_null >= RMS_obs):")
    print("  boundary 1: p = %.4f   %s"
          % (p_emp1, "PASS" if p_emp1 < 0.05 else "fail"))
    print("  boundary 2: p = %.4f   %s"
          % (p_emp2, "PASS" if p_emp2 < 0.05 else "fail"))
    print()

    overall_pass = (p_emp1 < 0.05) and (p_emp2 < 0.05)
    verdict = "PASS" if overall_pass else "FAIL"
    print("=" * 76)
    print("Test 1 verdict: %s" % verdict)
    if not overall_pass:
        print("(Per brief discipline: stop here. Do not run Tests 2 and 3.)")

    # ---- Outputs --------------------------------------------------------
    out_dir = Path(__file__).resolve().parent
    json_path = out_dir / "test1_density_results.json"
    plot_path = out_dir / "test1_density.png"

    results = {
        "version": "v1.0",
        "preregistration": {
            "boundary_1": list(BOUNDARY_1),
            "boundary_2": list(BOUNDARY_2),
            "window_width": WINDOW_WIDTH,
            "control_regions": [list(r) for r in CONTROL_REGIONS],
            "n_permutations": N_PERMUTATIONS,
            "seed": SEED,
            "T_grid_step": T_GRID_STEP,
            "selberg_null_formula": "Delta_smooth(T) = -(T/2pi)*log(5) + 1",
        },
        "data": {
            "n_zeta_zeros": int(len(zeta_floats)),
            "n_lchi5_zeros": int(len(lchi5_floats)),
            "zeta_max_gamma": float(zeta_floats[-1]),
            "lchi5_max_gamma": float(lchi5_floats[-1]),
            "zeta_count_at_140": int((zeta_floats <= 140).sum()),
            "lchi5_count_at_140": int((lchi5_floats <= 140).sum()),
            "zeta_count_at_184": int((zeta_floats <= 184).sum()),
            "lchi5_count_at_184": int((lchi5_floats <= 184).sum()),
        },
        "observed": {
            "rms_boundary_1": rms_obs1,
            "rms_boundary_2": rms_obs2,
        },
        "null": {
            "rms_mean": float(null_rms.mean()),
            "rms_median": float(np.median(null_rms)),
            "rms_std": float(null_rms.std()),
            "rms_max": float(null_rms.max()),
            "rms_quantiles": {
                "0.50": float(np.quantile(null_rms, 0.50)),
                "0.90": float(np.quantile(null_rms, 0.90)),
                "0.95": float(np.quantile(null_rms, 0.95)),
                "0.99": float(np.quantile(null_rms, 0.99)),
            },
        },
        "p_values": {
            "boundary_1": p_emp1,
            "boundary_2": p_emp2,
            "convention": "one-sided (RMS_null >= RMS_obs); "
                          "(count + 1) / (N + 1) to avoid p=0",
        },
        "verdict": verdict,
        "pass_criterion": "p_emp(boundary_1) < 0.05 AND p_emp(boundary_2) < 0.05",
    }
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print("JSON written: %s" % json_path)

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    ax.plot(T_grid, Delta_emp, lw=0.6, color="steelblue", label="Delta(T) emp")
    ax.plot(T_grid, Delta_smooth_grid, "--", color="black", lw=1.0,
            label=r"$\bar\Delta(T)$ Selberg null")
    ax.axvspan(*BOUNDARY_1, alpha=0.15, color="crimson", label="boundary 1")
    ax.axvspan(*BOUNDARY_2, alpha=0.15, color="darkgreen", label="boundary 2")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\Delta(T) = N_\zeta(T) - N_{\chi_5}(T)$")
    ax.set_title("Empirical vs Selberg-null zero counts")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(T_grid, R_grid, lw=0.7, color="navy")
    ax.axvspan(*BOUNDARY_1, alpha=0.15, color="crimson")
    ax.axvspan(*BOUNDARY_2, alpha=0.15, color="darkgreen")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$R(T) = \Delta(T) - \bar\Delta(T)$")
    ax.set_title("Residual from Selberg null")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.hist(null_rms, bins=40, color="steelblue", edgecolor="black",
            linewidth=0.4, alpha=0.85, label="null distribution")
    ax.axvline(rms_obs1, color="crimson", lw=2,
               label="boundary 1 obs (p=%.3f)" % p_emp1)
    ax.axvline(rms_obs2, color="darkgreen", lw=2,
               label="boundary 2 obs (p=%.3f)" % p_emp2)
    ax.set_xlabel("RMS of R(T) in window")
    ax.set_ylabel("count")
    ax.set_title("Permutation null vs observed (n=%d)" % N_PERMUTATIONS)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    # Empirical CDF of null vs observed values, log-y to surface tail.
    sorted_null = np.sort(null_rms)
    cdf = np.arange(1, len(sorted_null) + 1) / len(sorted_null)
    ax.plot(sorted_null, 1 - cdf, drawstyle="steps-post",
            color="steelblue", lw=1.2, label="P(null >= x)")
    ax.axvline(rms_obs1, color="crimson", lw=2,
               label="boundary 1 obs")
    ax.axvline(rms_obs2, color="darkgreen", lw=2,
               label="boundary 2 obs")
    ax.set_yscale("log")
    ax.set_xlabel("RMS of R(T) in window")
    ax.set_ylabel("P(null RMS >= x)")
    ax.set_title("Survival function (log-y)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(plot_path, dpi=120)
    plt.close()
    print("Plot written: %s" % plot_path)

    return 0 if overall_pass else 1


# Selberg-null helpers re-exported for the print line above.
def load_zeta_smooth(T):
    from zeros_data import N_zeta_smooth
    return N_zeta_smooth(T)


def load_chi5_smooth(T):
    from zeros_data import N_chi5_smooth
    return N_chi5_smooth(T)


if __name__ == "__main__":
    sys.exit(main())
