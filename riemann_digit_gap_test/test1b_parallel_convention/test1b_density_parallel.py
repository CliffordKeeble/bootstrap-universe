# test1b_density_parallel.py - Test 1b, parallel-convention density test
# Bootstrap Universe Programme - Brief 4 May 2026 (CinC), Test 1b of 1 in
# the parallel-convention arc.
#
# This is an INDEPENDENT pre-registered test. NOT a re-run of Test 1, NOT
# a continuation. Hypothesis is the second of two defensible heuristic
# conventions linking digit-gap scale 10^d to zero-height test windows.
# Pre-registration locked in pre_registration.md (committed at the same
# git timestamp as this file, before this script is run).
#
# UTILITY IMPORTS - declared explicitly per CinC's audit-by-inspection
# requirement. The shared functions are the implementations of the lemmas
# (Pattern 34 separation of concerns); locked-parameters block at the top
# of this file is what's pre-registered, NOT the utilities.
#
#   from zeros_data:
#       Delta_smooth(T)               -- Selberg null (DERIVED, audited)
#       N_zeta_smooth, N_chi5_smooth  -- component formulas
#       DEFAULT_DPS                   -- precision constant
#
#   from test1_density:
#       empirical_count(zeros, T)               -- bisect-based N(T)
#       evaluate_R(T, zeta, lchi5)              -- residual at T (unused
#                                                  here; we vectorise)
#       rms_in_window(R_grid, T_grid, lo, hi)   -- windowed RMS
#       random_window_in_control(rng, w, regions) -- permutation sampler
#                                                    (window-fits-inside)
#
# A future reviewer can inspect zeros_data.py and test1_density.py to
# verify these functions match Test 1's. The pre-registered constants
# in this file are independent.
#
# CACHE: reads from lchi5_zeros_dps30_n180.csv (n=166 actual zeros from
# Stage 1 coverage extension, gamma_max = 244.36). Does NOT modify or
# read from the n=150 cache committed for Test 1.

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mpmath
from mpmath import mp, mpf

# Import shared utilities. Path setup so this script can be run directly.
_TEST1B_DIR = Path(__file__).resolve().parent
_RDGT_DIR = _TEST1B_DIR.parent
sys.path.insert(0, str(_RDGT_DIR))
from zeros_data import (
    Delta_smooth,
    N_zeta_smooth,
    N_chi5_smooth,
    DEFAULT_DPS,
)
from test1_density import (
    empirical_count,
    rms_in_window,
    random_window_in_control,
)

# ----------------------------------------------------------------------
# PRE-REGISTERED CONSTANTS - frozen at commit time. Source of truth is
# pre_registration.md (same commit as this file). Any change here
# invalidates the test.
# ----------------------------------------------------------------------

# Resonance-period convention: gamma_d = 2*pi * d / ln(10).
# d=61 -> gamma ~ 166.4524; d=80 -> gamma ~ 218.2982.
BOUNDARY_1 = (161.4524, 171.4524)
BOUNDARY_2 = (213.2982, 223.2982)
WINDOW_WIDTH = BOUNDARY_1[1] - BOUNDARY_1[0]   # = 10.0

# Control regions (post-Stage-1: gamma_max=244.36, C4 valid).
CONTROL_REGIONS = [
    (55.0, 130.0),               # C1
    (150.0, 161.4524),           # C2 - fits one width-10 window
    (189.0, 213.2982),           # C3 - fits two non-overlapping windows
    (223.2982, 239.3594),        # C4 - upper = gamma_max - 5
]

N_PERMUTATIONS = 1000
SEED = 2 * 137 + 42              # = 316; trial=2 in framework seed scheme
T_GRID_STEP = 0.05
T_GRID_MIN = 50.0
T_GRID_MAX = 240.0               # below gamma_max (244.36) by safety margin

# Cache for L(chi_5) zeros - the n=180 file from Stage 1.
LCHI5_CSV = _RDGT_DIR / "zeros_cache" / "lchi5_zeros_dps30_n180.csv"
ZETA_CSV  = _RDGT_DIR / "zeros_cache" / "riemann_zeros_dps30_n100.csv"

OUT_DIR = _TEST1B_DIR


# ----------------------------------------------------------------------
# Helpers (custom CSV loader; everything else imported from above).
# ----------------------------------------------------------------------

def load_zeros_csv(path):
    import csv
    zeros = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            zeros.append(mpf(row["gamma"]))
    return zeros


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    mp.dps = DEFAULT_DPS + 5

    print("test1b_density_parallel - Riemann x digit-gap, resonance-period convention")
    print("=" * 78)
    print("Convention: gamma_d = 2*pi * d / ln(10)")
    print("  d=61 -> gamma ~ 166.4524   window [%g, %g]" % BOUNDARY_1)
    print("  d=80 -> gamma ~ 218.2982   window [%g, %g]" % BOUNDARY_2)
    print()
    print("Control regions: %s" % CONTROL_REGIONS)
    print("Permutations: %d, seed = %d (trial 2 in framework seed scheme)"
          % (N_PERMUTATIONS, SEED))
    print()

    # ---- Load zeros ----------------------------------------------------
    print("Loading zero data...")
    if not ZETA_CSV.exists() or not LCHI5_CSV.exists():
        print("  ERROR: cache files not found. Run Stage 1 (extend_lchi5_sweep.py)"
              " and zeros_data.py self-test first.")
        return 2
    zeta_zeros  = load_zeros_csv(ZETA_CSV)
    lchi5_zeros = load_zeros_csv(LCHI5_CSV)
    zeta_floats  = np.sort(np.array([float(z) for z in zeta_zeros]))
    lchi5_floats = np.sort(np.array([float(z) for z in lchi5_zeros]))
    print("  zeta zeros:    n=%d, max gamma=%.4f"
          % (len(zeta_floats), zeta_floats[-1]))
    print("  L(chi_5) zeros: n=%d, max gamma=%.4f"
          % (len(lchi5_floats), lchi5_floats[-1]))
    print()

    # ---- Compute R(T) on the analysis grid -----------------------------
    print("Computing R(T) on grid step=%g over [%g, %g]..."
          % (T_GRID_STEP, T_GRID_MIN, T_GRID_MAX))
    T_grid = np.arange(T_GRID_MIN, T_GRID_MAX + T_GRID_STEP / 2, T_GRID_STEP)
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
    import random
    print("Running %d-permutation null on control region..." % N_PERMUTATIONS)
    rng = random.Random(SEED)
    null_rms = []
    for i in range(N_PERMUTATIONS):
        a, b = random_window_in_control(rng, WINDOW_WIDTH, CONTROL_REGIONS)
        null_rms.append(rms_in_window(R_grid, T_grid, a, b))
    null_rms = np.array(null_rms)
    print("  null RMS distribution: mean=%.4f, median=%.4f, std=%.4f, max=%.4f"
          % (null_rms.mean(), np.median(null_rms),
             null_rms.std(), null_rms.max()))
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
    print("=" * 78)
    print("Test 1b verdict: %s" % verdict)

    # Map to the four pre-registered outcome interpretations.
    if p_emp1 < 0.05 and p_emp2 < 0.05:
        outcome = "B"
        outcome_label = ("Both boundaries pass - resonance-period convention "
                         "shows predicted signature.")
    elif p_emp1 >= 0.05 and p_emp2 >= 0.05:
        outcome = "A"
        outcome_label = ("Both boundaries fail - density-transition framing "
                         "null at BOTH conventions tested.")
    else:
        outcome = "C"
        outcome_label = ("Asymmetric - one passes, one fails. Reports "
                         "honestly per pre-reg; no third-convention test.")
    if (0.05 <= p_emp1 < 0.15) or (0.05 <= p_emp2 < 0.15):
        outcome += " (with marginal flag D)"
        outcome_label += (
            " WARNING: at least one p in marginal band [0.05, 0.15) - "
            "treat as null per pre-reg."
        )
    print("Pre-registered outcome interpretation: %s" % outcome)
    print("  -> %s" % outcome_label)

    if not overall_pass:
        print()
        print("(Per brief discipline: this is the only test in the parallel-")
        print(" convention arc. No third-convention test launched.)")

    # ---- Outputs --------------------------------------------------------
    json_path = OUT_DIR / "test1b_density_parallel_results.json"
    plot_path = OUT_DIR / "test1b_density_parallel.png"

    results = {
        "version": "v1.0",
        "preregistration_doc": "pre_registration.md",
        "preregistration": {
            "convention": "gamma_d = 2*pi * d / ln(10)",
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
            "convention": ("one-sided (RMS_null >= RMS_obs); "
                           "(count + 1) / (N + 1) Laplace correction"),
        },
        "verdict": verdict,
        "pre_registered_outcome": outcome,
        "pass_criterion": "p_emp(1) < 0.05 AND p_emp(2) < 0.05",
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
    ax.axvspan(*BOUNDARY_1, alpha=0.15, color="crimson", label="boundary 1 (gamma=166.45)")
    ax.axvspan(*BOUNDARY_2, alpha=0.15, color="darkgreen", label="boundary 2 (gamma=218.30)")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\Delta(T) = N_\zeta(T) - N_{\chi_5}(T)$")
    ax.set_title("Empirical vs Selberg-null zero counts (Test 1b)")
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
    sorted_null = np.sort(null_rms)
    cdf = np.arange(1, len(sorted_null) + 1) / len(sorted_null)
    ax.plot(sorted_null, 1 - cdf, drawstyle="steps-post",
            color="steelblue", lw=1.2, label="P(null >= x)")
    ax.axvline(rms_obs1, color="crimson", lw=2, label="boundary 1 obs")
    ax.axvline(rms_obs2, color="darkgreen", lw=2, label="boundary 2 obs")
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


if __name__ == "__main__":
    sys.exit(main())
