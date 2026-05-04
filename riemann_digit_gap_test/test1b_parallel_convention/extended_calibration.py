# extended_calibration.py - 9-point Selberg-null calibration for Test 1b
# Bootstrap Universe Programme - Brief 4 May 2026 (CinC), Test 1b
#
# Mr Adversary recommendation upgraded the 5-point audit (Test 1) to a
# 9-point audit covering the broader test region of Test 1b. Five points
# from Test 1 are subsumed (60, 80, 100, 120, 160); four new points at
# 145, 195, 210, 230 extend coverage into the higher-T regime where
# Test 1b's boundary 2 + control region C4 sit.
#
# Pass conditions (pre-registered):
#   All 9 |R(T)| < 1.5                 (Test 1 envelope)
#   Mean R within +/- 0.2              (drift indicator, CinC tightened)
#   High-T points (T = 195, 210, 230): |R| < 1.0
#       -> CinC: "S(T) fluctuation scale grows slowly with T;
#          residuals should stay within roughly +/- 1.0 even at T = 230.
#          If they're systematically larger, the Selberg null may need a
#          higher-order correction term and that's a structural finding
#          worth pausing for, not improvising past."
#
# Output:
#   stdout: audit table + verdict
#   extended_calibration_results.json (machine-readable record)

import csv
import json
import sys
from pathlib import Path

import mpmath
from mpmath import mp, mpf

# Reuse Test 1's Selberg-null functions verbatim. The pre-registered
# convention is that calibration uses the exact same Delta_smooth(T) as
# the windowed test - so importing from the parent zeros_data.py keeps
# the calibration honest. CinC: "DRY-ing them across Test 1 and Test 1b
# strengthens the audit by ensuring identical implementation."
_RDGT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_RDGT))
from zeros_data import Delta_smooth, DEFAULT_DPS

CACHE_DIR = _RDGT / "zeros_cache"
ZETA_CSV = CACHE_DIR / "riemann_zeros_dps30_n100.csv"
LCHI5_CSV = CACHE_DIR / "lchi5_zeros_dps30_n180.csv"

# Pre-registered audit points.
AUDIT_T = [60, 80, 100, 120, 145, 160, 195, 210, 230]
HIGH_T  = {195, 210, 230}      # CinC pause-and-report set

ENVELOPE_BOUND = 1.5
HIGH_T_BOUND = 1.0
DRIFT_BOUND = 0.2


def load_zeros_csv(path):
    zeros = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            zeros.append(mpf(row["gamma"]))
    return zeros


def main():
    mp.dps = DEFAULT_DPS + 5

    print("extended_calibration - 9-point Selberg null audit for Test 1b")
    print("=" * 76)
    print("Pass conditions (pre-registered):")
    print("  All |R| < %.1f  (Test 1 envelope)" % ENVELOPE_BOUND)
    print("  Mean R within +/- %.1f  (drift indicator)" % DRIFT_BOUND)
    print("  High-T |R| (T in %s) < %.1f"
          % (sorted(HIGH_T), HIGH_T_BOUND))
    print("    (S(T) grows slowly; if high-T residuals exceed this,")
    print("     pause and report - may indicate higher-order correction needed)")
    print()

    print("Loading zero data...")
    if not ZETA_CSV.exists():
        print("  ERROR: %s not found." % ZETA_CSV)
        return 1
    if not LCHI5_CSV.exists():
        print("  ERROR: %s not found. Run extend_lchi5_sweep.py first."
              % LCHI5_CSV)
        return 1
    zeta_zeros = load_zeros_csv(ZETA_CSV)
    lchi5_zeros = load_zeros_csv(LCHI5_CSV)
    print("  zeta zeros:    n=%d, max gamma=%s"
          % (len(zeta_zeros), mpmath.nstr(zeta_zeros[-1], 8)))
    print("  L(chi_5) zeros: n=%d, max gamma=%s"
          % (len(lchi5_zeros), mpmath.nstr(lchi5_zeros[-1], 8)))
    print()

    # Compute audit table.
    print("Audit at non-boundary heights:")
    print()
    print("  | T   | N_zeta | N_chi5 | Delta_emp | Delta_pred  | R       | tag       |")
    print("  |-----|--------|--------|-----------|-------------|---------|-----------|")
    rows = []
    for T in AUDIT_T:
        n_z = sum(1 for g in zeta_zeros if g <= T)
        n_l = sum(1 for g in lchi5_zeros if g <= T)
        delta_emp = n_z - n_l
        delta_pred = float(Delta_smooth(T))
        R = delta_emp - delta_pred
        tag = "high-T" if T in HIGH_T else "low-T"
        print("  | %3d | %6d | %6d | %9d | %11.4f | %+7.4f | %-9s |"
              % (T, n_z, n_l, delta_emp, delta_pred, R, tag))
        rows.append({
            "T": T, "n_zeta": n_z, "n_chi5": n_l,
            "delta_emp": delta_emp, "delta_pred": float(delta_pred),
            "R": R, "is_high_T": T in HIGH_T,
        })
    print()

    # Verdict.
    R_values = [r["R"] for r in rows]
    R_mean = sum(R_values) / len(R_values)
    R_min, R_max = min(R_values), max(R_values)
    high_T_R = [r["R"] for r in rows if r["is_high_T"]]

    envelope_pass = all(abs(r["R"]) < ENVELOPE_BOUND for r in rows)
    drift_pass = abs(R_mean) < DRIFT_BOUND
    high_T_pass = all(abs(R) < HIGH_T_BOUND for R in high_T_R)

    print("Summary:")
    print("  R mean (drift):        %+.4f   (bound: +/- %.2f)   %s"
          % (R_mean, DRIFT_BOUND, "PASS" if drift_pass else "FAIL"))
    print("  R range:               [%+.3f, %+.3f]   (envelope: +/- %.2f)   %s"
          % (R_min, R_max, ENVELOPE_BOUND, "PASS" if envelope_pass else "FAIL"))
    print("  High-T |R| range:      [%.3f, %.3f]    (bound: %.2f)         %s"
          % (min(abs(R) for R in high_T_R), max(abs(R) for R in high_T_R),
             HIGH_T_BOUND, "PASS" if high_T_pass else "FAIL  *PAUSE*"))
    print()

    overall = envelope_pass and drift_pass and high_T_pass
    print("=" * 76)
    if overall:
        print("CALIBRATION: PASS")
        print("  Linear-plus-constant Selberg null is calibrated through T = 230.")
        print("  Safe to proceed to pre-registration commit and Test 1b run.")
    else:
        print("CALIBRATION: FAIL  -  PAUSE AND REPORT")
        print("  Do not proceed to pre-registration. Selberg null may need a")
        print("  higher-order correction term, or there is a defect in the data.")

    # Write JSON.
    out = {
        "audit_T": AUDIT_T,
        "high_T_set": sorted(HIGH_T),
        "envelope_bound": ENVELOPE_BOUND,
        "high_T_bound": HIGH_T_BOUND,
        "drift_bound": DRIFT_BOUND,
        "rows": rows,
        "summary": {
            "R_mean": R_mean,
            "R_min": R_min,
            "R_max": R_max,
            "high_T_R": high_T_R,
            "envelope_pass": envelope_pass,
            "drift_pass": drift_pass,
            "high_T_pass": high_T_pass,
            "overall": overall,
        },
        "data": {
            "zeta_n": len(zeta_zeros),
            "zeta_max_gamma": float(zeta_zeros[-1]),
            "lchi5_n": len(lchi5_zeros),
            "lchi5_max_gamma": float(lchi5_zeros[-1]),
        },
    }
    out_path = Path(__file__).parent / "extended_calibration_results.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print()
    print("JSON written: %s" % out_path)

    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
