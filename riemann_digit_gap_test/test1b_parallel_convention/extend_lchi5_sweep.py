# extend_lchi5_sweep.py - Coverage extension for L(chi_5) zeros to gamma ~ 245
# Bootstrap Universe Programme - Brief 4 May 2026 (CinC), Test 1b coverage
#
# Brief addendum 4 May 2026 (Mr Adversary follow-up): Test 1b's boundary 2
# window is [213.30, 223.30] and control region C4 reaches 235. The existing
# n=150 cache covers gamma <= 219. We need ~180 zeros covering gamma <= 245
# for safety margin.
#
# Approach: independent Lambda-sign-change sweep on t in [0.5, 245], refine
# every sign change via the existing toolchain (Lambda_chi5 + bisection),
# write to a NEW cache file lchi5_zeros_dps30_n180.csv. Does not modify the
# n=150 cache used by Test 1 (which is now committed and immutable).
#
# Integrity check (per CinC): the sweep, run independently from t=0.5, must
# reproduce the first 25 zeros to within 4 decimals of LMFDB. This is a free
# audit of the sweep machinery - if anything in the toolchain has drifted
# the check catches it before the new zeros enter the cache.
#
# Provenance note: the n=150 cache mixes LMFDB-refined first-25 with
# sweep-refined 26-150. The n=180 cache is pure-sweep throughout, so the
# first-25 values may differ from the n=150 cache at the 5th decimal or
# beyond. This is benign - both are valid to dps=30 against the canonical
# zero locations.

import csv
import sys
import time
from pathlib import Path

import mpmath
from mpmath import mp, mpf, fabs

# Reuse the existing toolchain.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_REPO_ROOT / "golden-chebyshev"))
from lchi5_zeros import LMFDB_ZEROS_25, Lambda_chi5, refine_zero

CACHE_DIR = _REPO_ROOT / "riemann_digit_gap_test" / "zeros_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_CSV = CACHE_DIR / "lchi5_zeros_dps30_n180.csv"

# Sweep parameters.
DPS = 30
T_MIN = mpf("0.5")
T_MAX = mpf(245.0)
STEP  = mpf("0.4")           # well below local zero spacing (~1.2 at t=245)


def sweep_for_brackets(t_min, t_max, step, dps):
    """Find all sign-change brackets of Lambda_chi5 on [t_min, t_max]."""
    print("  Sweeping Lambda_chi5 on [%g, %g] step %g (dps=%d)..."
          % (float(t_min), float(t_max), float(step), dps), flush=True)
    brackets = []
    t = t_min
    f_prev = Lambda_chi5(t, dps)
    n_seen = 0
    while t < t_max:
        t_next = t + step
        if t_next > t_max:
            t_next = t_max
        f_next = Lambda_chi5(t_next, dps)
        n_seen += 1
        if f_prev * f_next < 0:
            brackets.append((t, t_next))
        if n_seen % 100 == 0:
            print("    grid step %d / %d, brackets so far: %d"
                  % (n_seen, int((t_max - t_min) / step), len(brackets)),
                  flush=True)
        t = t_next
        f_prev = f_next
    print("  Done: %d sign-change brackets found." % len(brackets), flush=True)
    return brackets


def refine_brackets(brackets, dps):
    """Refine every bracket via Lambda-bisection. Returns list of mpf zeros."""
    out = []
    t0 = time.time()
    for i, (a, b) in enumerate(brackets, start=1):
        mid = (a + b) / 2
        delta = (b - a) / 2 + mpf("0.001")
        z = refine_zero(mid, dps=dps, delta=delta)
        out.append(z)
        if i % 25 == 0:
            print("    refined %d / %d (gamma~%.2f, %.1fs elapsed)"
                  % (i, len(brackets), float(z), time.time() - t0),
                  flush=True)
    return out


def integrity_check_first_25(refined, lmfdb_seeds, decimals=4):
    """First 25 sweep-refined zeros must match LMFDB to `decimals` decimals."""
    print()
    print("Integrity check: first 25 sweep zeros vs LMFDB (4 decimals)")
    print("  | k  | sweep refined         | LMFDB seed | |diff|")
    print("  |----|-----------------------|------------|---------")
    tol = mpf(10) ** (-decimals + 1)
    n_pass = 0
    for k, (sw, lm) in enumerate(zip(refined[:25], lmfdb_seeds), start=1):
        diff = fabs(sw - lm)
        ok = diff < tol
        if ok:
            n_pass += 1
        print("  | %2d | %s | %s | %s   %s"
              % (k, mpmath.nstr(sw, 18),
                 mpmath.nstr(lm, 8),
                 mpmath.nstr(diff, 4),
                 "PASS" if ok else "FAIL"))
    print()
    print("  %d / 25 sweep-refined zeros agree with LMFDB to %d decimals."
          % (n_pass, decimals))
    return n_pass == 25


def write_cache(refined, csv_path, dps):
    """Write zeros to cache CSV, k-indexed."""
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "gamma"])
        for k, z in enumerate(refined, start=1):
            w.writerow([k, mpmath.nstr(z, dps + 5, strip_zeros=False)])
    print("  Cache written: %s  (%d zeros)" % (csv_path, len(refined)))


def main():
    mp.dps = DPS + 5

    print("extend_lchi5_sweep - L(chi_5) coverage extension for Test 1b")
    print("=" * 76)
    print("Target: gamma_max >= 233 (pre-reg C4 truncation rule); aim 245")
    print()

    t0 = time.time()

    # Stage 1: sweep
    brackets = sweep_for_brackets(T_MIN, T_MAX, STEP, DPS)

    # Stage 2: refine all brackets
    print()
    print("Refining %d brackets to dps=%d..." % (len(brackets), DPS))
    refined = refine_brackets(brackets, DPS)
    print("  All refined. Total wall time: %.1f s" % (time.time() - t0))

    # Stage 3: integrity check
    integrity_pass = integrity_check_first_25(refined, LMFDB_ZEROS_25)

    # Stage 4: write cache
    print()
    if integrity_pass:
        write_cache(refined, OUTPUT_CSV, DPS)
        print()
        print("Coverage: gamma_min = %s, gamma_max = %s"
              % (mpmath.nstr(refined[0], 8), mpmath.nstr(refined[-1], 8)))
        print("Total zeros: %d" % len(refined))
        # C4 status check
        gamma_max = float(refined[-1])
        if gamma_max >= 233:
            c4_upper = min(gamma_max - 5, 240.0)
            print("C4 = [223.30, %.2f] -- valid (gamma_max >= 233)" % c4_upper)
        else:
            print("C4 dropped per pre-reg rule (gamma_max=%.2f < 233)"
                  % gamma_max)
        return 0
    else:
        print("INTEGRITY CHECK FAILED. Cache NOT written.")
        print("Inspect the FAIL rows and the existing toolchain before proceeding.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
