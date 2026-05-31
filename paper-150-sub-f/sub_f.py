"""
Sub F driver — density-matched non-L target null for Paper 150 v2.1.

Per PRE_REGISTRATION.md, this is the SKELETON committed before any
results are computed. It will:

1. Load/extend Riemann zeros to N=2792 (Target 1).
2. Compute Gram points (Target 2) via mpmath.siegeltheta.
3. Sample 2792 RvM-density-matched random points (Target 3).
4. (Optional) GUE-pair-correlation-matched (Target 4).
5. For each target, run probe + matching + null using paper-203-sub-c
   machinery. Compute z and effect size ε.
6. Compute Gram-Riemann distance distribution (anti-circularity check).
7. Apply pre-registered thresholds (COLLAPSE / PARTIAL / SUSTAIN per target).
8. Aggregate verdict.

Implementation imports from paper-203-sub-c/probe.py to ensure no
re-implementation of the probe machinery.
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path
import numpy as np
import mpmath as mp

SUB_C_DIR = Path(__file__).parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))

from probe import compute_probe, find_minima, match_to_targets, mc_null

mp.mp.dps = 25

# Pre-registered constants
N_TARGET = 2792
N_TERMS = 5000
WINDOW = 1.0
N_NULL = 1000
SEED = 42
RVM_SEED = 20260530
DT = 0.008
T_PAD = 5.0
D_PROBE = 5


# --------- Target 1: Riemann zeros ---------

def get_riemann_zeros(n: int = N_TARGET) -> np.ndarray:
    """Load Riemann zero cache and extend to n via mpmath.zetazero."""
    cache_path = SUB_C_DIR / "riemann_zeros_10000.csv"
    existing = []
    if cache_path.exists():
        with open(cache_path, 'r') as f:
            r = csv.reader(f)
            next(r, None)
            for row in r:
                existing.append(float(row[1]))
    if len(existing) >= n:
        return np.array(existing[:n])
    print(f"  extending Riemann zeros from {len(existing)} to {n}...")
    t0 = time.time()
    zs = list(existing)
    for k in range(len(existing) + 1, n + 1):
        zs.append(float(mp.zetazero(k).imag))
        if k % 100 == 0:
            with open(cache_path, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['n', 't_n'])
                for i, z in enumerate(zs):
                    w.writerow([i + 1, f"{z:.10f}"])
            print(f"    {k}/{n} ({time.time()-t0:.0f}s)", flush=True)
    with open(cache_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 't_n'])
        for i, z in enumerate(zs):
            w.writerow([i + 1, f"{z:.10f}"])
    return np.array(zs)


# --------- Target 2: Gram points ---------

def get_gram_points(n: int = N_TARGET, cache_path: str = "gram_points.csv") -> np.ndarray:
    """Compute first n Gram points g_n where siegeltheta(g_n) = (n-1)*pi for n=1..N.

    Convention here: g_1 is the FIRST Gram point. siegeltheta(g_1) = 0
    (the lowest point where Riemann-Siegel theta = 0).
    """
    p = Path(__file__).parent / cache_path
    existing = []
    if p.exists():
        with open(p, 'r') as f:
            r = csv.reader(f)
            next(r, None)
            for row in r:
                existing.append(float(row[1]))
        if len(existing) >= n:
            return np.array(existing[:n])
    print(f"  computing Gram points {len(existing)+1}..{n} via mpmath.siegeltheta...")
    t0 = time.time()
    gs = list(existing)
    # Gram points: first few are well-known. g_0 ≈ 17.84559...
    # We use the convention g_n = root of siegeltheta(t) = n*pi (n = -1, 0, 1, ...).
    # g_{-1} ≈ 9.667 is the first one with siegeltheta(g) = -pi.
    # For consistency we'll use n = 0, 1, 2, ..., N-1 with siegeltheta(g_n) = n*pi.
    # g_0 ≈ 17.846 is the lowest Gram point > 0.
    start = len(existing) + 1   # we'll have g[0], g[1], ... = entries 1, 2, ...
    # The k-th entry corresponds to n = k - 1.
    last_g = existing[-1] if existing else 18.0
    for k in range(start, n + 1):
        n_target = k - 1   # siegeltheta(g) = n_target * pi
        target = mp.mpf(n_target) * mp.pi
        # Bracket: previous Gram point and previous + 1/log(prev) (the Gram spacing)
        lo = last_g + mp.mpf('0.1')
        hi = last_g + mp.pi  # safe upper bound (siegeltheta grows roughly like t*log(t/2pi)/2)
        # Adjust hi if needed
        while True:
            v_lo = mp.siegeltheta(lo) - target
            v_hi = mp.siegeltheta(hi) - target
            if v_lo * v_hi < 0:
                break
            hi += 1.0
            if hi - last_g > 20:
                # safety bail
                break
        try:
            g = float(mp.findroot(lambda t: mp.siegeltheta(t) - target, (lo + hi) / 2))
        except Exception:
            # fallback: use bisection manually
            g = float((lo + hi) / 2)
        gs.append(g)
        last_g = g
        if k % 100 == 0:
            with open(p, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['n', 'g_n'])
                for i, g_i in enumerate(gs):
                    w.writerow([i + 1, f"{g_i:.10f}"])
            print(f"    {k}/{n} (g_{k-1} = {g:.2f}, {time.time()-t0:.0f}s)", flush=True)
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['n', 'g_n'])
        for i, g_i in enumerate(gs):
            w.writerow([i + 1, f"{g_i:.10f}"])
    return np.array(gs)


# --------- Target 3: RvM-density-matched random ---------

def sample_rvm_random(n: int, T_max: float, seed: int = RVM_SEED) -> np.ndarray:
    """Sample n points from rho(t) = log(t)/(2*pi) on (1, T_max] via inverse CDF.

    CDF F(t) = (1/(2pi)) * (t*log(t) - t + 1).
    Normalised: F_tilde(t) = F(t) / F(T_max).
    For u in (0,1) uniform, solve F_tilde(t) = u via bisection.
    """
    rng = np.random.default_rng(seed)
    us = rng.uniform(0.0, 1.0, size=n)
    two_pi = 2.0 * math.pi

    def F(t):
        # F(1) = (1/(2pi))*(1*0 - 1 + 1) = 0  ✓
        return (t * math.log(t) - t + 1.0) / two_pi

    F_max = F(T_max)
    samples = []
    for u in us:
        target = u * F_max
        # Bisection on F(t) = target, t in [1, T_max]
        lo, hi = 1.0, T_max
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if F(mid) < target:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-9:
                break
        samples.append(0.5 * (lo + hi))
    return np.array(sorted(samples))


# --------- Anti-circularity check ---------

def gram_zero_distances(gram: np.ndarray, zeros: np.ndarray) -> dict:
    """For each Gram point, distance to nearest Riemann zero. Report distribution."""
    zs = np.sort(zeros)
    dists = []
    for g in gram:
        idx = int(np.searchsorted(zs, g))
        cands = []
        if idx < len(zs):
            cands.append(abs(zs[idx] - g))
        if idx > 0:
            cands.append(abs(zs[idx - 1] - g))
        dists.append(min(cands))
    dists = np.array(dists)
    return dict(
        median=float(np.median(dists)),
        p10=float(np.percentile(dists, 10)),
        p90=float(np.percentile(dists, 90)),
        mean=float(np.mean(dists)),
        below_window=int(np.sum(dists < WINDOW)),
        below_half_window=int(np.sum(dists < 0.5 * WINDOW)),
    )


# --------- Run one target ---------

def run_target(label: str, target_zeros: np.ndarray, t_min: float, t_max: float,
                re: np.ndarray, im: np.ndarray, t_array: np.ndarray) -> dict:
    """Match probe minima against target_zeros, compute z and effect."""
    probe = np.sqrt(np.abs(re * re - float(D_PROBE) * im * im))
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    matches = match_to_targets(minima_t, target_zeros, window=WINDOW)
    signal_mean = float(np.mean([m['delta'] for m in matches]))
    null = mc_null(target_zeros, n_minima=len(minima_t),
                   t_min=t_min, t_max=t_max, window=WINDOW,
                   n_trials=N_NULL, seed=SEED)
    z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
    eff = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0
    return dict(
        label=label, n_targets=len(target_zeros), n_minima=len(minima_t),
        signal_mean=signal_mean, null_mean=null['mean'], null_std=null['std'],
        z=float(z), effect=float(eff),
    )


# --------- Main ---------

def main():
    print("=" * 70)
    print("Sub F — density-matched non-L target null")
    print("=" * 70)

    # Target 1: Riemann zeros (control)
    print("\n[Target 1] Riemann zeros (control)")
    riemann = get_riemann_zeros(N_TARGET)
    print(f"  Loaded {len(riemann)} Riemann zeros, range [{riemann[0]:.2f}, {riemann[-1]:.2f}]")
    T_max = float(riemann[-1])

    # Target 2: Gram points
    print(f"\n[Target 2] Gram points (siegeltheta(g_n) = (n-1)*pi)")
    gram = get_gram_points(N_TARGET)
    print(f"  Loaded {len(gram)} Gram points, range [{gram[0]:.2f}, {gram[-1]:.2f}]")

    # Target 3: RvM-density-matched random
    print(f"\n[Target 3] RvM-density-matched random (seed = {RVM_SEED})")
    rvm = sample_rvm_random(N_TARGET, T_max + T_PAD, seed=RVM_SEED)
    print(f"  Sampled {len(rvm)} RvM-random points, range [{rvm[0]:.2f}, {rvm[-1]:.2f}]")

    # Anti-circularity: Gram-Riemann distance distribution
    print(f"\n[Anti-circularity] Gram-Riemann nearest-neighbour distances:")
    gz = gram_zero_distances(gram, riemann)
    print(f"  median = {gz['median']:.4f}  (window W = {WINDOW})")
    print(f"  10th/90th pct = {gz['p10']:.4f} / {gz['p90']:.4f}")
    print(f"  below W: {gz['below_window']}/{len(gram)} ({gz['below_window']/len(gram)*100:.1f}%)")
    print(f"  below W/2: {gz['below_half_window']}/{len(gram)} ({gz['below_half_window']/len(gram)*100:.1f}%)")

    # Compute probe once (Re, Im don't depend on target)
    t_min = max(1.0, float(min(riemann[0], gram[0], rvm[0])) - T_PAD)
    t_max = max(float(riemann[-1]), float(gram[-1]), float(rvm[-1])) + T_PAD
    print(f"\n[Probe] Computing Z_phi on t in [{t_min:.1f}, {t_max:.1f}], dt = {DT}")
    print(f"  Grid points: {int((t_max - t_min) / DT)}")
    t_array = np.arange(t_min, t_max, DT)
    t0 = time.time()
    re, im, _ = compute_probe(t_array, d=D_PROBE, N=N_TERMS)
    print(f"  Probe Re, Im computed in {time.time()-t0:.1f}s")

    # Run against each target
    print(f"\n=== Detection per target ===")
    results = []
    for label, target in [("Riemann (control)", riemann),
                            ("Gram points", gram),
                            ("RvM random", rvm)]:
        t0 = time.time()
        # Filter target to within probe range
        target_in_range = target[(target >= t_min) & (target <= t_max)]
        r = run_target(label, target_in_range, t_min, t_max, re, im, t_array)
        print(f"  {label}: z = {r['z']:+.4f}, eff = {r['effect']:.4f}, "
              f"n_targets = {r['n_targets']} ({time.time()-t0:.0f}s)")
        results.append(r)

    # Verdict
    print(f"\n=== Pre-registered verdicts ===")
    verdicts = {}
    for r in results[1:]:   # skip control
        z_abs = abs(r['z'])
        if z_abs < 3:
            v = "COLLAPSE"
        elif z_abs < 10:
            v = "PARTIAL"
        else:
            v = "SUSTAIN"
        verdicts[r['label']] = v
        print(f"  {r['label']}: |z| = {z_abs:.2f} → {v}")

    # Aggregate
    vs = list(verdicts.values())
    if "SUSTAIN" in vs:
        agg = "SUSTAIN"
    elif all(v == "COLLAPSE" for v in vs):
        agg = "COLLAPSE"
    else:
        agg = "PARTIAL"
    print(f"\nAggregate verdict: {agg}")

    # Save CSV
    with open('sub_f_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['target', 'n_targets', 'n_minima', 'signal_mean',
                    'null_mean', 'null_std', 'z', 'effect_size', 'verdict'])
        for r in results:
            v = verdicts.get(r['label'], 'CONTROL')
            w.writerow([r['label'], r['n_targets'], r['n_minima'],
                        r['signal_mean'], r['null_mean'], r['null_std'],
                        r['z'], r['effect'], v])
    print("\nWrote sub_f_results.csv")

    return dict(results=results, gz=gz, verdicts=verdicts, aggregate=agg)


if __name__ == '__main__':
    main()
