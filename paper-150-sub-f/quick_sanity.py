"""
Quick sanity preview at N=1000 using existing Riemann zero cache.

This is NOT the official Sub F result (which requires N=2792 to match
the pre-registration). It's a quick check to make sure the pipeline
behaves reasonably while the main run computes 2792 targets.

Uses cached Riemann zeros (first 1000) plus Gram + RvM at the same N.
"""

import sys, csv, math, time
from pathlib import Path
import numpy as np
import mpmath as mp

SUB_C_DIR = Path(__file__).parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))
from probe import compute_probe, find_minima, match_to_targets, mc_null

mp.mp.dps = 15

N = 1000
D_PROBE = 5
N_TERMS = 5000
DT = 0.008
T_PAD = 5.0
WINDOW = 1.0
N_NULL = 1000
SEED = 42


def load_riemann(n=N):
    p = SUB_C_DIR / "riemann_zeros_10000.csv"
    zs = []
    with open(p, 'r') as f:
        r = csv.reader(f); next(r)
        for row in r:
            zs.append(float(row[1]))
    return np.array(zs[:n])


def gram_quick(n=N):
    """Quick Gram point computation using local-bracket findroot."""
    gs = []
    last = 17.85
    for k in range(1, n + 1):
        n_target = k - 1
        target = mp.mpf(n_target) * mp.pi
        if last > 6.5:
            step = 2.0 * math.pi / math.log(last / (2.0 * math.pi))
        else:
            step = 1.0
        try:
            g = float(mp.findroot(lambda t: mp.siegeltheta(t) - target, last + step))
        except Exception:
            g = last + step
        gs.append(g)
        last = g
    return np.array(gs)


def rvm_random(n, T_max, seed):
    rng = np.random.default_rng(seed)
    us = rng.uniform(0, 1, n)
    two_pi = 2.0 * math.pi
    F = lambda t: (t * math.log(t) - t + 1.0) / two_pi
    F_max = F(T_max)
    out = []
    for u in us:
        target = u * F_max
        lo, hi = 1.0, T_max
        for _ in range(60):
            mid = 0.5*(lo + hi)
            if F(mid) < target:
                lo = mid
            else:
                hi = mid
        out.append(0.5*(lo + hi))
    return np.array(sorted(out))


def main():
    print(f"Quick sanity at N={N}")
    riemann = load_riemann(N)
    print(f"  Loaded {len(riemann)} Riemann zeros, range [{riemann[0]:.2f}, {riemann[-1]:.2f}]")
    T_max = float(riemann[-1])

    print(f"  Computing {N} Gram points...")
    t0 = time.time()
    gram = gram_quick(N)
    print(f"  Gram in {time.time()-t0:.1f}s; range [{gram[0]:.2f}, {gram[-1]:.2f}]")

    rvm = rvm_random(N, T_max + T_PAD, seed=20260530)
    print(f"  RvM range [{rvm[0]:.2f}, {rvm[-1]:.2f}]")

    # Probe
    t_min = max(1.0, float(min(riemann[0], gram[0], rvm[0])) - T_PAD)
    t_max = max(float(riemann[-1]), float(gram[-1]), float(rvm[-1])) + T_PAD
    t_array = np.arange(t_min, t_max, DT)
    print(f"  Probe grid: {len(t_array)} pts")
    t0 = time.time()
    re, im, _ = compute_probe(t_array, d=D_PROBE, N=N_TERMS)
    print(f"  Probe in {time.time()-t0:.1f}s")
    probe = np.sqrt(np.abs(re*re - float(D_PROBE)*im*im))
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    print(f"  Minima: {len(minima_t)} (density {len(minima_t)/(t_max-t_min):.3f}/u)")

    # Match against each target
    for label, target in [("Riemann (control)", riemann), ("Gram", gram), ("RvM", rvm)]:
        tin = target[(target >= t_min) & (target <= t_max)]
        matches = match_to_targets(minima_t, tin, window=WINDOW)
        sig = float(np.mean([m['delta'] for m in matches]))
        null = mc_null(tin, n_minima=len(minima_t),
                       t_min=t_min, t_max=t_max, window=WINDOW,
                       n_trials=N_NULL, seed=SEED)
        z = (sig - null['mean']) / null['std']
        eff = (null['mean'] - sig) / null['mean']
        print(f"  [{label}] z = {z:+.3f}, eff = {eff:+.4f}, n_targets = {len(tin)}")

    # Anti-circularity: Gram-Riemann distances
    zs = np.sort(riemann)
    dists = []
    for g in gram:
        i = int(np.searchsorted(zs, g))
        cands = []
        if i < len(zs):
            cands.append(abs(zs[i] - g))
        if i > 0:
            cands.append(abs(zs[i-1] - g))
        dists.append(min(cands))
    dists = np.array(dists)
    print(f"\n  Gram-Riemann distance: median = {np.median(dists):.4f}, "
          f"10/90 pct = {np.percentile(dists,10):.4f}/{np.percentile(dists,90):.4f}")
    print(f"  Fraction within W={WINDOW}: {np.sum(dists<WINDOW)/len(dists)*100:.1f}%")
    print(f"  Fraction within W/2: {np.sum(dists<0.5*WINDOW)/len(dists)*100:.1f}%")


if __name__ == '__main__':
    main()
