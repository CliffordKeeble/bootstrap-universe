"""
Sub G — gap-shuffle null on Riemann zeros.

Reuses Sub F follow-up machinery; only the target type changes.
Per PRE_REGISTRATION.md, committed before any gap-shuffle is performed.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path
import numpy as np

SUB_F_DIR = Path(__file__).parent.parent / "paper-150-sub-f"
SUB_C_DIR = Path(__file__).parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))
sys.path.insert(0, str(SUB_F_DIR))
from probe import compute_probe, find_minima, match_to_targets, mc_null


# Pre-registered constants (frozen, identical to Sub F follow-up)
N_TARGET = 2792
N_TERMS = 5000
WINDOW = 1.0
N_NULL = 1000
NULL_SEED = 42
DT = 0.008
T_PAD = 5.0
D_PROBE = 5

# Sub F cached numbers (for comparison)
Z_RIEMANN = -12.5629
EFF_RIEMANN = 0.2335
NULL_STD_RIEMANN = 0.0079
NULL_STD_RVM_MEAN = 0.0105

# 10 pre-registered gap-shuffle seeds (same as Sub F follow-up)
SHUFFLE_SEEDS = [20260530, 20260531, 20260601, 20260602, 20260603,
                 20260604, 20260605, 20260606, 20260607, 20260608]


def load_riemann_first_n(n: int = N_TARGET) -> np.ndarray:
    p = SUB_C_DIR / "riemann_zeros_10000.csv"
    zs = []
    with open(p, 'r') as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            zs.append(float(row[1]))
    return np.array(zs[:n])


def gap_shuffle(zeros: np.ndarray, seed: int) -> np.ndarray:
    """Shuffle the gap order of `zeros` using a seeded RNG.

    Result: same first point, same last point, same gap multiset,
    randomised gap order.
    """
    gaps = np.diff(zeros)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(gaps))
    shuffled_gaps = gaps[perm]
    return np.concatenate([[zeros[0]], zeros[0] + np.cumsum(shuffled_gaps)])


def main():
    print("=" * 70)
    print("Sub G — target gap-shuffle null on Riemann zeros")
    print(f"Riemann control (Sub F): z = {Z_RIEMANN:.4f}, "
          f"eff = {EFF_RIEMANN:.4f}, null_std = {NULL_STD_RIEMANN:.4f}")
    print(f"RvM-Poisson (Sub F follow-up mean): null_std = {NULL_STD_RVM_MEAN:.4f}")
    print("=" * 70)

    riemann = load_riemann_first_n(N_TARGET)
    T_max = float(riemann[-1])
    print(f"\nRiemann cache: {len(riemann)} zeros, T_max = {T_max:.2f}")

    # Verify gap-shuffle preserves first/last point and gap multiset
    test_shuf = gap_shuffle(riemann, seed=SHUFFLE_SEEDS[0])
    assert abs(test_shuf[0] - riemann[0]) < 1e-9, "first point should match"
    assert abs(test_shuf[-1] - riemann[-1]) < 1e-9, "last point should match"
    orig_gaps = np.sort(np.diff(riemann))
    shuf_gaps = np.sort(np.diff(test_shuf))
    assert np.allclose(orig_gaps, shuf_gaps), "gap multiset should match"
    print("Gap-shuffle verified: preserves start, end, and gap multiset.")

    # Compute probe Re/Im once
    t_min = 1.0
    t_max_grid = T_max + 10.0
    t_array = np.arange(t_min, t_max_grid, DT)
    print(f"\nProbe grid: t in [{t_min}, {t_max_grid:.1f}], {len(t_array)} pts")
    t0 = time.time()
    re, im, _ = compute_probe(t_array, d=D_PROBE, N=N_TERMS)
    print(f"Probe Re/Im in {time.time()-t0:.1f}s")
    probe = np.sqrt(np.abs(re * re - float(D_PROBE) * im * im))
    minima_t, _ = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    print(f"Minima: {len(minima_t)} (density {len(minima_t)/(t_max_grid-t_min):.3f}/u)")

    rows = []
    print(f"\nRunning {len(SHUFFLE_SEEDS)} gap-shuffle seeds...")
    for seed in SHUFFLE_SEEDS:
        t0 = time.time()
        target = gap_shuffle(riemann, seed)
        target_in = target[(target >= t_min) & (target <= t_max_grid)]
        matches = match_to_targets(minima_t, target_in, window=WINDOW)
        sig = float(np.mean([m['delta'] for m in matches]))
        null = mc_null(target_in, n_minima=len(minima_t),
                       t_min=t_min, t_max=t_max_grid, window=WINDOW,
                       n_trials=N_NULL, seed=NULL_SEED)
        z = (sig - null['mean']) / null['std'] if null['std'] > 0 else 0.0
        eff = (null['mean'] - sig) / null['mean'] if null['mean'] > 0 else 0.0
        n_matched = sum(1 for m in matches if m['matched'])
        row = dict(seed=seed, n_targets=len(target_in), n_matched=n_matched,
                   signal_mean=sig, null_mean=null['mean'], null_std=null['std'],
                   z=z, effect=eff)
        rows.append(row)
        print(f"  seed={seed}: z = {z:+.4f}, eff = {eff:.4f}, "
              f"null_std = {null['std']:.4f}, n_matched = {n_matched} "
              f"({time.time()-t0:.1f}s)")

    # Aggregate
    zs = np.array([r['z'] for r in rows])
    effs = np.array([r['effect'] for r in rows])
    null_stds = np.array([r['null_std'] for r in rows])
    abs_zs = np.abs(zs)

    print(f"\n=== z distribution (10 shuffles) ===")
    print(f"  mean = {zs.mean():+.4f}")
    print(f"  std  = {zs.std():+.4f}")
    print(f"  min  = {zs.min():+.4f}")
    print(f"  max  = {zs.max():+.4f}")
    print(f"  median = {np.median(zs):+.4f}")

    print(f"\n=== eff distribution (10 shuffles) ===")
    print(f"  mean = {effs.mean():.4f}")
    print(f"  std  = {effs.std():.4f}")

    print(f"\n=== null_std distribution (10 shuffles) ===")
    print(f"  mean = {null_stds.mean():.4f}")
    print(f"  std  = {null_stds.std():.4f}")
    print(f"  min  = {null_stds.min():.4f}")
    print(f"  max  = {null_stds.max():.4f}")
    print(f"  vs Riemann null_std = {NULL_STD_RIEMANN:.4f}")
    print(f"  vs RvM null_std mean = {NULL_STD_RVM_MEAN:.4f}")

    # K_shuffle
    K = (abs(Z_RIEMANN) - abs_zs.mean()) / abs_zs.std() if abs_zs.std() > 0 else float('inf')
    print(f"\n=== K_shuffle ===")
    print(f"  K = ({abs(Z_RIEMANN):.2f} - {abs_zs.mean():.2f}) / {abs_zs.std():.3f} = {K:.3f}")

    # Verdict per pre-reg
    null_std_mean = float(null_stds.mean())
    h1_null = abs(null_std_mean - NULL_STD_RIEMANN) / NULL_STD_RIEMANN <= 0.15
    h2_null = abs(null_std_mean - NULL_STD_RVM_MEAN) / NULL_STD_RVM_MEAN <= 0.30

    print(f"\n=== Verdict criteria ===")
    print(f"  K < 1.5: {K < 1.5}")
    print(f"  null_std within 15% of Riemann 0.0079 ({h1_null}): {null_std_mean:.4f}")
    print(f"  K > 3: {K > 3}")
    print(f"  null_std within 30% of RvM 0.0105 ({h2_null}): {null_std_mean:.4f}")

    if K < 1.5 and h1_null:
        verdict = "NEAREST-NEIGHBOR SUFFICIENT (H1)"
    elif K > 3 and h2_null:
        verdict = "LONGER-RANGE MATTERS (H2)"
    else:
        verdict = "PARTIAL (H1+H2 mixed)"
    print(f"\nVerdict: {verdict}")

    with open('sub_g_target_shuffle_results.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['seed', 'n_targets', 'n_matched', 'signal_mean',
                    'null_mean', 'null_std', 'z', 'effect_size'])
        for r in rows:
            w.writerow([r['seed'], r['n_targets'], r['n_matched'],
                        r['signal_mean'], r['null_mean'], r['null_std'],
                        r['z'], r['effect']])
    print("Wrote sub_g_target_shuffle_results.csv")


if __name__ == '__main__':
    main()
