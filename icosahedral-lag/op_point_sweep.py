"""
Pre-data operating-point sweep per pre-registration: small N=30, T=2000.
Find (eps, tau) where dynamics reaches a differentiated steady state.

Criteria (pre-registered):
(a) mean inter-unit correlation in (0.05, 0.95)
(b) std of correlation distribution high (most differentiated)
(c) first-half-vs-second-half mean correlation difference < 0.05 (stability)
"""

from __future__ import annotations

import csv
import numpy as np
from dynamics import RunConfig, run_dynamics

EPS_GRID = [0.025, 0.05, 0.1, 0.2, 0.4, 0.6]
TAU_GRID = [0, 1, 2, 3, 4, 5, 7, 10, 20, 50]
N_SWEEP = 30
T_SWEEP = 2000
TRANSIENT_SWEEP = 500


def quick_stability(record: np.ndarray) -> float:
    """Compare first-half vs second-half mean correlation."""
    n_steps = record.shape[0]
    half = n_steps // 2

    def mean_off(M):
        iu = np.triu_indices(M.shape[0], k=1)
        return float(np.mean(M[iu]))

    C1 = np.zeros((record.shape[1], record.shape[1]))
    for t in range(half):
        C1 += record[t] @ record[t].T
    C1 /= half
    C2 = np.zeros_like(C1)
    for t in range(half, n_steps):
        C2 += record[t] @ record[t].T
    C2 /= (n_steps - half)
    return abs(mean_off(C2) - mean_off(C1))


def sweep():
    rows = []
    for eps in EPS_GRID:
        for tau in TAU_GRID:
            cfg = RunConfig(N=N_SWEEP, T=T_SWEEP, transient=TRANSIENT_SWEEP,
                            eps=eps, tau=tau, seed=20260528)
            res = run_dynamics(cfg)
            m = res['meta']
            stab = quick_stability(res['record'])
            ok_diff = 0.05 < m['mean_corr'] < 0.95 or -0.95 < m['mean_corr'] < -0.05 or \
                      (abs(m['mean_corr']) < 0.05 and m['std_corr'] > 0.2)
            ok_stab = stab < 0.05
            rows.append({
                'eps': eps, 'tau': tau,
                'mean_corr': m['mean_corr'], 'std_corr': m['std_corr'],
                'min_corr': m['min_corr'], 'max_corr': m['max_corr'],
                'stability_diff': stab,
                'differentiated': ok_diff,
                'stable': ok_stab,
                'criterion_C': m['std_corr'] if (ok_diff and ok_stab) else 0.0,
            })
            print(f"  eps={eps:.3f} tau={tau:2d}: mean={m['mean_corr']:+.3f} "
                  f"std={m['std_corr']:.3f} stab={stab:.4f} diff={ok_diff} "
                  f"stable={ok_stab}")

    # Pick winner: maximises std_corr among (differentiated AND stable)
    valid = [r for r in rows if r['differentiated'] and r['stable']]
    if not valid:
        print("No operating point satisfies both criteria. Need wider grid.")
        best = max(rows, key=lambda r: r['criterion_C'])
    else:
        best = max(valid, key=lambda r: r['std_corr'])
    print(f"\nSelected operating point: eps={best['eps']} tau={best['tau']} "
          f"std_corr={best['std_corr']:.3f}")
    # Write CSV
    with open('op_point_sweep.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return best


if __name__ == '__main__':
    best = sweep()
    print(f"\nWinner: {best}")
