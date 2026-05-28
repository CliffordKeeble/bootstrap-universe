"""
v0.1 verdict run: icosahedral+lag dynamics at fixed operating point
(eps=0.025, tau=3) from the pre-data sweep, plus the three pre-registered
nulls.

Computes Observables 1 (spectral dimension), 2 (volume saturation), 3
(triangle curvature), and 4 (spectral gap, supplementary).
"""

from __future__ import annotations

import csv
import time
import numpy as np

from dynamics import (
    RunConfig, run_dynamics, correlation_to_distance,
    weighted_graph_from_corr, spectral_dimension, volume_growth,
    fit_volume_dim_and_saturation, mean_triangle_defect,
    C_2I_DRIFT, PHI,
)

# Pre-registered operating point (from op_point_sweep.py output)
EPS = 0.025
TAU = 3
N = 240
T = 8000
TRANSIENT = 2000
SEED = 20260528


def analyse_corr(C: np.ndarray, label: str, n_total: int) -> dict:
    """Compute all observables on a correlation matrix."""
    D = correlation_to_distance(C)
    W, L = weighted_graph_from_corr(C, drop_negative=True)

    # Sanity: nonzero off-diagonal entries
    iu = np.triu_indices(n_total, k=1)
    n_pos_edges = int(np.sum(W[iu] > 0))

    # Spectral dimension
    try:
        d_s, ts, p_t = spectral_dimension(L)
    except Exception as e:
        d_s, ts, p_t = float('nan'), None, None

    # Volume growth
    r_grid, Nr = volume_growth(D)
    vol = fit_volume_dim_and_saturation(r_grid, Nr, n_total)

    # Triangle curvature
    tri = mean_triangle_defect(D, n_samples=5000, seed=SEED)

    # Observable 4: spectral gap (only meaningful if O1-O3 broadly pass)
    evals = np.linalg.eigvalsh(L)
    evals_sorted = np.sort(evals)
    lam1_raw = float(evals_sorted[1]) if len(evals_sorted) > 1 else float('nan')
    lam1_x_N23 = lam1_raw * (n_total ** (2.0 / 3.0))

    return dict(
        label=label,
        N=n_total,
        n_pos_edges=n_pos_edges,
        mean_corr=float(np.mean(C[iu])),
        std_corr=float(np.std(C[iu])),
        min_corr=float(np.min(C[iu])),
        max_corr=float(np.max(C[iu])),
        d_s=d_s,
        slope_volume=vol['slope'],
        saturated=vol['saturated'],
        sat_frac=vol['sat_frac'],
        mean_defect=tri['mean_defect'],
        frac_pos_defect=tri['frac_positive'],
        lam1_raw=lam1_raw,
        lam1_x_N23=lam1_x_N23,
    )


def shuffled_K(N: int, rng: np.random.Generator) -> np.ndarray:
    """Per pre-reg: K_ij ~ Uniform[0, 2] (mean 1) then shuffled per pair."""
    K = rng.uniform(0.0, 2.0, size=(N, N))
    # Symmetrise
    K = 0.5 * (K + K.T)
    np.fill_diagonal(K, 0.0)
    return K


def run_all() -> list[dict]:
    rows = []

    # ---- Icosahedral + lag (PRIMARY) ----
    print("Running icosahedral + lag (primary)...")
    t0 = time.time()
    cfg = RunConfig(N=N, T=T, transient=TRANSIENT, eps=EPS, tau=TAU,
                    drift_q=C_2I_DRIFT, K_mode="all-to-all", seed=SEED)
    res = run_dynamics(cfg)
    print(f"  ran in {time.time()-t0:.1f}s; mean_corr={res['meta']['mean_corr']:.4f} "
          f"std={res['meta']['std_corr']:.4f}")
    a = analyse_corr(res['corr'], "icosahedral+lag", N)
    rows.append(a)
    print(f"  d_s={a['d_s']:.3f} sat_frac={a['sat_frac']:.3f} sat={a['saturated']} "
          f"slope_vol={a['slope_volume']:.3f} frac_pos_def={a['frac_pos_defect']:.3f} "
          f"lam1*N^(2/3)={a['lam1_x_N23']:.3f}")

    # ---- Null 1: non-icosahedral drift ----
    print("\nRunning null 1: non-icosahedral drift...")
    t0 = time.time()
    # Irrational angle, axis along (1, 0, 0, 0) (no special symmetry)
    nonico_angle_half = np.pi / 3.7   # irrational with respect to A_5
    c_nonico = np.array([np.cos(nonico_angle_half),
                         np.sin(nonico_angle_half), 0.0, 0.0], dtype=np.float64)
    c_nonico /= np.linalg.norm(c_nonico)
    cfg = RunConfig(N=N, T=T, transient=TRANSIENT, eps=EPS, tau=TAU,
                    drift_q=c_nonico, K_mode="all-to-all", seed=SEED)
    res = run_dynamics(cfg)
    print(f"  ran in {time.time()-t0:.1f}s; mean_corr={res['meta']['mean_corr']:.4f} "
          f"std={res['meta']['std_corr']:.4f}")
    a = analyse_corr(res['corr'], "null:non-icosahedral", N)
    rows.append(a)
    print(f"  d_s={a['d_s']:.3f} sat_frac={a['sat_frac']:.3f} sat={a['saturated']}")

    # ---- Null 2: no delay ----
    print("\nRunning null 2: no delay (tau=0)...")
    t0 = time.time()
    cfg = RunConfig(N=N, T=T, transient=TRANSIENT, eps=EPS, tau=0,
                    drift_q=C_2I_DRIFT, K_mode="all-to-all", seed=SEED)
    res = run_dynamics(cfg)
    print(f"  ran in {time.time()-t0:.1f}s; mean_corr={res['meta']['mean_corr']:.4f} "
          f"std={res['meta']['std_corr']:.4f}")
    a = analyse_corr(res['corr'], "null:no-delay", N)
    rows.append(a)
    print(f"  d_s={a['d_s']:.3f} sat_frac={a['sat_frac']:.3f} sat={a['saturated']}")

    # ---- Null 3: shuffled coupling ----
    print("\nRunning null 3: shuffled coupling...")
    t0 = time.time()
    rng = np.random.default_rng(SEED + 7)
    K_shuf = shuffled_K(N, rng)
    cfg = RunConfig(N=N, T=T, transient=TRANSIENT, eps=EPS, tau=TAU,
                    drift_q=C_2I_DRIFT, K_mode="shuffled", K_matrix=K_shuf,
                    seed=SEED)
    res = run_dynamics(cfg)
    print(f"  ran in {time.time()-t0:.1f}s; mean_corr={res['meta']['mean_corr']:.4f} "
          f"std={res['meta']['std_corr']:.4f}")
    a = analyse_corr(res['corr'], "null:shuffled-K", N)
    rows.append(a)
    print(f"  d_s={a['d_s']:.3f} sat_frac={a['sat_frac']:.3f} sat={a['saturated']}")

    # Write CSV
    with open('observables_lag_v0_1.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print("\nWrote observables_lag_v0_1.csv")
    return rows


if __name__ == '__main__':
    rows = run_all()
    print("\n=== Summary ===")
    for r in rows:
        print(f"{r['label']:30s} d_s={r['d_s']:+.3f} sat={r['saturated']!s:5s} "
              f"slope={r['slope_volume']:+.3f} frac_pos_def={r['frac_pos_defect']:.3f}")
