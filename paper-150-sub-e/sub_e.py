"""
Sub E: mechanism verification for Paper 150 v2.1.

Question: are minima locations of P_d(t) = |Re^2 - d*Im^2| shared
across d in {3, 5, 7, 13}?

Method:
1. Compute Z_phi (Re, Im) once on grid t in (1, 1000] at dt=0.01.
2. For each d, form P_d(t) = sqrt(|Re^2 - d * Im^2|).
3. Find 100 lowest local minima per d with non-overlap window 2.
4. For each pair (d, d'): count d's minima within 0.02 of any d' minimum.
5. Report 4x4 sharing-percentage matrix.

Pre-registration committed at cfd7599.
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path
import numpy as np

# Bring in the Paper 150 / Sub C-1 probe machinery (Z_phi components)
SUB_C_DIR = Path(__file__).parent.parent / "paper-203-sub-c"
sys.path.insert(0, str(SUB_C_DIR))
from probe import compute_probe  # uses alpha=phi, golden-norm parameterised by d


PHI = (1 + math.sqrt(5)) / 2

T_MIN = 1.0
T_MAX = 1000.0
DT = 0.01
N_TERMS = 5000
D_VALUES = [3, 5, 7, 13]
TOP_N = 100
NONOVERLAP_WINDOW = 2.0     # per brief: non-overlap window of width 2
SHARING_TOLERANCE = 2 * DT  # +/- 2*dt = 0.02 per brief


def compute_re_im_once(t_array: np.ndarray, N: int = N_TERMS,
                       alpha: float = PHI):
    """Compute Re(Z_phi), Im(Z_phi) ONCE for all d (they share this)."""
    t = np.asarray(t_array, dtype=np.float64)
    re = np.zeros_like(t)
    im = np.zeros_like(t)
    for n in range(1, N + 1):
        theta_n = 2.0 * math.pi * ((n * alpha) % 1.0)
        sigma_n = 1.0 - n / N
        amp = sigma_n / math.sqrt(n)
        phase = t * math.log(n) + theta_n
        re += amp * np.cos(phase)
        im += amp * np.sin(phase)
    return re, im


def probe_from_re_im(re: np.ndarray, im: np.ndarray, d: int) -> np.ndarray:
    return np.sqrt(np.abs(re * re - float(d) * im * im))


def find_local_minima(t_array: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Three-point local minima detection. Return (indices, t-values, probe-values)."""
    n = len(values)
    is_min = np.zeros(n, dtype=bool)
    is_min[1:-1] = (values[1:-1] < values[:-2]) & (values[1:-1] < values[2:])
    # Allow ties: <= on one side
    # Per brief: "P_d(t* - dt) >= P_d(t*) <= P_d(t* + dt)" — use strict < to avoid plateau noise.
    idx = np.where(is_min)[0]
    return idx, t_array[idx], values[idx]


def select_top_n_with_nonoverlap(t_min: np.ndarray, v_min: np.ndarray,
                                  top_n: int = 100,
                                  window: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Greedy selection: pick lowest, exclude all within +/- window/2 of it
    (so non-overlap windows of width `window`), pick next lowest, etc.
    Returns (selected t-values, selected probe-values) of length top_n
    (or fewer if not enough minima).
    """
    # Sort by value
    order = np.argsort(v_min)
    t_sorted = t_min[order]
    v_sorted = v_min[order]
    selected_t = []
    selected_v = []
    # Use a simple in-list check; OK since top_n <= 100 and total minima ~10k
    for ti, vi in zip(t_sorted, v_sorted):
        too_close = False
        for tj in selected_t:
            if abs(ti - tj) < window:
                too_close = True
                break
        if not too_close:
            selected_t.append(ti)
            selected_v.append(vi)
            if len(selected_t) >= top_n:
                break
    return np.array(selected_t), np.array(selected_v)


def sharing_pct(t_d: np.ndarray, t_dprime: np.ndarray, tol: float) -> float:
    """For each t in t_d, is there a t' in t_dprime with |t - t'| <= tol?"""
    if len(t_d) == 0 or len(t_dprime) == 0:
        return 0.0
    t_dprime_sorted = np.sort(t_dprime)
    matched = 0
    for t in t_d:
        idx = np.searchsorted(t_dprime_sorted, t)
        candidates = []
        if idx < len(t_dprime_sorted):
            candidates.append(t_dprime_sorted[idx])
        if idx > 0:
            candidates.append(t_dprime_sorted[idx - 1])
        if any(abs(t - c) <= tol for c in candidates):
            matched += 1
    return 100.0 * matched / len(t_d)


def main():
    print("=" * 70)
    print("Sub E: mechanism verification (minima sharing)")
    print(f"Grid: t in ({T_MIN}, {T_MAX}] at dt={DT}; N_terms={N_TERMS}")
    print(f"Per d, top {TOP_N} lowest minima with non-overlap window {NONOVERLAP_WINDOW}")
    print(f"Sharing tolerance: +/-{SHARING_TOLERANCE} t-units")
    print("=" * 70)

    t_array = np.arange(T_MIN + DT, T_MAX + DT/2, DT)
    print(f"\nGrid has {len(t_array)} points")

    print(f"\nComputing Z_phi (Re, Im) once for shared use...")
    t0 = time.time()
    re, im = compute_re_im_once(t_array, N=N_TERMS, alpha=PHI)
    print(f"  Done in {time.time()-t0:.1f}s")

    # For each d, compute probe, find minima, select top 100
    minima_per_d = {}
    for d in D_VALUES:
        t0 = time.time()
        probe = probe_from_re_im(re, im, d)
        idx, t_min, v_min = find_local_minima(t_array, probe)
        sel_t, sel_v = select_top_n_with_nonoverlap(
            t_min, v_min, top_n=TOP_N, window=NONOVERLAP_WINDOW)
        print(f"  d={d}: total local minima = {len(idx)}; "
              f"top-{len(sel_t)} selected (lowest v = {sel_v[0]:.6f}, "
              f"highest selected = {sel_v[-1]:.6f}) in {time.time()-t0:.1f}s")
        minima_per_d[d] = sel_t

    # Compute 4x4 sharing matrix
    print(f"\nSharing matrix (% of row's 100 lowest minima within {SHARING_TOLERANCE} of any column minimum):")
    matrix = np.zeros((len(D_VALUES), len(D_VALUES)))
    for i, d in enumerate(D_VALUES):
        for j, dp in enumerate(D_VALUES):
            matrix[i, j] = sharing_pct(minima_per_d[d], minima_per_d[dp],
                                        SHARING_TOLERANCE)

    # Print matrix
    print(f"\n  {'':>8s} | " + " | ".join(f"d'={dp:2d}" for dp in D_VALUES))
    print("-" * 50)
    for i, d in enumerate(D_VALUES):
        row = " | ".join(f"{matrix[i, j]:5.1f}%" for j in range(len(D_VALUES)))
        print(f"  d={d:2d}    | {row}")

    # Mean off-diagonal sharing
    n = len(D_VALUES)
    off_diag = []
    for i in range(n):
        for j in range(n):
            if i != j:
                off_diag.append(matrix[i, j])
    mean_off = np.mean(off_diag)
    print(f"\nMean off-diagonal sharing: {mean_off:.2f}%")

    if mean_off >= 80:
        verdict = "DECORATION-CONFIRMED"
    elif mean_off >= 40:
        verdict = "DECORATION-PARTIAL"
    else:
        verdict = "DECORATION-REFUTED"
    print(f"\nVerdict: {verdict}")

    # Save results
    with open('sharing_matrix.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['d_row'] + [f'd_col_{dp}' for dp in D_VALUES])
        for i, d in enumerate(D_VALUES):
            w.writerow([d] + [f"{matrix[i, j]:.2f}" for j in range(len(D_VALUES))])
    print("\nWrote sharing_matrix.csv")

    # Save the selected minima per d for traceability
    with open('minima_per_d.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['d', 'rank', 't_value'])
        for d in D_VALUES:
            for rank, t in enumerate(minima_per_d[d], start=1):
                w.writerow([d, rank, f"{t:.6f}"])
    print("Wrote minima_per_d.csv")

    return dict(matrix=matrix, mean_off=mean_off, verdict=verdict,
                minima=minima_per_d)


if __name__ == '__main__':
    main()
