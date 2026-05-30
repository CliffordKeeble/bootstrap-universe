"""
Parameterised golden-angle Dirichlet probe with discriminant d.

zeta_phi_d(t) = sqrt(|Re^2(Z_phi(t)) - d * Im^2(Z_phi(t))|)

where Z_phi(t) = sum_{n=1}^N sigma_n / sqrt(n) * exp(i (t log n + theta_n))
with theta_n = 2*pi*{n*phi} and sigma_n = 1 - n/N (Fejer smoothing).

This is Paper 150 v2.0 zeta_gen.compute_zeta_gen with d parameterised.
"""

from __future__ import annotations

import math
import numpy as np

PHI = (1.0 + math.sqrt(5.0)) / 2.0


def compute_probe(t_array: np.ndarray, d: int, N: int = 5000,
                  alpha: float = PHI) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the parameterised probe.

    Returns (re, im, probe_value) where probe_value = sqrt(|re^2 - d*im^2|).
    """
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
    probe = np.sqrt(np.abs(re * re - float(d) * im * im))
    return re, im, probe


def find_minima(t_array: np.ndarray, values: np.ndarray,
                percentile: float = 5.0,
                dedup_window: float = 0.3) -> tuple[np.ndarray, np.ndarray]:
    """Three-point local minima with percentile threshold and dedup window."""
    threshold = np.percentile(values, percentile)
    cand_idx = []
    for i in range(1, len(values) - 1):
        if values[i] < values[i - 1] and values[i] < values[i + 1]:
            if values[i] < threshold:
                cand_idx.append(i)
    if not cand_idx:
        return np.array([]), np.array([])
    cand_t = t_array[cand_idx]
    cand_v = values[cand_idx]
    deduped_t = []
    deduped_v = []
    i = 0
    while i < len(cand_t):
        group_t = [cand_t[i]]
        group_v = [cand_v[i]]
        j = i + 1
        while j < len(cand_t) and cand_t[j] - cand_t[i] < dedup_window:
            group_t.append(cand_t[j])
            group_v.append(cand_v[j])
            j += 1
        best = int(np.argmin(group_v))
        deduped_t.append(group_t[best])
        deduped_v.append(group_v[best])
        i = j
    return np.array(deduped_t), np.array(deduped_v)


def match_to_targets(minima_t: np.ndarray, target_zeros: np.ndarray,
                      window: float = 1.0) -> list[dict]:
    """For each target zero, find nearest probe minimum."""
    matches = []
    for i, z in enumerate(target_zeros):
        if len(minima_t) == 0:
            matches.append(dict(n=i + 1, zero=z, min_t=float('nan'),
                                delta=float('nan'), matched=False))
            continue
        dists = np.abs(minima_t - z)
        best = int(np.argmin(dists))
        delta = dists[best]
        matches.append(dict(n=i + 1, zero=z, min_t=float(minima_t[best]),
                            delta=float(delta), matched=delta <= window))
    return matches


def mc_null(target_zeros: np.ndarray, n_minima: int, t_min: float,
            t_max: float, window: float = 1.0,
            n_trials: int = 1000, seed: int = 42) -> dict:
    """Monte Carlo null: draw fake minima uniformly, match to targets."""
    rng = np.random.default_rng(seed)
    null_means = []
    null_per_zero = np.zeros(len(target_zeros))
    for _ in range(n_trials):
        fake = np.sort(rng.uniform(t_min, t_max, size=n_minima))
        fm = match_to_targets(fake, target_zeros, window=window)
        fd = np.array([m['delta'] for m in fm])
        null_means.append(np.mean(fd))
        null_per_zero += fd
    null_per_zero /= n_trials
    arr = np.array(null_means)
    return dict(
        mean=float(np.mean(arr)),
        std=float(np.std(arr)),
        samples=arr,
        per_zero=null_per_zero,
    )


def run_cell(d_probe: int, target_zeros: np.ndarray, N_terms: int = 5000,
             t_pad: float = 5.0, dt: float = 0.008,
             window: float = 1.0, n_null: int = 1000,
             seed: int = 42) -> dict:
    """Run a single cell: probe-d against given target zeros, with MC null."""
    if len(target_zeros) == 0:
        return dict(error='no target zeros')
    t_min, t_max = max(1.0, float(target_zeros[0]) - t_pad), float(target_zeros[-1]) + t_pad
    t_array = np.arange(t_min, t_max, dt)
    re, im, probe = compute_probe(t_array, d=d_probe, N=N_terms)
    minima_t, minima_v = find_minima(t_array, probe, percentile=5.0, dedup_window=0.3)
    matches = match_to_targets(minima_t, target_zeros, window=window)
    signal_deltas = np.array([m['delta'] for m in matches])
    signal_mean = float(np.mean(signal_deltas))
    null = mc_null(target_zeros, n_minima=len(minima_t),
                   t_min=t_min, t_max=t_max, window=window,
                   n_trials=n_null, seed=seed)
    z = (signal_mean - null['mean']) / null['std'] if null['std'] > 0 else 0.0
    effect = (null['mean'] - signal_mean) / null['mean'] if null['mean'] > 0 else 0.0
    return dict(
        d_probe=d_probe,
        n_targets=len(target_zeros),
        n_minima=len(minima_t),
        density=len(minima_t) / (t_max - t_min),
        signal_mean=signal_mean,
        null_mean=null['mean'],
        null_std=null['std'],
        z=float(z),
        effect_size=float(effect),
        matched_count=int(sum(m['matched'] for m in matches)),
        t_min=t_min,
        t_max=t_max,
    )
