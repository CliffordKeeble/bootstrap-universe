# zeta_phi.py
# Golden zeta function ζ_φ(t): construction, minima detection, Riemann zero matching
# Reproduces Paper 150 §2 baseline
# Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 12 2026

"""
Standalone: python zeta_phi.py
Outputs baseline results to stdout.
Also importable: from zeta_phi import compute_zeta_phi, find_minima, match_zeros
"""

import math
import numpy as np
from mpmath import zetazero

# ── Constants ─────────────────────────────────────────────────────────────────

PHI = (1 + math.sqrt(5)) / 2  # golden ratio

# ── ζ_φ construction (Paper 150 §2) ──────────────────────────────────────────

def compute_zeta_phi(t_array, N=5000):
    """
    Compute ζ_φ(t) for an array of t values.

    Construction:
      ζ_φ(t) = Σ_{n=1}^{N} σ_n · cos(t·log(n) + θ_n) / n^{1/2}   (real part)
      ζ_φ_im(t) = Σ_{n=1}^{N} σ_n · sin(t·log(n) + θ_n) / n^{1/2}  (imag part)

    where:
      θ_n = 2π·{n·φ}      (golden-angle phase, {x} = fractional part)
      σ_n = 1 - n/N        (Fejér smoothing)

    Golden norm: ||ζ_φ|| = sqrt(|Re² - 5·Im²|)

    Returns: t_array, golden_norm array
    """
    t = np.asarray(t_array, dtype=np.float64)
    re = np.zeros_like(t)
    im = np.zeros_like(t)

    for n in range(1, N + 1):
        theta_n = 2 * math.pi * ((n * PHI) % 1.0)  # golden-angle phase
        sigma_n = 1.0 - n / N                        # Fejér smoothing
        amp = sigma_n / math.sqrt(n)

        log_n = math.log(n)
        phase = t * log_n + theta_n

        re += amp * np.cos(phase)
        im += amp * np.sin(phase)

    # Golden norm: sqrt(|Re² - 5·Im²|)
    golden_norm = np.sqrt(np.abs(re**2 - 5 * im**2))

    return re, im, golden_norm


def find_minima(t_array, values, percentile=5, dedup_window=0.3):
    """
    Find local minima below the given percentile threshold.
    Three-point local comparison, deduplication within window.

    Returns: array of t values at minima, array of values at minima
    """
    threshold = np.percentile(values, percentile)

    # Three-point local minima
    candidates_idx = []
    for i in range(1, len(values) - 1):
        if values[i] < values[i-1] and values[i] < values[i+1]:
            if values[i] < threshold:
                candidates_idx.append(i)

    if not candidates_idx:
        return np.array([]), np.array([])

    # Deduplication: keep deepest minimum within each window
    candidates_t = t_array[candidates_idx]
    candidates_v = values[candidates_idx]

    deduped_t = []
    deduped_v = []

    i = 0
    while i < len(candidates_t):
        # Collect all candidates within dedup_window of current
        group_t = [candidates_t[i]]
        group_v = [candidates_v[i]]
        j = i + 1
        while j < len(candidates_t) and candidates_t[j] - candidates_t[i] < dedup_window:
            group_t.append(candidates_t[j])
            group_v.append(candidates_v[j])
            j += 1

        # Keep the deepest
        best = np.argmin(group_v)
        deduped_t.append(group_t[best])
        deduped_v.append(group_v[best])
        i = j

    return np.array(deduped_t), np.array(deduped_v)


def get_riemann_zeros(count=100):
    """Get first `count` nontrivial Riemann zeta zeros (imaginary parts)."""
    zeros = []
    for n in range(1, count + 1):
        z = float(zetazero(n).imag)
        zeros.append(z)
    return np.array(zeros)


def match_zeros(minima_t, riemann_zeros, window=1.0):
    """
    Match each Riemann zero to its nearest ζ_φ minimum.
    Returns list of dicts with zero, matched minimum, delta, and whether within window.
    """
    matches = []
    for i, z in enumerate(riemann_zeros):
        if len(minima_t) == 0:
            matches.append({
                'n': i + 1, 'zero': z, 'min_t': np.nan,
                'delta': np.nan, 'matched': False
            })
            continue

        dists = np.abs(minima_t - z)
        best_idx = np.argmin(dists)
        delta = dists[best_idx]

        matches.append({
            'n': i + 1,
            'zero': z,
            'min_t': minima_t[best_idx],
            'delta': delta,
            'matched': delta <= window
        })

    return matches


def bin_matches(matches, bins=None):
    """Bin match statistics by t-range."""
    if bins is None:
        bins = [(0, 80), (80, 160), (160, 237)]

    results = {}
    for lo, hi in bins:
        in_bin = [m for m in matches if lo <= m['zero'] < hi]
        if not in_bin:
            continue
        deltas = [m['delta'] for m in in_bin if not np.isnan(m['delta'])]
        matched = [m for m in in_bin if m['matched']]
        results[(lo, hi)] = {
            'count': len(in_bin),
            'matched': len(matched),
            'match_rate': len(matched) / len(in_bin) if in_bin else 0,
            'mean_delta': np.mean(deltas) if deltas else np.nan,
            'median_delta': np.median(deltas) if deltas else np.nan,
            'max_delta': np.max(deltas) if deltas else np.nan,
        }
    return results


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 70)
    print("GOLDEN ZETA BASELINE REPRODUCTION (Paper 150 §2)")
    print("=" * 70)

    # Parameters
    N = 5000
    t_min, t_max, dt = 1.0, 237.0, 0.008
    t_array = np.arange(t_min, t_max, dt)

    print(f"\nParameters: N={N}, t=[{t_min}, {t_max}], dt={dt}")
    print(f"Grid points: {len(t_array)}")

    # Step 1: Compute ζ_φ
    print("\nComputing zeta_phi...")
    re, im, golden_norm = compute_zeta_phi(t_array, N=N)

    print(f"  Golden norm: mean={golden_norm.mean():.4f}, "
          f"std={golden_norm.std():.4f}, "
          f"min={golden_norm.min():.6f}, max={golden_norm.max():.4f}")

    # Step 2: Find minima
    print("\nFinding minima (5th percentile, dedup window 0.3)...")
    minima_t, minima_v = find_minima(t_array, golden_norm,
                                      percentile=5, dedup_window=0.3)
    print(f"  Minima found: {len(minima_t)}")
    print(f"  Density: {len(minima_t) / (t_max - t_min):.3f} per unit t")

    # Step 3: Get Riemann zeros
    print("\nFetching first 100 Riemann zeros...")
    riemann_zeros = get_riemann_zeros(100)
    print(f"  Range: [{riemann_zeros[0]:.4f}, {riemann_zeros[-1]:.4f}]")
    # How many fall within our t range?
    in_range = riemann_zeros[riemann_zeros <= t_max]
    print(f"  In range [1, {t_max}]: {len(in_range)}")

    # Step 4: Match at W=1.0 (Paper 150 default)
    print("\nMatching zeros to minima (window W=1.0)...")
    matches = match_zeros(minima_t, in_range, window=1.0)

    total_matched = sum(1 for m in matches if m['matched'])
    all_deltas = [m['delta'] for m in matches if not np.isnan(m['delta'])]

    print(f"  Total matched: {total_matched}/{len(matches)}")
    print(f"  Overall mean delta: {np.mean(all_deltas):.4f}")
    print(f"  Overall median delta: {np.median(all_deltas):.4f}")
    print(f"  Overall max delta: {np.max(all_deltas):.4f}")

    # Step 5: Binned results
    print("\n-- Binned match statistics (W=1.0) --")
    binned = bin_matches(matches)
    print(f"  {'Range':>12s}  {'N':>4s}  {'Match':>5s}  {'Rate':>6s}  "
          f"{'Mean_D':>7s}  {'Med_D':>7s}  {'Max_D':>7s}")
    for (lo, hi), stats in sorted(binned.items()):
        print(f"  [{lo:3.0f},{hi:3.0f})  {stats['count']:4d}  "
              f"{stats['matched']:5d}  {stats['match_rate']:6.1%}  "
              f"{stats['mean_delta']:7.4f}  {stats['median_delta']:7.4f}  "
              f"{stats['max_delta']:7.4f}")

    # Step 6: First 20 vs rest
    print("\n-- First 20 zeros vs zeros 21-100 --")
    first20 = [m for m in matches if m['n'] <= 20]
    rest = [m for m in matches if m['n'] > 20]

    d20 = [m['delta'] for m in first20 if not np.isnan(m['delta'])]
    d_rest = [m['delta'] for m in rest if not np.isnan(m['delta'])]

    if d20:
        print(f"  First 20:  mean_delta={np.mean(d20):.4f}, "
              f"median={np.median(d20):.4f}, max={np.max(d20):.4f}")
    if d_rest:
        print(f"  Zeros 21+: mean_delta={np.mean(d_rest):.4f}, "
              f"median={np.median(d_rest):.4f}, max={np.max(d_rest):.4f}")

    print("\n" + "=" * 70)
    print("BASELINE COMPLETE")
    print("=" * 70)
