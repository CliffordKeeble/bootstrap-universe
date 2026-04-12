# zeta_gen.py
# Generalised golden zeta builder — parameterised by alpha and norm
# Shared utility for compound investigation
# Bootstrap Universe Programme
# Written by Mr Code, April 12 2026

"""
Importable module. Generalises zeta_phi.py:
  - alpha: phase angle (phi, sqrt(2), e, pi, etc.)
  - norm_type: 'golden' (|Re^2 - 5*Im^2|) or 'circular' (Re^2 + Im^2)
"""

import math
import numpy as np

PHI = (1 + math.sqrt(5)) / 2

def compute_zeta_gen(t_array, N=5000, alpha=PHI, norm_type='golden'):
    """
    Compute generalised zeta function.

    Parameters:
        t_array: array of t values
        N: number of terms in Dirichlet sum
        alpha: irrational number for phase theta_n = 2*pi*{n*alpha}
        norm_type: 'golden' for sqrt(|Re^2 - 5*Im^2|)
                   'circular' for sqrt(Re^2 + Im^2)

    Returns: re, im, norm_values
    """
    t = np.asarray(t_array, dtype=np.float64)
    re = np.zeros_like(t)
    im = np.zeros_like(t)

    for n in range(1, N + 1):
        theta_n = 2 * math.pi * ((n * alpha) % 1.0)
        sigma_n = 1.0 - n / N  # Fejer smoothing
        amp = sigma_n / math.sqrt(n)

        phase = t * math.log(n) + theta_n
        re += amp * np.cos(phase)
        im += amp * np.sin(phase)

    if norm_type == 'golden':
        norm_values = np.sqrt(np.abs(re**2 - 5 * im**2))
    elif norm_type == 'circular':
        norm_values = np.sqrt(re**2 + im**2)
    else:
        raise ValueError(f"Unknown norm_type: {norm_type}")

    return re, im, norm_values


def run_test(riemann_zeros, t_array, N_terms=5000, alpha=PHI,
             norm_type='golden', n_mc=1000, rng_seed=42, label="test"):
    """
    Run a complete test: compute zeta, find minima, match, MC null, return z.

    Returns dict with all results.
    """
    from zeta_phi import find_minima, match_zeros

    t_max = t_array[-1] + (t_array[1] - t_array[0])

    # Compute
    re, im, norm_values = compute_zeta_gen(t_array, N=N_terms,
                                            alpha=alpha, norm_type=norm_type)

    # Find minima
    minima_t, minima_v = find_minima(t_array, norm_values,
                                      percentile=5, dedup_window=0.3)

    # Signal matching
    matches = match_zeros(minima_t, riemann_zeros, window=1.0)
    signal_deltas = np.array([m['delta'] for m in matches])
    signal_mean = np.mean(signal_deltas)

    # MC null
    rng = np.random.default_rng(rng_seed)
    null_means = []
    null_per_zero = np.zeros(len(riemann_zeros))

    for trial in range(n_mc):
        fake = np.sort(rng.uniform(t_array[0], t_max, size=len(minima_t)))
        fm = match_zeros(fake, riemann_zeros, window=1.0)
        fd = np.array([m['delta'] for m in fm])
        null_means.append(np.mean(fd))
        null_per_zero += fd

    null_per_zero /= n_mc
    null_arr = np.array(null_means)
    null_mean = np.mean(null_arr)
    null_se = np.std(null_arr)

    z = (signal_mean - null_mean) / null_se if null_se > 0 else 0
    p = np.mean(null_arr <= signal_mean)

    return {
        'label': label,
        'alpha': alpha,
        'norm_type': norm_type,
        'n_terms': N_terms,
        'n_zeros': len(riemann_zeros),
        'n_minima': len(minima_t),
        'density': len(minima_t) / (t_max - t_array[0]),
        'signal_mean_delta': signal_mean,
        'null_mean_delta': null_mean,
        'null_se': null_se,
        'z_overall': z,
        'p_value': p,
        'matches': matches,
        'signal_deltas': signal_deltas,
        'null_per_zero': null_per_zero,
        'minima_t': minima_t,
    }


def print_result(r):
    """Print a test result summary."""
    print(f"\n  [{r['label']}]")
    print(f"  alpha={r['alpha']:.6f}, norm={r['norm_type']}, N={r['n_terms']}")
    print(f"  Zeros: {r['n_zeros']}, Minima: {r['n_minima']}, "
          f"density: {r['density']:.3f}/unit")
    print(f"  Signal mean delta: {r['signal_mean_delta']:.6f}")
    print(f"  Null mean delta:   {r['null_mean_delta']:.6f} +/- {r['null_se']:.6f}")
    print(f"  z_overall:         {r['z_overall']:.4f}")
    print(f"  Empirical p:       {r['p_value']:.4f}")
