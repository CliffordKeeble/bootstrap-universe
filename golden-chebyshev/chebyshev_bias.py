# chebyshev_bias.py - Layer 1: direct count of E(x) for q=5
# Bootstrap Universe Programme - Paper 196 candidate
# Companion to addendum 1 (Granville-Martin extract) of Mr Code brief
#
# Standalone: python chebyshev_bias.py
#
# Computes E(x) = #{p <= x : p == +/-2 mod 5}  -  #{p <= x : p == +/-1 mod 5}
# at log-spaced x up to X_MAX. Stubborns minus splitters; positive E(x) = stubborn lead.
#
# Outputs:
#   chebyshev_bias_data.csv  - x, E(x), normalised E(x)/(sqrt(x)/ln(x)), counts per class
#   chebyshev_bias.png       - plot of E(x)/(sqrt(x)/ln(x)) vs log10(x)
#
# Sanity check: at x = 10^6, the count by residue class is reported and compared
# to Granville-Martin 2006 Table 4 (mod-10 reorganisation; see addendum).

import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np

# ---------- Sieve ----------

def sieve_of_eratosthenes(n):
    """Return numpy array of primes <= n."""
    if n < 2:
        return np.array([], dtype=np.int64)
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if is_prime[i]:
            is_prime[i * i::i] = False
    return np.flatnonzero(is_prime).astype(np.int64)


def segmented_sieve(n, segment=10_000_000):
    """
    Memory-friendlier sieve for larger n. Returns numpy array of primes <= n.
    Uses base-prime sieve up to sqrt(n), then sieves segments by trial division
    with the base primes. Returns primes in ascending order.
    """
    if n < segment:
        return sieve_of_eratosthenes(n)
    base = sieve_of_eratosthenes(int(math.isqrt(n)) + 1)
    out = [base.copy()]
    lo = base[-1] + 1 if base.size else 2
    while lo <= n:
        hi = min(lo + segment - 1, n)
        size = hi - lo + 1
        mark = np.ones(size, dtype=bool)
        for p in base:
            if p * p > hi:
                break
            start = ((lo + p - 1) // p) * p
            if start < p * p:
                start = p * p
            if start > hi:
                continue
            mark[start - lo::p] = False
        seg_primes = np.flatnonzero(mark) + lo
        out.append(seg_primes.astype(np.int64))
        lo = hi + 1
    return np.concatenate(out)


# ---------- Classification ----------

# Residues mod 5: splitters {1, 4} (squares), stubborns {2, 3} (non-squares).
# p = 5 itself is ramified - excluded from both.

SPLITTER_RES = (1, 4)
STUBBORN_RES = (2, 3)


def classify(primes):
    """
    Given a sorted numpy array of primes, return arrays:
      is_splitter, is_stubborn  - boolean masks (mutually exclusive; p=5 in neither)
    """
    r = primes % 5
    is_splitter = (r == 1) | (r == 4)
    is_stubborn = (r == 2) | (r == 3)
    # p=5 (r==0) is in neither.
    return is_splitter, is_stubborn


def cumulative_counts(primes, mask):
    """Cumulative count of primes satisfying mask, indexed by prime position."""
    return np.cumsum(mask.astype(np.int64))


# ---------- E(x) at log-spaced sample points ----------

def evaluate_E(primes, is_split, is_stub, x_samples):
    """
    For each x in x_samples (sorted ascending), return:
      n_split[x], n_stub[x], E(x) = n_stub[x] - n_split[x]
    Vectorised via searchsorted.
    """
    cum_split = cumulative_counts(primes, is_split)
    cum_stub = cumulative_counts(primes, is_stub)
    # Index of last prime <= x.
    idx = np.searchsorted(primes, x_samples, side='right') - 1
    valid = idx >= 0
    n_split = np.where(valid, cum_split[idx.clip(min=0)], 0)
    n_stub = np.where(valid, cum_stub[idx.clip(min=0)], 0)
    return n_split, n_stub, n_stub - n_split


def normalised_E(x_samples, E):
    """E(x) / (sqrt(x) / ln(x)) - the Granville-Martin normalisation."""
    x = np.asarray(x_samples, dtype=np.float64)
    denom = np.sqrt(x) / np.log(x)
    return E / denom


# ---------- Sanity check vs Granville-Martin Table 4 (x = 10^6) ----------

def sanity_count_at_million(primes, is_split, is_stub):
    """
    Print per-residue-class prime counts up to 10^6 for sanity-checking against
    Granville-Martin 2006 Table 4 (their mod-10 ordering reorganises to mod-5).

    Mod-5 residues: 1, 2, 3, 4 (and 5 itself).
    Last digit (mod-10) reorganisation:
      last digit 1, 11, 21, ... -> mod-5 = 1
      last digit 9 -> mod-5 = 4
      last digit 3 -> mod-5 = 3
      last digit 7 -> mod-5 = 2
    For primes > 5 we have last digit in {1, 3, 7, 9}, so the mod-5 classes
    correspond cleanly: splitters = {last 1, last 9}; stubborns = {last 3, last 7}.
    """
    X = 1_000_000
    idx = np.searchsorted(primes, X, side='right')
    p_under = primes[:idx]
    r = p_under % 5
    counts = {res: int(np.sum(r == res)) for res in (0, 1, 2, 3, 4)}
    n_split = counts[1] + counts[4]
    n_stub = counts[2] + counts[3]
    E_at_M = n_stub - n_split
    print(f"  pi(10^6, 5, r) by residue:")
    for res in (1, 2, 3, 4):
        kind = "splitter" if res in SPLITTER_RES else "stubborn"
        print(f"    r = {res} ({kind:>8s}): {counts[res]:6d}")
    print(f"    r = 0 (p = 5):          {counts[0]:6d}")
    print(f"  splitters total: {n_split:6d}")
    print(f"  stubborns total: {n_stub:6d}")
    print(f"  E(10^6) = stub - split = {E_at_M:+d}")
    # Granville-Martin Table 4 anchor (paraphrased — direct count is source of truth).
    return E_at_M


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="Layer 1: Chebyshev bias direct count for q=5")
    parser.add_argument('--xmax', type=int, default=10**8,
                        help='Sieve up to this bound (default 10^8; 10^9 is ~30s and ~6 GB RAM)')
    parser.add_argument('--segment', type=int, default=20_000_000,
                        help='Segmented-sieve segment size in bytes-of-bool')
    parser.add_argument('--samples', type=int, default=200,
                        help='Number of log-spaced sample points')
    parser.add_argument('--outdir', type=str, default=str(Path(__file__).parent),
                        help='Output directory')
    args = parser.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("LAYER 1: Chebyshev bias direct count for q = 5")
    print(f"  X_MAX = {args.xmax:,}")
    print("=" * 70)

    print("\nSieving primes...")
    t0 = time.time()
    primes = segmented_sieve(args.xmax, segment=args.segment)
    print(f"  {len(primes):,} primes <= {args.xmax:,} in {time.time() - t0:.1f}s")

    is_split, is_stub = classify(primes)
    n_split_total = int(is_split.sum())
    n_stub_total = int(is_stub.sum())
    print(f"  splitters: {n_split_total:,}")
    print(f"  stubborns: {n_stub_total:,}")
    print(f"  E(X_MAX) = {n_stub_total - n_split_total:+d}")

    print("\nSanity vs Granville-Martin Table 4 anchor (x = 10^6):")
    if args.xmax >= 1_000_000:
        sanity_count_at_million(primes, is_split, is_stub)
    else:
        print("  X_MAX < 10^6, skipping sanity check")

    print("\nEvaluating E(x) at log-spaced sample points...")
    x_samples = np.unique(np.round(
        np.geomspace(100, args.xmax, args.samples)
    ).astype(np.int64))
    n_split, n_stub, E = evaluate_E(primes, is_split, is_stub, x_samples)
    E_norm = normalised_E(x_samples, E)

    csv_path = out / 'chebyshev_bias_data.csv'
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['x', 'n_splitters', 'n_stubborns', 'E', 'E_normalised'])
        for xs, ns, nb, e, en in zip(x_samples, n_split, n_stub, E, E_norm):
            w.writerow([int(xs), int(ns), int(nb), int(e), f"{en:.6f}"])
    print(f"  wrote {csv_path}")

    # Sign-frequency over log(x).
    pos_frac = float(np.mean(E > 0))
    neg_frac = float(np.mean(E < 0))
    zero_frac = float(np.mean(E == 0))
    print(f"\n  Sign frequency over {len(x_samples)} log-spaced samples:")
    print(f"    E(x) > 0 (stubborn lead): {pos_frac:.4f}")
    print(f"    E(x) < 0 (splitter lead): {neg_frac:.4f}")
    print(f"    E(x) = 0:                 {zero_frac:.4f}")
    print("  (Empirical density - to compare to R-S logarithmic density in density_qch5.py)")

    # Plot.
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4.5))
        ax.axhline(0, color='gray', linewidth=0.5)
        ax.semilogx(x_samples, E_norm, 'b-', linewidth=0.8, label=r'$E(x) / (\sqrt{x}/\ln x)$')
        ax.set_xlabel('x')
        ax.set_ylabel(r'$E(x) / (\sqrt{x}/\ln x)$')
        ax.set_title(r'Chebyshev bias for $q=5$: stubborns ($\pm 2$) minus splitters ($\pm 1$)')
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(loc='best')
        plot_path = out / 'chebyshev_bias.png'
        fig.tight_layout()
        fig.savefig(plot_path, dpi=120)
        plt.close(fig)
        print(f"  wrote {plot_path}")
    except ImportError:
        print("  (matplotlib not installed - skipping plot)")

    print("\n" + "=" * 70)
    print("LAYER 1 COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
