"""
Sub B: zeta-value-denominator chance-overlap null.

Population Z (per pre-registration): denom(zeta(-2k+1)) for k = 1..50.
Population I (per pre-registration, exactly per brief):
    {1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, 120, 168, 240}.

Poisson null: mu = |Z| * |I| / M with M = max(Z u I). z = (actual - mu)/sqrt(mu).
Tight null: replace |Z|/M with the empirical density of Z values in {1..M}.
"""

from __future__ import annotations

import math
from sympy import bernoulli, Rational, S

# Population I, frozen at pre-registration
I_POPULATION = sorted({1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, 120, 168, 240})


def zeta_value_denominator(n: int) -> int:
    """
    denom(zeta(-n)) for n >= 1.
    zeta(-n) = -B_{n+1} / (n+1).
    Reduce the fraction and return the denominator.
    """
    B = bernoulli(n + 1)             # exact rational
    val = Rational(-B, n + 1)        # reduces automatically
    return int(val.q)


def build_z_population(N: int = 50) -> list[int]:
    """Z_k = denom(zeta(-2k+1)) for k = 1..N. Negative odd integers => Bernoulli at even index."""
    z = []
    for k in range(1, N + 1):
        n = 2 * k - 1
        denom = zeta_value_denominator(n)
        z.append(denom)
    return z


def run_sub_b(N: int = 50, verbose: bool = True) -> dict:
    Z = build_z_population(N)
    I = I_POPULATION

    Z_set = set(Z)
    I_set = set(I)
    matches = sorted(Z_set & I_set)
    n_matches = len(matches)

    M = max(max(Z), max(I))

    # Coarse Poisson null
    mu_coarse = len(Z) * len(I) / M
    z_coarse = (n_matches - mu_coarse) / math.sqrt(mu_coarse) if mu_coarse > 0 else float('inf')

    # Tight null: density of Z in {1..M}
    unique_Z_in_M = len(Z_set)
    density_Z = unique_Z_in_M / M
    expected_under_tight = unique_Z_in_M * len(I_set) / M
    # equivalent computation: same as coarse if |Z|=unique_Z_in_M. So tight = coarse minus duplicates.
    # Actually a more meaningful tight null: P(a uniformly-chosen element of {1..M} hits I) = |I|/M.
    # Probability that any element of Z hits I = 1 - (1 - |I|/M)^|unique_Z_in_M|.
    p_single = len(I_set) / M
    expected_tight = unique_Z_in_M * p_single
    z_tight = (n_matches - expected_tight) / math.sqrt(expected_tight) if expected_tight > 0 else float('inf')

    if verbose:
        print(f"Population Z ({len(Z)} elements, {len(Z_set)} unique):")
        print(f"  First 10 (k=1..10, n=1..19): {Z[:10]}")
        print(f"  Sample (k=1, 2, 3, 4, 5, 10, 20, 30, 40, 50):")
        for k in [1, 2, 3, 4, 5, 10, 20, 30, 40, 50]:
            print(f"    k={k:3d}, n=2k-1={2*k-1:3d}: denom(zeta(-{2*k-1})) = {Z[k-1]}")
        print(f"\nPopulation I ({len(I)} elements): {I}")
        print(f"\nMatches: {n_matches} out of {len(Z)} ({n_matches/len(Z)*100:.1f}%)")
        print(f"  Matched values: {matches}")
        print(f"\nM = max(Z u I) = {M}")
        print(f"\nCoarse Poisson null:")
        print(f"  mu = |Z| * |I| / M = {len(Z)} * {len(I)} / {M} = {mu_coarse:.4f}")
        print(f"  actual = {n_matches}, expected = {mu_coarse:.4f}")
        print(f"  z = ({n_matches} - {mu_coarse:.4f}) / sqrt({mu_coarse:.4f}) = {z_coarse:.4f}")
        print(f"\nTight null (unique-Z density):")
        print(f"  expected = |unique Z| * |I| / M = {unique_Z_in_M} * {len(I)} / {M} = {expected_tight:.4f}")
        print(f"  z = ({n_matches} - {expected_tight:.4f}) / sqrt({expected_tight:.4f}) = {z_tight:.4f}")

        # Pre-registered verdict
        max_z = max(abs(z_coarse), abs(z_tight))
        if max_z >= 3:
            verdict = "STRIKING"
        elif max_z >= 2:
            verdict = "CONSISTENT WITH FRAMEWORK ONLY"
        else:
            verdict = "CONSISTENT WITH CHANCE"
        print(f"\nVerdict (using max |z| = {max_z:.2f}): {verdict}")

    # SUPPLEMENTARY null: restrict to Z values <= max(I).
    # This is the more meaningful chance-overlap null: only Z values in
    # the size range of I can possibly match I, so condition on that.
    M_I = max(I_set)
    Z_small = sorted({z for z in Z_set if z <= M_I})
    n_small_Z = len(Z_small)
    matches_small = sorted(set(Z_small) & I_set)
    n_matches_small = len(matches_small)
    # Probability a "small" Z (uniform in {1..M_I}) lands in I:
    p_in_I = len(I_set) / M_I
    expected_small = n_small_Z * p_in_I
    z_small = (n_matches_small - expected_small) / math.sqrt(expected_small) if expected_small > 0 else float('inf')

    if verbose:
        print(f"\n--- Supplementary tight null (Z restricted to <= max(I) = {M_I}) ---")
        print(f"  unique Z values <= {M_I}: {n_small_Z} ({Z_small})")
        print(f"  matches with I: {n_matches_small} ({matches_small})")
        print(f"  expected = {n_small_Z} * {len(I_set)}/{M_I} = {expected_small:.4f}")
        print(f"  z_small = {z_small:.4f}")
        if abs(z_small) >= 3:
            v2 = "STRIKING"
        elif abs(z_small) >= 2:
            v2 = "CONSISTENT WITH FRAMEWORK ONLY"
        else:
            v2 = "CONSISTENT WITH CHANCE"
        print(f"  supplementary verdict (|z_small|={abs(z_small):.2f}): {v2}")

    return dict(
        Z=Z,
        Z_unique=sorted(Z_set),
        I=I,
        matches=matches,
        n_matches=n_matches,
        M=M,
        mu_coarse=mu_coarse,
        z_coarse=z_coarse,
        expected_tight=expected_tight,
        z_tight=z_tight,
        # supplementary
        M_I=M_I,
        Z_small=Z_small,
        n_small_Z=n_small_Z,
        matches_small=matches_small,
        n_matches_small=n_matches_small,
        expected_small=expected_small,
        z_small=z_small,
    )


if __name__ == '__main__':
    print("=" * 70)
    print("Sub B: zeta-value-denominator chance-overlap null")
    print("=" * 70)
    res = run_sub_b(N=50, verbose=True)
    # Write CSV
    import csv
    with open('z_population.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['k', 'n', 'denom_zeta_minus_n'])
        for k, denom in enumerate(res['Z'], start=1):
            w.writerow([k, 2*k - 1, denom])
    print("\nWrote z_population.csv")
