#!/usr/bin/env python3
"""
Knuth Circle Product — Implementation and Analysis
====================================================

Implements Fibonacci multiplication (the "circle product" ⊚) as defined by
Knuth (1988), operating on Zeckendorf representations via base-φ polynomial
multiplication with full normalisation.

Key implementation detail:
    The carry rule 2·φ^p → φ^{p+1} + φ^{p-2} can push digits to negative
    base-φ positions. These positions have real φ-values (φ^{-1} = φ-1 ≈ 0.618,
    φ^{-2} = 2-φ ≈ 0.382, etc.) and MUST be tracked through intermediate
    computations. Dropping or redirecting them (e.g. treating F_0=0 as "discard",
    or F_1=F_2 as "merge to index 2") destroys associativity.

    This was verified through nine iterations of debugging. The product is
    provably associative (Knuth 1988); any implementation that fails the
    associativity test has a normalisation bug.

Reference:
    D. E. Knuth, "Fibonacci Multiplication", Applied Mathematics Letters,
    Vol. 1, No. 1, pp. 57-60, 1988.

Part of the Bootstrap Universe Programme — Paper 168 computational appendix.
Dr. Clifford Keeble, Woodbridge, UK. March 2026.
ORCID: 0009-0003-6828-2155
"""

import math
from collections import Counter

# =====================================================================
# FIBONACCI NUMBERS (extended to negative indices)
# =====================================================================

MAX_FIB = 100

_fib_cache = {0: 0, 1: 1, 2: 1}

def fib(n):
    """Extended Fibonacci: F_{-k} = (-1)^{k+1} F_k."""
    if n in _fib_cache:
        return _fib_cache[n]
    if n < 0:
        val = ((-1) ** ((-n) + 1)) * fib(-n)
    else:
        val = fib(n - 1) + fib(n - 2)
    _fib_cache[n] = val
    return val

# Precompute
for _i in range(-30, MAX_FIB + 1):
    fib(_i)


# =====================================================================
# CORE OPERATIONS
# =====================================================================

def int_to_positions(n):
    """Convert positive integer n to a set of base-φ positions.

    The Zeckendorf representation n = Σ F_{k_i} (greedy, non-adjacent, k≥2)
    maps to base-φ positions p_i = k_i - 2, so F_{k} → φ^{k-2}.

    Returns a frozenset of positions (integers, all ≥ 0 for valid input).
    """
    if n <= 0:
        return frozenset()
    r, pos = n, []
    for p in range(MAX_FIB - 2, -1, -1):
        f = fib(p + 2)  # F_{p+2} at position p
        if f <= r and f > 0:
            pos.append(p)
            r -= f
            if r == 0:
                break
    assert r == 0, f"Zeckendorf decomposition failed for {n}"
    return frozenset(pos)


def positions_to_int(pos):
    """Convert a position set back to an integer via Σ F_{p+2}."""
    return sum(fib(p + 2) for p in pos)


def normalise(counts):
    """Normalise a count dictionary of base-φ positions to valid Zeckendorf form.

    Input: {position: count} where count may be > 1 and position may be negative.
    Output: frozenset of positions (each appearing once, no two consecutive).

    Rules applied repeatedly until stable:
        (A) Duplicate:    2·φ^p → φ^{p+1} + φ^{p-2}  (ALL p, including negative)
        (B) Consecutive:  φ^p + φ^{p+1} → φ^{p+2}

    Duplicates are processed highest-first (pushes carries upward, preventing
    infinite loops). Consecutives are processed lowest-first.
    """
    c = dict(counts)
    c = {k: v for k, v in c.items() if v > 0}

    for _ in range(100000):
        changed = False

        # Phase 1: Eliminate any count ≥ 2 (highest position first)
        for k in sorted(c.keys(), reverse=True):
            if c.get(k, 0) >= 2:
                c[k] -= 2
                if c[k] == 0:
                    del c[k]
                c[k + 1] = c.get(k + 1, 0) + 1
                c[k - 2] = c.get(k - 2, 0) + 1  # ALWAYS, even for negative k-2
                changed = True
                break

        if changed:
            continue

        # Phase 2: Eliminate consecutive pairs (lowest position first)
        keys = sorted(c.keys())
        for i in range(len(keys) - 1):
            k = keys[i]
            if keys[i + 1] == k + 1 and c.get(k, 0) > 0 and c.get(k + 1, 0) > 0:
                c[k] -= 1
                c[k + 1] -= 1
                if c[k] == 0:
                    del c[k]
                if c.get(k + 1, 0) == 0:
                    del c[k + 1]
                c[k + 2] = c.get(k + 2, 0) + 1
                changed = True
                break

        if not changed:
            break

    return frozenset(k for k, v in c.items() if v > 0)


def rep_multiply(pos_a, pos_b):
    """Multiply two position sets as polynomials in φ.

    Product rule: position i × position j → position i + j
    (since φ^i · φ^j = φ^{i+j}).

    Returns normalised position set.
    """
    counts = {}
    for i in pos_a:
        for j in pos_b:
            p = i + j
            counts[p] = counts.get(p, 0) + 1
    return normalise(counts)


def knuth_product(m, n):
    """Compute the Knuth circle product m ⊚ n.

    The circle product interprets integers as polynomials in φ via their
    Zeckendorf representations, multiplies the polynomials, and converts
    back to an integer.

    Properties (verified computationally, proven by Knuth):
        - Identity: 1 ⊚ n = n
        - Commutativity: m ⊚ n = n ⊚ m
        - Associativity: (m ⊚ n) ⊚ k = m ⊚ (n ⊚ k)
        - NOT distributive over addition
    """
    if m == 0 or n == 0:
        return 0
    if m == 1:
        return n
    if n == 1:
        return m
    return positions_to_int(rep_multiply(int_to_positions(m), int_to_positions(n)))


# =====================================================================
# CIRCLE PRIMES
# =====================================================================

def is_prime(n):
    """Standard primality test."""
    if n < 2:
        return False
    for d in range(2, int(n ** 0.5) + 1):
        if n % d == 0:
            return False
    return True


def find_circle_primes(limit):
    """Find all circle primes up to limit.

    A circle prime is an integer n > 1 that cannot be expressed as a ⊚ b
    for any 1 < a, b < n.
    """
    composite = set()
    for a in range(2, limit + 1):
        for b in range(a, limit + 1):
            p = knuth_product(a, b)
            if p <= limit:
                composite.add(p)
            else:
                break
    return [n for n in range(2, limit + 1) if n not in composite]


def circle_divide(n, d):
    """Find q such that d ⊚ q = n, or return None if no such q exists."""
    for q in range(1, n + 1):
        p = knuth_product(d, q)
        if p == n:
            return q
        if p > n:
            return None
    return None


def find_all_decompositions(n, circle_primes):
    """Find ALL factorisations of n into circle primes under ⊚."""
    cp_set = set(circle_primes)
    if n in cp_set:
        return [[n]]

    decomps = set()

    def _search(remaining, min_cp, current, depth):
        if depth > 25:
            return
        if remaining == 1:
            decomps.add(tuple(sorted(current)))
            return
        if remaining in cp_set and remaining >= min_cp:
            decomps.add(tuple(sorted(current + [remaining])))
        for cp in circle_primes:
            if cp < min_cp or cp > remaining:
                continue
            q = circle_divide(remaining, cp)
            if q is not None:
                if q == 1:
                    decomps.add(tuple(sorted(current + [cp])))
                elif q >= 2:
                    _search(q, cp, current + [cp], depth + 1)

    _search(n, 2, [], 0)
    return [list(d) for d in decomps]


# =====================================================================
# ANALYSIS
# =====================================================================

def verify_properties(max_assoc=5000, assoc_range=100):
    """Verify algebraic properties of the circle product.

    Returns True if all tests pass.
    """
    import random

    print("=" * 70)
    print("PROPERTY VERIFICATION")
    print("=" * 70)

    # Identity
    ok_id = all(knuth_product(1, n) == n for n in range(1, 501))
    print(f"\n1. Identity (1 ⊚ n = n, n ∈ [1, 500]):  {'✓' if ok_id else '✗'}")

    # Commutativity
    ok_comm = True
    for a in range(2, 101):
        for b in range(a + 1, 101):
            if knuth_product(a, b) != knuth_product(b, a):
                ok_comm = False
                break
        if not ok_comm:
            break
    print(f"2. Commutativity ([2, 100]):              {'✓' if ok_comm else '✗'}")

    # Associativity
    random.seed(42)
    assoc_failures = 0
    for _ in range(max_assoc):
        a = random.randint(2, assoc_range)
        b = random.randint(2, assoc_range)
        c = random.randint(2, assoc_range)
        pa, pb, pc = int_to_positions(a), int_to_positions(b), int_to_positions(c)
        left = positions_to_int(rep_multiply(rep_multiply(pa, pb), pc))
        right = positions_to_int(rep_multiply(pa, rep_multiply(pb, pc)))
        if left != right:
            assoc_failures += 1
    ok_assoc = assoc_failures == 0
    print(f"3. Associativity ({max_assoc} triples, [2, {assoc_range}]): "
          f"{'✓' if ok_assoc else '✗'} ({assoc_failures} failures)")

    # Distributivity
    d_pass = d_fail = 0
    for a in range(2, 25):
        for b in range(2, 25):
            for c in range(2, 25):
                if knuth_product(a, b + c) == knuth_product(a, b) + knuth_product(a, c):
                    d_pass += 1
                else:
                    d_fail += 1
    print(f"4. Distributivity ([2, 24]):               {'✓' if d_fail == 0 else '✗'} "
          f"({d_fail} non-distributive cases)")

    return ok_id and ok_comm and ok_assoc


def full_analysis(limit=500):
    """Run complete circle prime analysis up to limit."""
    print("\n" + "=" * 70)
    print(f"CIRCLE PRIME ANALYSIS (up to {limit})")
    print("=" * 70)

    cprimes = find_circle_primes(limit)
    oprimes = [p for p in range(2, limit + 1) if is_prime(p)]

    print(f"\nOrdinary primes π({limit}): {len(oprimes)}")
    print(f"Circle primes  C({limit}): {len(cprimes)}")

    o_set, c_set = set(oprimes), set(cprimes)
    both = sorted(o_set & c_set)
    only_o = sorted(o_set - c_set)
    only_c = sorted(c_set - o_set)

    print(f"\nBoth ordinary AND circle prime:     {len(both):3d}")
    print(f"  {both}")
    print(f"\nOrdinary prime, ⊚-composite:        {len(only_o):3d}")
    for p in only_o[:10]:
        for a in range(2, p):
            for b in range(a, p):
                if knuth_product(a, b) == p:
                    print(f"  {p:4d} = {a} ⊚ {b}")
                    break
            else:
                continue
            break
    if len(only_o) > 10:
        print(f"  ... and {len(only_o) - 10} more")

    print(f"\n⊚-irreducible ordinary composite:   {len(only_c):3d}")
    for n in only_c[:10]:
        f, m = [], n
        for p in range(2, m + 1):
            while m % p == 0:
                f.append(p)
                m //= p
        z = sorted(int_to_positions(n))
        print(f"  {n:4d} = {'×'.join(map(str, f)):12s}  φ-positions: {z}")
    if len(only_c) > 10:
        print(f"  ... and {len(only_c) - 10} more")

    # Unique factorisation
    test_lim = min(limit, 200)
    print(f"\n{'='*70}")
    print(f"UNIQUE FACTORISATION (2..{test_lim})")
    print(f"{'='*70}")

    non_unique = []
    no_decomp = []
    unique = 0

    for n in range(2, test_lim + 1):
        d = find_all_decompositions(n, cprimes)
        if len(d) == 0:
            no_decomp.append(n)
        elif len(d) > 1:
            non_unique.append((n, d))
        else:
            unique += 1

    print(f"\n  Unique:    {unique}")
    print(f"  No decomp: {len(no_decomp)}")
    print(f"  Non-unique: {len(non_unique)}")

    if no_decomp:
        print(f"\n  No ⊚-prime decomposition: {no_decomp}")
    if non_unique:
        print(f"\n  Non-unique factorisations:")
        for n, d in non_unique[:15]:
            print(f"    {n}: {d}")

    # Distribution
    print(f"\n{'='*70}")
    print("DISTRIBUTION")
    print(f"{'='*70}")
    print(f"\n{'N':>6} | {'π(N)':>5} | {'C(N)':>5} | {'π/C':>6}")
    print("-" * 35)
    for b in [10, 20, 50, 100, 200, 300, 500]:
        if b > limit:
            break
        pi = len([p for p in oprimes if p <= b])
        ci = len([p for p in cprimes if p <= b])
        print(f"{b:6d} | {pi:5d} | {ci:5d} | {pi / ci if ci else 0:6.3f}")

    # Mod 5
    print(f"\nMod 5 distribution:")
    for r in range(5):
        oc = len([p for p in oprimes if p % 5 == r])
        cc = len([p for p in cprimes if p % 5 == r])
        label = {0: "ramified", 1: "split", 2: "inert", 3: "inert", 4: "split"}[r]
        print(f"  ≡{r}: ordinary={oc:3d}, circle={cc:3d}  ({label} in Z[φ])")

    # Zeckendorf structure
    print(f"\nFirst 20 circle primes — Zeckendorf:")
    for cp in cprimes[:20]:
        z = sorted(int_to_positions(cp))
        zeck_idx = [p + 2 for p in z]
        print(f"  {cp:4d} = {' + '.join(f'F_{k}' for k in zeck_idx)}")

    has_f2 = sum(1 for cp in cprimes if 0 in int_to_positions(cp))
    print(f"\nCircle primes containing F_2: {has_f2}/{len(cprimes)}")

    return cprimes, oprimes


def product_table(lo=2, hi=12):
    """Print the circle product multiplication table."""
    print(f"\n   ⊚  |", end="")
    for b in range(lo, hi + 1):
        print(f" {b:5d}", end="")
    print()
    print("   " + "-" * (6 * (hi - lo + 1) + 4))
    for a in range(lo, hi + 1):
        print(f"   {a:2d} |", end="")
        for b in range(lo, hi + 1):
            print(f" {knuth_product(a, b):5d}", end="")
        print()


# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║   KNUTH CIRCLE PRODUCT — Full Analysis                             ║")
    print("║   Bootstrap Universe Programme — Paper 168 Appendix                 ║")
    print("║   Dr. Clifford Keeble, Woodbridge, UK                              ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    ok = verify_properties()
    product_table()
    if ok:
        full_analysis(500)
    else:
        print("\n⚠️  Property verification failed — results may be unreliable")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
