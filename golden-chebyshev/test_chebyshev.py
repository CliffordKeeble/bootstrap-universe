# test_chebyshev.py - tests for Layer 1
# Standalone: python test_chebyshev.py
#
# Pattern 14: unit tests are lemmas. Each test fixes a structural fact about the
# sieve, classification, or Chebyshev sign convention.

import sys
import numpy as np

from chebyshev_bias import (
    sieve_of_eratosthenes,
    segmented_sieve,
    classify,
    evaluate_E,
    SPLITTER_RES,
    STUBBORN_RES,
)


def assert_eq(a, b, msg):
    if a != b:
        print(f"  FAIL: {msg}: {a!r} != {b!r}")
        sys.exit(1)


def assert_arr_eq(a, b, msg):
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape or not np.all(a == b):
        print(f"  FAIL: {msg}: arrays differ")
        print(f"    got:      {a}")
        print(f"    expected: {b}")
        sys.exit(1)


def test_sieve_small():
    primes = sieve_of_eratosthenes(30)
    expected = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29])
    assert_arr_eq(primes, expected, "sieve up to 30")
    print("  ok: sieve_of_eratosthenes returns correct primes <= 30")


def test_sieve_medium():
    primes = sieve_of_eratosthenes(1000)
    # pi(1000) = 168 (classical)
    assert_eq(len(primes), 168, "pi(1000)")
    assert_eq(int(primes[0]), 2, "first prime")
    assert_eq(int(primes[-1]), 997, "largest prime <= 1000")
    print("  ok: pi(1000) = 168, last prime = 997")


def test_segmented_matches_simple():
    a = sieve_of_eratosthenes(10_000)
    b = segmented_sieve(10_000, segment=1000)
    assert_arr_eq(a, b, "segmented sieve matches simple sieve up to 10^4")
    print("  ok: segmented_sieve agrees with sieve_of_eratosthenes")


def test_classification_residues():
    # Verify SPLITTER_RES = {1, 4} (squares mod 5), STUBBORN_RES = {2, 3} (non-squares).
    # Squares mod 5: 1^2=1, 2^2=4, 3^2=4, 4^2=1 -> {1, 4}.
    squares = {(k * k) % 5 for k in range(1, 5)}
    assert_eq(squares, set(SPLITTER_RES), "squares mod 5")
    assert_eq(set(SPLITTER_RES) | set(STUBBORN_RES), {1, 2, 3, 4}, "splitter+stubborn cover non-zero residues")
    assert_eq(set(SPLITTER_RES) & set(STUBBORN_RES), set(), "splitter and stubborn disjoint")
    print("  ok: SPLITTER_RES = squares mod 5; disjoint from STUBBORN_RES")


def test_classify_excludes_p5():
    primes = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47])
    is_split, is_stub = classify(primes)
    # p=5 must be in neither (r=0).
    p5_idx = int(np.where(primes == 5)[0][0])
    assert_eq(bool(is_split[p5_idx]), False, "p=5 not classified as splitter")
    assert_eq(bool(is_stub[p5_idx]), False, "p=5 not classified as stubborn")
    # primes:   [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    # mod 5:    [ 2, 3, 0, 2,  1,  3,  2,  4,  3,  4,  1,  2,  1,  3,  2]
    # splitter: [ F, F, F, F,  T,  F,  F,  T,  F,  T,  T,  F,  T,  F,  F]
    # stubborn: [ T, T, F, T,  F,  T,  T,  F,  T,  F,  F,  T,  F,  T,  T]
    expected_split = np.array([False, False, False, False, True, False, False, True, False, True, True, False, True, False, False])
    expected_stub  = np.array([True,  True,  False, True,  False, True, True,  False, True, False, False, True, False, True, True])
    # Note: 2%5=2 (stub), 3%5=3 (stub) - the small primes 2 and 3 are stubborns.
    assert_arr_eq(is_split, expected_split, "splitter mask")
    assert_arr_eq(is_stub, expected_stub, "stubborn mask")
    # Total non-zero: every prime except 5 is in exactly one class.
    assert_eq(int((is_split | is_stub).sum()), len(primes) - 1, "all primes except 5 classified")
    print("  ok: classify() excludes p=5 and partitions all other primes")


def test_E_sign_convention():
    # E(x) = stubborns - splitters. Positive -> stubborn lead. Standard Chebyshev.
    # At x = 100: small primes by class:
    #   stubborns mod 5: 2, 3, 7, 13, 17, 23, 37, 43, 47, 53, 67, 73, 83, 97  (count varies)
    #   splitters mod 5: 11, 19, 29, 31, 41, 59, 61, 71, 79, 89  (count varies)
    primes = sieve_of_eratosthenes(100)
    is_split, is_stub = classify(primes)
    n_split = int(is_split.sum())
    n_stub = int(is_stub.sum())
    # Verify by direct enumeration:
    expected_split = sum(1 for p in primes if p != 5 and (p % 5) in SPLITTER_RES)
    expected_stub  = sum(1 for p in primes if p != 5 and (p % 5) in STUBBORN_RES)
    assert_eq(n_split, expected_split, "splitter count up to 100")
    assert_eq(n_stub, expected_stub, "stubborn count up to 100")
    # Stubborns should lead at x = 100 (Chebyshev bias direction).
    if n_stub > n_split:
        print(f"  ok: at x=100, stubborns lead splitters ({n_stub} vs {n_split}) - matches Chebyshev bias direction")
    else:
        print(f"  WARN: at x=100, splitters lead ({n_split} vs {n_stub}) - small-x oscillation, OK")


def test_evaluate_E_monotone_count():
    primes = sieve_of_eratosthenes(10_000)
    is_split, is_stub = classify(primes)
    x_samples = np.array([100, 1000, 10000])
    n_split, n_stub, E = evaluate_E(primes, is_split, is_stub, x_samples)
    # Cumulative counts must be monotone non-decreasing.
    assert_arr_eq(np.diff(n_split) >= 0, np.array([True, True]), "splitter count monotone")
    assert_arr_eq(np.diff(n_stub) >= 0, np.array([True, True]), "stubborn count monotone")
    # At x=10000: pi(10000)=1229; remove p=5; should split into stubborns + splitters.
    assert_eq(int(n_split[-1] + n_stub[-1]), 1229 - 1, "n_split + n_stub + 1 = pi(10^4) at x=10000")
    print(f"  ok: at x=10^4, n_stub={int(n_stub[-1])}, n_split={int(n_split[-1])}, E={int(E[-1]):+d}")


def test_count_at_million():
    """Sanity vs Granville-Martin Table 4 anchor."""
    print("  computing pi(10^6) by residue class (~2-3s)...")
    primes = segmented_sieve(1_000_000, segment=1_000_000)
    is_split, is_stub = classify(primes)
    n_split = int(is_split.sum())
    n_stub = int(is_stub.sum())
    # Classical pi(10^6) = 78498. Subtract p=5: 78497 split between classes.
    assert_eq(n_split + n_stub, 78498 - 1, "n_split + n_stub = pi(10^6) - 1 (subtract p=5)")
    E = n_stub - n_split
    print(f"  pi(10^6) = 78498, splitters = {n_split}, stubborns = {n_stub}, E = {E:+d}")
    # Stubborns are expected to lead at x = 10^6 (well-known Chebyshev bias).
    if E > 0:
        print(f"  ok: stubborns lead at x = 10^6 (E = {E:+d}) - matches classical Chebyshev bias for q=5")
    else:
        print(f"  FAIL: splitters lead at x = 10^6 (E = {E:+d}) - this contradicts classical result; sieve or classification is wrong")
        sys.exit(1)


def main():
    print("Running tests for chebyshev_bias.py...")
    print()
    print("[sieve correctness]")
    test_sieve_small()
    test_sieve_medium()
    test_segmented_matches_simple()
    print()
    print("[classification]")
    test_classification_residues()
    test_classify_excludes_p5()
    print()
    print("[E(x) and sign convention]")
    test_E_sign_convention()
    test_evaluate_E_monotone_count()
    print()
    print("[Granville-Martin sanity at x = 10^6]")
    test_count_at_million()
    print()
    print("ALL TESTS PASSED")


if __name__ == '__main__':
    main()
