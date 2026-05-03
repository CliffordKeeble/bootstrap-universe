# hurwitz_120_mod5.py - N=120 Hurwitz mod-5 decomposition
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a
#
# Per Mr Code Brief Addendum 2 (3 May 2026, CinC). Provides the four
# residue-class partial sums
#
#   S(s, r) = sum over a in [1, 120] with a mod 5 == r of  zeta(s, a/120)
#
# for r in {1, 2, 3, 4}. Each S has exactly 24 Hurwitz-zeta terms; the
# 24 multiples of 5 in [1, 120] are excluded as non-coprime to 5. The four
# sums together span the 96 Hurwitz terms coprime to 5.
#
# Foundation: the Hurwitz identity
#
#   sum_{a=1}^{N} zeta(s, a/N) = N^s * zeta(s)
#
# applied with N=120 and the multiples-of-5 a-values excluded gives the
# principal-character Dirichlet L-function L(s, chi_0 mod 5):
#
#   120^(-s) * (S[1] + S[2] + S[3] + S[4]) = (1 - 5^(-s)) * zeta(s)
#
# This is the structural sanity test in __main__ below. (The brief's loose
# phrasing "reproduces zeta(s)" refers to the full N=120 sum including
# multiples of 5; with multiples excluded we land on L(s, chi_0), which
# differs from zeta(s) by the (1 - 5^(-s)) Euler factor at p=5. Same identity,
# stated precisely.)

import mpmath as mp
from mpmath import mp as _mp_ctx, mpf, mpc, zeta

# -----------------------------------------------------------------------------
# Precision precaution - load-bearing comment, do NOT "simplify" away.
# -----------------------------------------------------------------------------
# The Hurwitz argument MUST be constructed as mpf(a) / 120, NOT as the Python
# expression a / 120. The latter is float64 (~16 digits) and would silently
# cap the precision of every zeta(s, a/120) call at ~10^-16, regardless of how
# high mp.dps is set. The dps=50 pass threshold of Layer 3a (residual < 10^-35)
# is meaningless if this is wrong.
#
# Future Mr Code or anyone else: leave the mpf(a) / 120 construction alone
# unless you have re-verified the identity
#     120^(-s) * sum(S(s, r) for r in 1..4) = (1 - 5^(-s)) * zeta(s)
# at dps=50 with the change.
# -----------------------------------------------------------------------------

# Precompute the residue-class indexing once.
_RESIDUE_CLASSES = {r: [a for a in range(1, 121) if a % 5 == r] for r in (1, 2, 3, 4)}

# Structural sanity: 24 elements per class, 96 total.
assert all(len(_RESIDUE_CLASSES[r]) == 24 for r in (1, 2, 3, 4)), \
    "residue-class indexing broken"
assert sum(len(v) for v in _RESIDUE_CLASSES.values()) == 96, \
    "total coprime-to-5 count broken"


def S(s, r):
    """Hurwitz partial sum for residue class r mod 5 at N = 120.

    S(s, r) = sum of zeta(s, a/120) over the 24 values of a in [1, 120]
    with a == r (mod 5).

    Parameters
    ----------
    s : complex (mpmath mpc preferred) or real
        For Layer 3a typically s = 0.5 + i*gamma at zeta zeros.
    r : int in {1, 2, 3, 4}
        Residue class mod 5.

    Returns
    -------
    mpc
        The 24-term Hurwitz sum at the working mp.dps precision.
    """
    if r not in (1, 2, 3, 4):
        raise ValueError("r must be in {1, 2, 3, 4}; got %r" % (r,))
    s = mpc(s)
    total = mpc(0)
    for a in _RESIDUE_CLASSES[r]:
        # mpf(a) / 120 - see Precision precaution above.
        total += zeta(s, mpf(a) / 120)
    return total


def _principal_character_residual(s):
    """Return (lhs, rhs, |lhs - rhs|) for the principal-character identity.

        lhs = 120^(-s) * (S[1] + S[2] + S[3] + S[4])
        rhs = (1 - 5^(-s)) * zeta(s)

    Pass condition at dps=50 with tolerance 10^-30 for s away from zeros;
    near a zeta zero both sides go to 0, residual |lhs - rhs| stays small
    in absolute terms.
    """
    s = mpc(s)
    lhs = mpf(120) ** (-s) * sum(S(s, r) for r in (1, 2, 3, 4))
    rhs = (1 - mpf(5) ** (-s)) * zeta(s)
    return lhs, rhs, abs(lhs - rhs)


if __name__ == "__main__":
    _mp_ctx.dps = 50

    print("hurwitz_120_mod5 sanity - mp.dps = %d" % _mp_ctx.dps)
    print("=" * 72)
    print()
    print("Identity tested: 120^(-s) * (S[1]+S[2]+S[3]+S[4]) = (1 - 5^(-s)) * zeta(s)")
    print("(This is the principal-character Dirichlet L-function L(s, chi_0 mod 5).)")
    print()

    # Structural report.
    print("Residue-class indexing:")
    for r in (1, 2, 3, 4):
        idx = _RESIDUE_CLASSES[r]
        print("  r = %d: %d terms, a in [%d, %d, %d, ..., %d, %d]"
              % (r, len(idx), idx[0], idx[1], idx[2], idx[-2], idx[-1]))
    print()

    # Numerical sanity battery.
    test_points = [
        ("s = 2 (real, away from zeros)", mpf(2)),
        ("s = 3 (real, away from zeros)", mpf(3)),
        ("s = 2 + i (complex, away from zeros)", mpc(2, 1)),
        ("s = rho_1 = 0.5 + i*gamma_1 (first zeta zero)", _mp_ctx.zetazero(1)),
    ]

    tol = mpf(10) ** (-30)
    all_passed = True
    for label, s_test in test_points:
        lhs, rhs, residual = _principal_character_residual(s_test)
        passed = residual < tol
        all_passed &= passed
        print("  %s" % label)
        print("    LHS      = %s" % mp.nstr(lhs, 25))
        print("    RHS      = %s" % mp.nstr(rhs, 25))
        print("    |residual| = %s   [%s at tol 1e-30]"
              % (mp.nstr(residual, 6), "PASS" if passed else "FAIL"))
        print()

    print("=" * 72)
    print("Overall: %s" % ("ALL PASS" if all_passed else "FAILURES PRESENT"))
