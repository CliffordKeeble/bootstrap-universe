# character_projections.py - mod-5 Dirichlet character projections P(chi_k)
# Bootstrap Universe Programme - Paper 196 candidate, Layer 3a
#
# Per Mr Code Brief Addendum 2 (3 May 2026, CinC). Builds on
# hurwitz_120_mod5.S(s, r) to define
#
#   P(chi_k)(s) = sum_{r=1}^{4} conj(chi_k(r)) * S(s, r)
#
# for the four Dirichlet characters mod 5:
#   chi_0  trivial
#   chi_1  Legendre (n/5), real, order 2
#   chi_2  order 4 with chi_2(2) = +i
#   chi_3  conjugate of chi_2, order 4 with chi_3(2) = -i
#
# (Conjugation in the projection formula is per Paper 125 sec.3 / standard
# Dirichlet projection convention. Without it, signs of the order-4 phase
# would be flipped and the phase-lock test would land on +arctan(1/phi)
# instead of -arctan(1/phi).)
#
# Sanities in __main__:
#   1. Closed form: P(chi_0)(2) = 2304 * pi^2  (derivation in source).
#   2. Real-s typing at s=3: P(chi_0), P(chi_1) real;
#      P(chi_2), P(chi_3) complex-conjugate pair.
#   3. Zero behaviour: |P(chi_0)(rho_1)| < 1e-40 (vanishes with zeta(s)).
#   4. KEY (gating): P(chi_2)(conj(s)) = conj(P(chi_3)(s)) at rho_1, rho_2.
#      The phase-lock test depends on this relation; if it fails, stop.

import mpmath as mp
from mpmath import mp as _mp_ctx, mpc, mpf, conj as mp_conj

from hurwitz_120_mod5 import S


# Character table mod 5. Values are exact (1, -1, +i, -i) so Python complex
# literals are precision-safe; mpmath promotes to mpc on multiplication.
#
#   r=1  r=2  r=3  r=4
#   ---- ---- ---- ----
#  chi_0:   1    1    1    1
#  chi_1:   1   -1   -1    1     (Legendre symbol (r/5))
#  chi_2:   1   +i   -i   -1     (order 4, chi_2(2) = +i)
#  chi_3:   1   -i   +i   -1     (order 4, chi_3 = conj(chi_2))
CHI_MOD5 = {
    0: {1: 1,    2:  1,   3:  1,   4:  1},
    1: {1: 1,    2: -1,   3: -1,   4:  1},
    2: {1: 1,    2:  1j,  3: -1j,  4: -1},
    3: {1: 1,    2: -1j,  3:  1j,  4: -1},
}


def P(s, k):
    """Dirichlet character projection P(chi_k)(s) at N=120, mod 5.

        P(chi_k)(s) = sum_{r=1}^{4}  conj(chi_k(r)) * S(s, r)

    Parameters
    ----------
    s : complex (mpc preferred); for Layer 3a typically s = 1/2 + i*gamma.
    k : int in {0, 1, 2, 3}, character index.

    Returns
    -------
    mpc
    """
    if k not in (0, 1, 2, 3):
        raise ValueError("k must be in {0, 1, 2, 3}; got %r" % (k,))
    chi = CHI_MOD5[k]
    total = mpc(0)
    for r in (1, 2, 3, 4):
        total += complex(chi[r]).conjugate() * S(s, r)
    return total


# -----------------------------------------------------------------------------
# Sanity machinery
# -----------------------------------------------------------------------------

def _sanity_chi0_closed_form(tol=mpf(10) ** -30):
    """P(chi_0)(2) = 2304 * pi^2.

    Derivation:
        P(chi_0)(s) = sum_r S(s, r)  (since chi_0(r) = 1 for r coprime to 5)
        From hurwitz_120_mod5: 120^(-s) * sum_r S(s, r) = (1 - 5^(-s)) * zeta(s)
        => P(chi_0)(s) = 120^s * (1 - 5^(-s)) * zeta(s)
        At s=2: P(chi_0)(2) = 14400 * (24/25) * pi^2/6
                            = (14400 * 24) / (25 * 6) * pi^2
                            = 345600 / 150 * pi^2
                            = 2304 * pi^2
    """
    val = P(mpf(2), 0)
    expected_real = mpf(2304) * mp.pi ** 2
    re_err = abs(val.real - expected_real)
    im_err = abs(val.imag)
    return val, expected_real, re_err, im_err, (re_err < tol and im_err < tol)


def _sanity_real_s_typing(s_real=mpf(3), tol=mpf(10) ** -30):
    """At real s, P(chi_0), P(chi_1) are real; P(chi_2), P(chi_3) are conjugates."""
    p0 = P(s_real, 0)
    p1 = P(s_real, 1)
    p2 = P(s_real, 2)
    p3 = P(s_real, 3)
    im0 = abs(p0.imag)
    im1 = abs(p1.imag)
    re_diff = abs(p2.real - p3.real)
    im_sum = abs(p2.imag + p3.imag)
    passed = (im0 < tol and im1 < tol and re_diff < tol and im_sum < tol)
    return (p0, p1, p2, p3), (im0, im1, re_diff, im_sum), passed


def _sanity_chi0_at_zero(rho, tol=mpf(10) ** -40):
    """|P(chi_0)(rho)| < tol  (rho a zeta zero)."""
    val = P(rho, 0)
    mag = abs(val)
    return val, mag, (mag < tol)


def _sanity_conjugate_pair(s, tol=mpf(10) ** -30):
    """KEY GATING SANITY: P(chi_2)(conj(s)) = conj(P(chi_3)(s)).

    Holds because chi_2 = conj(chi_3) and the Hurwitz argument a/120 is real
    (so S(conj(s), r) = conj(S(s, r)) by Schwarz reflection). With both
    factors conjugating exactly, residuals are *algebraically* zero, not
    just at the dps roundoff floor - this is checked numerically below.

    The phase-lock test depends on this relation; if it fails, the arg ratio
    of P(chi_2)/P(chi_3) at zeta zeros has no defined symmetry to lock to.

    Test points include both critical-line zeros (rho_1, rho_2) and an
    off-line point (s = 2 + i). The off-line check is deliberate: a passing
    on-line test could be coincidence specific to the critical line; the
    off-line check confirms the relation is structural, not positional.
    """
    s = mpc(s)
    s_bar = mp_conj(s)
    lhs = P(s_bar, 2)
    rhs = mp_conj(P(s, 3))
    residual = abs(lhs - rhs)
    return lhs, rhs, residual, (residual < tol)


if __name__ == "__main__":
    _mp_ctx.dps = 50
    print("character_projections sanity - mp.dps = %d" % _mp_ctx.dps)
    print("=" * 72)
    print()

    # ---- Sanity 1 -----------------------------------------------------------
    print("Sanity 1: P(chi_0)(2) = 2304 * pi^2 closed form")
    val1, expected, re_err, im_err, passed1 = _sanity_chi0_closed_form()
    print("    P(chi_0)(2)   = %s" % mp.nstr(val1, 25))
    print("    2304 * pi^2   = %s" % mp.nstr(expected, 25))
    print("    |Re error|    = %s" % mp.nstr(re_err, 6))
    print("    |Im part|     = %s" % mp.nstr(im_err, 6))
    print("    [%s at tol 1e-30]" % ("PASS" if passed1 else "FAIL"))
    print()

    # ---- Sanity 2 -----------------------------------------------------------
    print("Sanity 2: real-s typing at s = 3")
    (p0, p1, p2, p3), (im0, im1, re_diff, im_sum), passed2 = _sanity_real_s_typing()
    print("    P(chi_0)(3)   = %s" % mp.nstr(p0, 20))
    print("    P(chi_1)(3)   = %s" % mp.nstr(p1, 20))
    print("    P(chi_2)(3)   = %s" % mp.nstr(p2, 20))
    print("    P(chi_3)(3)   = %s" % mp.nstr(p3, 20))
    print("    |Im P(chi_0)(3)|              = %s" % mp.nstr(im0, 6))
    print("    |Im P(chi_1)(3)|              = %s" % mp.nstr(im1, 6))
    print("    |Re P(chi_2) - Re P(chi_3)|  = %s" % mp.nstr(re_diff, 6))
    print("    |Im P(chi_2) + Im P(chi_3)|  = %s" % mp.nstr(im_sum, 6))
    print("    [%s at tol 1e-30]" % ("PASS" if passed2 else "FAIL"))
    print()

    # ---- Sanity 3 -----------------------------------------------------------
    rho1 = _mp_ctx.zetazero(1)
    print("Sanity 3: P(chi_0)(rho_1) vanishes")
    val3, mag3, passed3 = _sanity_chi0_at_zero(rho1)
    print("    rho_1         = %s" % mp.nstr(rho1, 25))
    print("    P(chi_0)(rho_1) = %s" % mp.nstr(val3, 12))
    print("    |P(chi_0)(rho_1)| = %s" % mp.nstr(mag3, 6))
    print("    [%s at tol 1e-40]" % ("PASS" if passed3 else "FAIL"))
    print()

    # ---- Sanity 4 (gating) --------------------------------------------------
    print("Sanity 4 (KEY/GATING): P(chi_2)(conj(s)) = conj(P(chi_3)(s))")
    rho2 = _mp_ctx.zetazero(2)
    sanity4_results = []
    for label, s_test in [("rho_1", rho1), ("rho_2", rho2),
                          ("2 + i (off-line check)", mpc(2, 1))]:
        lhs, rhs, residual, passed = _sanity_conjugate_pair(s_test)
        sanity4_results.append((label, residual, passed))
        print("    s = %s" % label)
        print("      LHS = P(chi_2)(conj(s)) = %s" % mp.nstr(lhs, 18))
        print("      RHS = conj(P(chi_3)(s)) = %s" % mp.nstr(rhs, 18))
        print("      |residual| = %s   [%s at tol 1e-30]"
              % (mp.nstr(residual, 6), "PASS" if passed else "FAIL"))
    passed4 = all(r[2] for r in sanity4_results)
    print()

    # ---- Summary ------------------------------------------------------------
    print("=" * 72)
    print("Summary:")
    print("  Sanity 1 (P(chi_0)(2) closed form):       %s" % ("PASS" if passed1 else "FAIL"))
    print("  Sanity 2 (real-s typing):                  %s" % ("PASS" if passed2 else "FAIL"))
    print("  Sanity 3 (P(chi_0) vanishes at rho_1):     %s" % ("PASS" if passed3 else "FAIL"))
    print("  Sanity 4 (KEY conjugate-pair identity):    %s" % ("PASS" if passed4 else "FAIL"))
    overall = passed1 and passed2 and passed3 and passed4
    print()
    print("Overall: %s" % ("ALL PASS - cleared for commit 3" if overall
                            else "FAILURES PRESENT - DO NOT PROCEED to phase-lock test"))
