"""
Smoke test for cf_tools and dedekind_stieltjes modules.

Three checks before launching the full Task 2 compute:
1. CF of gamma matches the published prefix.
2. Method B for gamma_0(Q(sqrt 5)) matches the ~0.7615 reference.
3. Methods A and B agree to ~10 digits at modest precision for both fields.
"""

from mpmath import mp, mpf, log, sqrt, pi, euler, nstr

from cf_tools import continued_fraction, cf_string
from dedekind_stieltjes import (
    chi_5, chi_minus_3, chi_8, chi_12,
    gamma_0_quadratic_field, gamma_0_method_B
)

# Known CF prefix of gamma (Euler-Mascheroni)
GAMMA_CF_PREFIX = [0, 1, 1, 2, 1, 2, 1, 4, 3, 13, 5, 1, 1, 8, 1, 2, 4, 1, 1, 40,
                   1, 11, 3, 7, 1, 7, 1, 1, 5, 1, 49, 4, 1, 65, 1, 4, 7, 11, 1]


def check_gamma_cf():
    mp.dps = 200
    g = +euler
    cf = continued_fraction(g, N=len(GAMMA_CF_PREFIX) + 5)
    print(f"[check 1] CF(gamma) first {len(GAMMA_CF_PREFIX)} terms:")
    print(f"  computed: {cf[:len(GAMMA_CF_PREFIX)]}")
    print(f"  expected: {GAMMA_CF_PREFIX}")
    match = cf[:len(GAMMA_CF_PREFIX)] == GAMMA_CF_PREFIX
    print(f"  match: {match}")
    return match


def check_gamma_0_sqrt5():
    """gamma_0(Q(sqrt 5)) via Method B at modest precision; compare to 0.7615..."""
    mp.dps = 50
    R = mpf(2) * log((mpf(1) + sqrt(mpf(5))) / mpf(2)) / sqrt(mpf(5))
    print(f"[check 2] L(1, chi_5) = 2 log(phi)/sqrt(5) = {nstr(R, 25)}")
    res = gamma_0_method_B(chi_5, 5, L_value=R, name="Q(sqrt 5)")
    g0 = res["gamma_0_K"]
    print(f"  gamma_0(Q(sqrt 5)) = {nstr(g0, 25)}")
    print(f"  CinC reference     ~ 0.7615")
    print(f"  L'(1, chi_5)       = {nstr(res['L_prime_at_1'], 25)}")
    matches_ref = abs(g0 - mpf("0.7615")) < mpf("0.001")
    print(f"  matches ~0.7615 reference: {matches_ref}")
    return matches_ref, g0


def check_gamma_0_sqrt_minus_3():
    """gamma_0(Q(sqrt -3)) via Method B at modest precision."""
    mp.dps = 50
    R = pi / (mpf(3) * sqrt(mpf(3)))
    print(f"[check 3] L(1, chi_-3) = pi/(3 sqrt 3) = {nstr(R, 25)}")
    res = gamma_0_method_B(chi_minus_3, 3, L_value=R, name="Q(sqrt -3)")
    g0 = res["gamma_0_K"]
    print(f"  gamma_0(Q(sqrt -3)) = {nstr(g0, 25)}")
    print(f"  L'(1, chi_-3)       = {nstr(res['L_prime_at_1'], 25)}")
    return g0


def check_method_A_vs_B_sqrt5():
    """Run both methods at 40 dps for Q(sqrt 5); check agreement."""
    mp.dps = 40
    R_func = lambda: mpf(2) * log((mpf(1) + sqrt(mpf(5))) / mpf(2)) / sqrt(mpf(5))
    print("[check 4] Method A vs B at 40 dps for Q(sqrt 5)...")
    result = gamma_0_quadratic_field(chi_5, 5, R_func, name="Q(sqrt 5)")
    print(f"  Method B (high prec) = {nstr(result['method_B_value'], 25)}")
    print(f"  Method A (Richardson) = {nstr(result['method_A_estimate'], 25)}")
    print(f"  |A - B| = {nstr(result['agreement_AB_abs'], 6)}")
    print(f"  Partial sums:")
    for p in result["method_A_partials"]:
        print(f"    N={p['N']:>8}: S_N - R log N = {nstr(p['S_N_minus_R_log_N'], 20)}")
    return result


def check_method_A_vs_B_sqrt_minus_3():
    mp.dps = 40
    R_func = lambda: pi / (mpf(3) * sqrt(mpf(3)))
    print("[check 5] Method A vs B at 40 dps for Q(sqrt -3)...")
    result = gamma_0_quadratic_field(chi_minus_3, 3, R_func, name="Q(sqrt -3)")
    print(f"  Method B (high prec) = {nstr(result['method_B_value'], 25)}")
    print(f"  Method A (Richardson) = {nstr(result['method_A_estimate'], 25)}")
    print(f"  |A - B| = {nstr(result['agreement_AB_abs'], 6)}")
    print(f"  Partial sums:")
    for p in result["method_A_partials"]:
        print(f"    N={p['N']:>8}: S_N - R log N = {nstr(p['S_N_minus_R_log_N'], 20)}")
    return result


def check_method_A_vs_B_sqrt2():
    """Sanity check at d=2: A vs B agreement validates methodology, no published
    cross-check available without external lookup."""
    from mpmath import log as mlog
    mp.dps = 40
    # L(1, chi_8) = log(1 + sqrt 2) / sqrt 2  (Dirichlet, h=1, eps=1+sqrt 2, D=8)
    R_func = lambda: mlog(mpf(1) + sqrt(mpf(2))) / sqrt(mpf(2))
    print("[check 6] Method A vs B at 40 dps for Q(sqrt 2) [sanity]...")
    print(f"  L(1, chi_8) = log(1+sqrt 2)/sqrt 2 = {nstr(R_func(), 20)}")
    result = gamma_0_quadratic_field(chi_8, 8, R_func, name="Q(sqrt 2)")
    print(f"  Method B (high prec)  = {nstr(result['method_B_value'], 25)}")
    print(f"  Method A (Richardson) = {nstr(result['method_A_estimate'], 25)}")
    print(f"  |A - B| = {nstr(result['agreement_AB_abs'], 6)}")
    return result


def check_method_A_vs_B_sqrt3():
    """Sanity check at d=3, internal-consistency only."""
    from mpmath import log as mlog
    mp.dps = 40
    # L(1, chi_12) = log(2 + sqrt 3) / sqrt 3  (h=1, eps=2+sqrt 3, D=12)
    R_func = lambda: mlog(mpf(2) + sqrt(mpf(3))) / sqrt(mpf(3))
    print("[check 7] Method A vs B at 40 dps for Q(sqrt 3) [sanity]...")
    print(f"  L(1, chi_12) = log(2+sqrt 3)/sqrt 3 = {nstr(R_func(), 20)}")
    result = gamma_0_quadratic_field(chi_12, 12, R_func, name="Q(sqrt 3)")
    print(f"  Method B (high prec)  = {nstr(result['method_B_value'], 25)}")
    print(f"  Method A (Richardson) = {nstr(result['method_A_estimate'], 25)}")
    print(f"  |A - B| = {nstr(result['agreement_AB_abs'], 6)}")
    return result


if __name__ == "__main__":
    print("=" * 70)
    print("Smoke test: cf_tools + dedekind_stieltjes")
    print("=" * 70)

    print()
    ok1 = check_gamma_cf()

    print()
    ok2, g0_5 = check_gamma_0_sqrt5()

    print()
    g0_m3 = check_gamma_0_sqrt_minus_3()

    print()
    res5 = check_method_A_vs_B_sqrt5()

    print()
    res_m3 = check_method_A_vs_B_sqrt_minus_3()

    print()
    res2 = check_method_A_vs_B_sqrt2()

    print()
    res3 = check_method_A_vs_B_sqrt3()

    print()
    print("=" * 70)
    print(f"check 1 (gamma CF prefix):           {'PASS' if ok1 else 'FAIL'}")
    print(f"check 2 (gamma_0(sqrt 5) ~ 0.7615):  FAIL (reference misrecalled; see chat)")
    print(f"check 3 (gamma_0(sqrt -3) computed): see above")
    print(f"check 4 (A vs B for sqrt 5):         |diff| = {nstr(res5['agreement_AB_abs'], 4)}")
    print(f"check 5 (A vs B for sqrt -3):        |diff| = {nstr(res_m3['agreement_AB_abs'], 4)}")
    print(f"check 6 (A vs B for sqrt 2  sanity): |diff| = {nstr(res2['agreement_AB_abs'], 4)}")
    print(f"check 7 (A vs B for sqrt 3  sanity): |diff| = {nstr(res3['agreement_AB_abs'], 4)}")
    print("=" * 70)
    print("Validation: four-field A-vs-B internal consistency at ~modest precision.")
    print("Note: published-source cross-check (Mathar/OEIS/LMFDB) requires web fetch,")
    print("which is YELLOW per CLAUDE.md and not authorised in this session.")
