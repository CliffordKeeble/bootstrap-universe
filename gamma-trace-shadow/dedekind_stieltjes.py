"""
Stieltjes constants of Dedekind zeta-functions of quadratic number fields.

Pre-reg: Mr_Code_Brief_Paper_190_Phase_1_v0_2.md (git commit 971938b).

Convention: gamma_0(K) is the constant term of [zeta_K(s) - Res/(s-1)] at s=1.
No (-1)^n/n! sign flip. This matches Mathar (arXiv:1106.5552).

For K = Q(sqrt d), with primitive Kronecker character chi_d of conductor q,
zeta_K(s) = zeta(s) * L(s, chi_d). Laurent expansion at s=1:
    zeta_K(s) = L(1,chi_d)/(s-1) + [gamma * L(1,chi_d) + L'(1,chi_d)] + O(s-1)
so
    gamma_0(K) = gamma * L(1,chi_d) + L'(1,chi_d).
    R          = L(1, chi_d).

This module provides:
- Method B (PRIMARY for high precision): closed-form
  gamma_0(K) = gamma * L(1,chi) + L'(1,chi)
  with L, L' evaluated via mpmath's Hurwitz zeta expansion.

- Method A (CROSS-CHECK at modest precision): direct partial-sum
  S_N - R*log(N) -> gamma_0(K), with Richardson extrapolation
  against the residual.

Method B is used to seed the CF computation (needs 2000 dps); Method A
provides ~20-digit confirmation that conventions/signs are consistent.
"""

from mpmath import mp, mpf, mpc, log, sqrt, pi, euler, fsum, zeta, diff, floor


# ---------------------------------------------------------------------------
# Kronecker characters (period q, indexed a=1..q with chi(q)=chi(0)=0)
# ---------------------------------------------------------------------------

def chi_5(a):
    """Kronecker symbol (a/5), period 5, real character of Q(sqrt 5)."""
    r = a % 5
    return {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}[r]


def chi_minus_3(a):
    """Kronecker symbol (a/-3), period 3, real character of Q(sqrt -3)."""
    r = a % 3
    return {0: 0, 1: 1, 2: -1}[r]


# Sanity-check characters at d=2, 3 (used to validate the Method A/B
# methodology against simpler real-quadratic fields).

def chi_8(a):
    """Kronecker char of discriminant 8, primitive real char mod 8 for Q(sqrt 2).
    +1 if a == +-1 mod 8; -1 if a == +-3 mod 8; 0 otherwise (gcd(a,8) > 1)."""
    r = a % 8
    return {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}[r]


def chi_12(a):
    """Kronecker char of discriminant 12, primitive real char mod 12 for Q(sqrt 3).
    +1 if a == +-1 mod 12; -1 if a == +-5 mod 12; 0 otherwise."""
    r = a % 12
    table = {0: 0, 1: 1, 2: 0, 3: 0, 4: 0, 5: -1,
             6: 0, 7: -1, 8: 0, 9: 0, 10: 0, 11: 1}
    return table[r]


# ---------------------------------------------------------------------------
# Method B: closed form via Hurwitz expansion of Dirichlet L
# ---------------------------------------------------------------------------

def L_chi(s, chi, q):
    """Dirichlet L(s, chi) via L(s, chi) = q^{-s} sum_{a=1..q} chi(a) zeta(s, a/q).

    Valid for s != 1; the Hurwitz-zeta poles at s=1 cancel because
    sum_a chi(a) = 0 for non-principal chi.
    """
    return mpf(q) ** (-s) * fsum(
        chi(a) * zeta(s, mpf(a) / q) for a in range(1, q + 1) if chi(a) != 0
    )


def L_prime_at_1(chi, q):
    """L'(1, chi) via mpmath's numerical differentiation of L_chi at s=1.

    For real non-principal chi the function is analytic at s=1, so direct
    central differentiation is safe at sufficient working precision.
    """
    return diff(lambda s: L_chi(s, chi, q), mpf(1))


def gamma_0_method_B(chi, q, L_value=None, name=""):
    """Compute gamma_0(K) = gamma * L(1, chi) + L'(1, chi) at ambient mp.dps."""
    L1 = L_chi(mpf(1), chi, q) if L_value is None else L_value
    Lp1 = L_prime_at_1(chi, q)
    g = +euler
    val = g * L1 + Lp1
    return {
        "method": "B (Dedekind factorisation, closed form)",
        "name": name,
        "L_at_1": L1,
        "L_prime_at_1": Lp1,
        "gamma": g,
        "gamma_0_K": val,
        "dps": mp.dps,
    }


# ---------------------------------------------------------------------------
# Method A: direct ideal-counting partial sum with Richardson extrapolation
# ---------------------------------------------------------------------------

def ideal_counting_partial_sum(chi, q, R, N_list):
    """Compute S_N - R log N at multiple N values for Richardson extrapolation.

    a_n = (1 * chi)(n) = sum_{d | n} chi(d) gives the ideal-counting
    function. The partial sum S_N = sum_{n=1..N} a_n / n converges to
    R*log(N) + gamma_0(K) + lower-order terms.

    For efficiency at modest precision, we work at the prevailing mp.dps
    (typically 30-50 dps for this cross-check) and sum a_n/n via the
    Dirichlet hyperbola trick: a_n/n = sum_{de=n} chi(d)/n = chi(d)/(de),
    so sum_{n<=N} a_n/n = sum_{d<=N} chi(d)/d * sum_{e<=N/d} 1/e
                       = sum_{d<=N} chi(d)/d * H_{floor(N/d)}.
    """
    N_max = max(N_list)
    # Precompute harmonic numbers H_1..H_{N_max} at ambient mp.dps
    H = [mpf(0)] * (N_max + 1)
    acc = mpf(0)
    for k in range(1, N_max + 1):
        acc += mpf(1) / mpf(k)
        H[k] = acc

    # Precompute partial sums of chi(d)/d up to N_max
    # We'll need partial sums for arbitrary N, but here we compute one big
    # one and read off per-N values.
    chi_over_d_partial = [mpf(0)] * (N_max + 1)
    acc2 = mpf(0)
    for d in range(1, N_max + 1):
        if chi(d) != 0:
            acc2 += mpf(chi(d)) / mpf(d)
        chi_over_d_partial[d] = acc2

    results = []
    for N in N_list:
        # S_N = sum_{d=1..N} chi(d)/d * H_{floor(N/d)}
        s = mpf(0)
        for d in range(1, N + 1):
            cd = chi(d)
            if cd == 0:
                continue
            s += mpf(cd) / mpf(d) * H[N // d]
        diff_val = s - R * log(mpf(N))
        results.append({"N": N, "S_N": s, "S_N_minus_R_log_N": diff_val})
    return results


def gamma_0_method_A(chi, q, R, N_list=None, name=""):
    """Method A cross-check at modest precision.

    Default N_list: [10^4, 10^5, 10^6]. Returns the (N, S_N - R log N)
    sequence so the user can inspect convergence and report Richardson.
    """
    if N_list is None:
        N_list = [10**4, 10**5, 10**6]
    partials = ideal_counting_partial_sum(chi, q, R, N_list)
    # Richardson-style: assume S_N - R log N = gamma_0(K) + c1/sqrt(N) + ...
    # With three data points, fit (gamma_0, c1, c2) via the model
    # diff_N = g + c1/sqrt(N) + c2/N, solving the 3x3 linear system.
    # Robust enough for ~10 digit cross-check.
    # (We do not over-claim; this is a sanity check, not a precision tool.)
    if len(partials) >= 3:
        # Set up linear system
        from mpmath import matrix, lu_solve
        A = matrix(3, 3)
        b = matrix(3, 1)
        for i, p in enumerate(partials[-3:]):
            Ni = p["N"]
            A[i, 0] = mpf(1)
            A[i, 1] = mpf(1) / sqrt(mpf(Ni))
            A[i, 2] = mpf(1) / mpf(Ni)
            b[i, 0] = p["S_N_minus_R_log_N"]
        x = lu_solve(A, b)
        gamma_0_est = x[0]
    else:
        gamma_0_est = partials[-1]["S_N_minus_R_log_N"]

    return {
        "method": "A (ideal-sum partial + Richardson)",
        "name": name,
        "partial_sums": partials,
        "gamma_0_K_estimate": gamma_0_est,
        "dps": mp.dps,
    }


# ---------------------------------------------------------------------------
# Wrapper: both methods, returning Method B's high-precision value as the
# CF input, plus the agreement check.
# ---------------------------------------------------------------------------

def gamma_0_quadratic_field(chi, q, R_func, name=""):
    """Compute gamma_0(K) for K = Q(sqrt d) via both methods.

    chi:    Kronecker character function (takes int, returns -1/0/1)
    q:      conductor
    R_func: lambda returning R = L(1, chi) at ambient mp.dps (analytic form)
    name:   human-readable label
    """
    # High-precision Method B at current dps
    R = R_func()
    method_B = gamma_0_method_B(chi, q, L_value=R, name=name)

    # Cross-check Method A at modest precision (save dps; restore after)
    old_dps = mp.dps
    mp.dps = 40
    R_modest = R_func()
    method_A = gamma_0_method_A(chi, q, R_modest, name=name)
    mp.dps = old_dps

    # Agreement check at Method A's precision
    g_A = method_A["gamma_0_K_estimate"]
    g_B = method_B["gamma_0_K_estimate"] if "gamma_0_K_estimate" in method_B else method_B["gamma_0_K"]
    diff_AB = abs(g_A - g_B)

    return {
        "name": name,
        "gamma_0_K_high_precision": method_B["gamma_0_K"],
        "L_at_1": method_B["L_at_1"],
        "L_prime_at_1": method_B["L_prime_at_1"],
        "method_A_estimate": g_A,
        "method_A_partials": method_A["partial_sums"],
        "method_B_value": method_B["gamma_0_K"],
        "agreement_AB_abs": diff_AB,
        "dps_high": mp.dps,
        "dps_A_check": 40,
    }
