"""
spectral_2I.py - 2I-spectral Stieltjes constant for Paper 190 Phase 1 v0.2.

S_n = sum_{k <= n, k in S(12,20,30)} m_k / lambda_k,   lambda_k = k(k+2)

with m_k = [z^k] P(z) and Molien series
    P(z) = (1 - z^60) / ((1-z^12)(1-z^20)(1-z^30)).

Analytic prediction (from CinC + first-principles check):
- After cyclotomic cancellation, P(z) has a pole of order 2 at z=1 with
  leading coefficient 1/120; thus m_k ~ (k+1)/120 averaged.
- m_k = 0 for odd k (gcd(12,20,30) = 2; S subset 2N).
- Within S, m_k ~ k/60 for k large (density of S in N is 1/2).
- m_k / (k(k+2)) ~ 1/(60 k) on S; density-weighted sum gives
        S_N ~ (1/120) log N + gamma_{2I-spectral} + o(1).

Pre-Task-2 sanity-check procedure: compute S_N at multiple N, fit the
log slope. If the empirical slope agrees with 1/120 to ~3 digits, we
have confirmation; we then extract gamma_{2I-spectral} via Richardson.
If slope is wrong, HALT-on-(h) and skip the row per brief sec 4.

Precision ceiling on row (h): Method A converges as N^{-1/2} (no
functional equation / analytic continuation available the way Dedekind
has), so achievable precision on gamma_{2I-spectral} is ~ log10(N)/2.
At N = 10^7 that's ~3-4 digits raw, or ~10-15 with Richardson against
the 1/sqrt(N) asymptotic. This caps the CF length on row (h) at
roughly 20-30 partial quotients. Documented and reported.
"""

import numpy as np
from mpmath import mp, mpf, log as mlog, sqrt as msqrt, nstr, matrix, lu_solve


def compute_molien_coefficients(N_max):
    """Compute m_k = [z^k] P(z) for k = 0..N_max via vectorised power-series
    division.

    P(z) = N(z) / D(z) with N(z) = 1 - z^60,
    D(z) = (1 - z^12)(1 - z^20)(1 - z^30) (degree 62).

    Returns: numpy int64 array of length N_max+1.
    """
    d1 = np.zeros(13, dtype=np.int64); d1[0] = 1; d1[12] = -1
    d2 = np.zeros(21, dtype=np.int64); d2[0] = 1; d2[20] = -1
    d3 = np.zeros(31, dtype=np.int64); d3[0] = 1; d3[30] = -1
    D = np.convolve(np.convolve(d1, d2), d3)   # length 63, D[0]=1
    deg_D = len(D) - 1                          # 62

    N_poly = np.zeros(61, dtype=np.int64); N_poly[0] = 1; N_poly[60] = -1
    deg_N = len(N_poly) - 1                     # 60

    m = np.zeros(N_max + 1, dtype=np.int64)
    D_tail = D[1:]  # length 62

    # Initial block: k = 0 .. deg_D, handled with explicit cap
    for k in range(min(deg_D + 1, N_max + 1)):
        nk = N_poly[k] if k <= deg_N else 0
        s = nk
        for i in range(1, k + 1):
            s -= D[i] * m[k - i]
        m[k] = s

    # Vectorised recurrence for k > deg_D
    # P[k] = N[k] - sum_{i=1..62} D[i] * P[k-i]   (N[k]=0 for k > 60)
    for k in range(deg_D + 1, N_max + 1):
        # dot product D[1:63] . m[k-1:k-63:-1]
        m[k] = -np.dot(D_tail, m[k-1:k-1-deg_D:-1])

    return m


def compute_partial_sum_S_N_float(m, N):
    """Float64 partial sum S_N. Fast vectorised numpy."""
    N_eff = min(N, len(m) - 1)
    k_arr = np.arange(1, N_eff + 1, dtype=np.float64)
    terms = m[1:N_eff + 1].astype(np.float64) / (k_arr * (k_arr + 2.0))
    return float(np.sum(terms))


def compute_partial_sum_S_N_mp(m, N, dps=50):
    """High-precision partial sum at given dps (slow, for the working value)."""
    old_dps = mp.dps
    mp.dps = dps
    s = mpf(0)
    N_eff = min(N, len(m) - 1)
    for k in range(1, N_eff + 1):
        mk = int(m[k])
        if mk != 0:
            s += mpf(mk) / (mpf(k) * mpf(k + 2))
    mp.dps = old_dps
    return s


def fit_log_slope_richardson(N_values, S_values, mp_dps=50):
    """Fit S_N = c log N + gamma + a/sqrt(N) + b/N + ...
    Use 3-parameter fit if 3 points; 4-parameter if >=4."""
    old_dps = mp.dps
    mp.dps = mp_dps
    n = len(N_values)
    if n == 3:
        cols = 3   # c, gamma, a/sqrt(N)
    elif n >= 4:
        cols = 4   # c, gamma, a/sqrt(N), b/N
    else:
        raise ValueError("need at least 3 N values to fit slope + gamma + 1 correction")

    A = matrix(n, cols)
    b_vec = matrix(n, 1)
    for i, N in enumerate(N_values):
        A[i, 0] = mlog(mpf(N))
        A[i, 1] = mpf(1)
        A[i, 2] = mpf(1) / msqrt(mpf(N))
        if cols >= 4:
            A[i, 3] = mpf(1) / mpf(N)
        b_vec[i, 0] = mpf(S_values[i])

    if n == cols:
        x = lu_solve(A, b_vec)
    else:
        AT = A.T
        x = lu_solve(AT * A, AT * b_vec)

    mp.dps = old_dps
    return {"c": x[0], "gamma": x[1], "a_inv_sqrt_N": x[2],
            "b_inv_N": x[3] if cols >= 4 else None}


def partial_sum_sanity_check(N_list=None, verbose=True):
    """The (h) sanity-check entry point.

    Returns dict with empirical c, predicted c = 1/120, fitted gamma,
    and a verdict on whether to proceed.
    """
    if N_list is None:
        N_list = [10**4, 10**5, 10**6, 10**7]

    N_max = max(N_list)
    if verbose:
        print(f"[2I-spec] Computing Molien coefficients m_k for k = 0..{N_max}...")
    import time
    t0 = time.time()
    m = compute_molien_coefficients(N_max)
    t1 = time.time()
    if verbose:
        print(f"[2I-spec]   done in {t1-t0:.1f} s")
        print(f"[2I-spec]   m[0..29] = {m[:30].tolist()}")
        odd_nonzero = sum(1 for k in range(1, min(N_max, 1000), 2) if m[k] != 0)
        print(f"[2I-spec]   odd-k m_k all zero for k < 1000: {odd_nonzero == 0}")
        in_S = [k for k in range(1, 120) if m[k] > 0]
        print(f"[2I-spec]   k in [1,120] with m_k > 0 (in S): {in_S}")

    # Partial sums at each N
    S_list = []
    for N in N_list:
        t0 = time.time()
        S = compute_partial_sum_S_N_float(m, N)
        t1 = time.time()
        S_list.append(S)
        if verbose:
            print(f"[2I-spec]   S_{N:>8} = {S:.15f}  ({t1-t0:.2f} s, float64)")

    # Fit at float-equivalent precision (we use mpmath for the linear algebra)
    fit = fit_log_slope_richardson(N_list, S_list, mp_dps=30)
    c_emp = fit["c"]
    g_emp = fit["gamma"]
    c_pred = mpf(1) / mpf(120)

    if verbose:
        print(f"[2I-spec] Fit results")
        print(f"[2I-spec]   c (empirical)  = {nstr(c_emp, 12)}")
        print(f"[2I-spec]   c (1/120 pred) = {nstr(c_pred, 12)}")
        print(f"[2I-spec]   |c_emp - c_pred| / c_pred = {nstr(abs(c_emp - c_pred) / c_pred, 4)}")
        print(f"[2I-spec]   gamma (fit)    = {nstr(g_emp, 12)}")
        print(f"[2I-spec]   a/sqrt(N) coef = {nstr(fit['a_inv_sqrt_N'], 8)}")
        if fit['b_inv_N'] is not None:
            print(f"[2I-spec]   b/N coef       = {nstr(fit['b_inv_N'], 8)}")

    relative_error = abs(c_emp - c_pred) / c_pred
    slope_ok = relative_error < mpf("0.01")   # 1% tolerance

    return {
        "N_list": N_list,
        "S_list": S_list,
        "c_empirical": c_emp,
        "c_predicted_1_over_120": c_pred,
        "c_relative_error": relative_error,
        "gamma_2I_spectral": g_emp,
        "a_correction": fit["a_inv_sqrt_N"],
        "b_correction": fit["b_inv_N"],
        "slope_ok": bool(slope_ok),
        "m_array_length": len(m),
        "first_30_m": m[:30].tolist(),
    }


if __name__ == "__main__":
    print("=" * 70)
    print("2I-spectral sanity check (Task 2 row (h))")
    print("=" * 70)
    result = partial_sum_sanity_check()
    print()
    print("=" * 70)
    if result["slope_ok"]:
        print(f"VERDICT: empirical c agrees with 1/120 to <1%.")
        print(f"         gamma_{{2I-spectral}} = {nstr(result['gamma_2I_spectral'], 10)}")
        print(f"         Proceed with CF computation for row (h), capped at ~20-30")
        print(f"         partial quotients due to N^{{-1/2}} precision limit.")
    else:
        print(f"VERDICT: empirical c does NOT match 1/120 within 1% tolerance.")
        print(f"         HALT-on-(h). Document divergence behaviour, skip row (h).")
    print("=" * 70)
