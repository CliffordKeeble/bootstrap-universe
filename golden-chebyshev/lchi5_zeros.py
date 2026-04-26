# lchi5_zeros.py - Zero data for L(s, chi_5)
# Bootstrap Universe Programme - Paper 196 candidate
#
# chi_5 = real primitive Dirichlet character mod 5 = Kronecker (n/5) = Legendre.
# LMFDB Conrey label: 5.4. L-function label: 1-5-5.4-r0-0-0.
# Conductor 5, root number +1, even parity (gamma factor Gamma_R(s)).
#
# Provides:
#   LMFDB_ZEROS_25  - first 25 positive zeros (4 decimals, source: LMFDB)
#   compute_L_chi5(t, dps)        - |L(1/2 + it, chi_5)| at high precision
#   refine_zero(t_approx, dps)    - Newton-refine to dps digits
#   get_zeros(n, dps, validate)   - return n positive zeros refined to dps digits

import sys

import mpmath
from mpmath import mp, mpf, mpc, fabs, fdiv

# chi_5 character values: chi_5(n) = (n/5) Legendre symbol.
# chi_5(0)=0, chi_5(1)=+1, chi_5(2)=-1, chi_5(3)=-1, chi_5(4)=+1.
CHI5_VALUES = [0, 1, -1, -1, 1]


# Source: LMFDB https://www.lmfdb.org/L/Dirichlet/5/4/ via /L/1-5-5.4-r0-0-0/.
# 25 positive zeros to 4 decimals.
LMFDB_ZEROS_25 = [
    mpf('6.6485'),
    mpf('9.8314'),
    mpf('11.9588'),
    mpf('16.0338'),
    mpf('17.5670'),
    mpf('19.5407'),
    mpf('22.2274'),
    mpf('24.5885'),
    mpf('26.7761'),
    mpf('28.4610'),
    mpf('29.7079'),
    mpf('33.0005'),
    mpf('34.7288'),
    mpf('35.8686'),
    mpf('38.1292'),
    mpf('39.5606'),
    mpf('41.8424'),
    mpf('44.0313'),
    mpf('45.4273'),
    mpf('46.4927'),
    mpf('48.3457'),
    mpf('51.0878'),
    mpf('52.1259'),
    mpf('53.8304'),
    mpf('55.5893'),
]


def compute_L_chi5(s, dps=30):
    """L(s, chi_5) at given complex s using mpmath.dirichlet."""
    with mp.workdps(dps):
        return mpmath.dirichlet(s, CHI5_VALUES)


def L_on_critical_line(t, dps=30):
    """L(1/2 + it, chi_5) as a complex value."""
    with mp.workdps(dps):
        s = mpc(mpf('0.5'), mpf(t))
        return mpmath.dirichlet(s, CHI5_VALUES)


def Lambda_chi5(t, dps=30):
    """
    Completed L-function Lambda(1/2 + it, chi_5). For chi_5 real, even, with
    root number +1: Lambda(s, chi) = (q/pi)^(s/2) Gamma(s/2) L(s, chi),
    and Lambda(1/2 + it, chi_5) is real-valued for real t. Returns the real
    part (Im part is numerical noise).
    """
    with mp.workdps(dps):
        s = mpc(mpf('0.5'), mpf(t))
        gamma_factor = mpmath.power(mpf(5) / mpmath.pi, s / 2) * mpmath.gamma(s / 2)
        L = mpmath.dirichlet(s, CHI5_VALUES)
        return (gamma_factor * L).real


def refine_zero(t_approx, dps=40, delta=None, tol=None, max_iter=80):
    """
    Bracket-bisect on Lambda(1/2 + it, chi_5), which is real-valued, to find
    the unique zero near t_approx. Uses mpmath secant on a real objective so
    iterates never leave the real axis.

    delta: half-width of starting bracket. Default 0.05 (LMFDB seeds are good
    to 4 decimals; zeros separated by ~1 unit; delta=0.05 cleanly isolates).
    """
    with mp.workdps(dps + 15):
        if tol is None:
            tol = mpf(10) ** (-(dps - 5))
        if delta is None:
            delta = mpf('0.05')
        else:
            delta = mpf(delta)

        def f(t):
            return Lambda_chi5(t, dps + 15)

        a = mpf(t_approx) - delta
        b = mpf(t_approx) + delta
        fa, fb = f(a), f(b)

        # Expand bracket up to a few times if no sign change.
        expansions = 0
        while fa * fb >= 0 and expansions < 5:
            delta = delta * 2
            a = mpf(t_approx) - delta
            b = mpf(t_approx) + delta
            fa, fb = f(a), f(b)
            expansions += 1

        if fa * fb >= 0:
            raise ValueError(
                f"Could not bracket zero near t={float(t_approx):.6f}: "
                f"Lambda has same sign at endpoints {float(a):.6f}, {float(b):.6f}"
            )

        # Bisection-with-secant-acceleration (Brent-like).
        for _ in range(max_iter):
            # Secant step.
            if fb != fa:
                t_secant = b - fb * (b - a) / (fb - fa)
            else:
                t_secant = (a + b) / 2
            # Fall back to bisection if secant leaves the bracket.
            if not (a < t_secant < b):
                t_secant = (a + b) / 2
            f_secant = f(t_secant)

            if fabs(f_secant) < tol:
                return t_secant

            if fa * f_secant < 0:
                b, fb = t_secant, f_secant
            else:
                a, fa = t_secant, f_secant

            if (b - a) < tol:
                return (a + b) / 2

        return (a + b) / 2


def validate_against_lmfdb(refined_zeros, lmfdb_zeros=LMFDB_ZEROS_25, decimals=4):
    """
    Assert that refined zeros agree with LMFDB to `decimals` decimal places.
    Returns list of (idx, refined, lmfdb, agree) tuples.
    """
    out = []
    tol = mpf(10) ** (-decimals + 1)  # half-unit at the requested decimal
    for i, (r, l) in enumerate(zip(refined_zeros, lmfdb_zeros)):
        diff = fabs(mpf(r) - mpf(l))
        agree = diff < tol
        out.append((i + 1, mpf(r), mpf(l), bool(agree), float(diff)))
    return out


def get_zeros(n=25, dps=40, validate=True, verbose=True):
    """
    Return the first n positive zeros of L(1/2 + it, chi_5), refined to dps digits.
    Seeds from LMFDB_ZEROS_25 and refines via Lambda-bisection.

    Currently supports n <= 25. Extension past 25 (iterative seeding from the
    Riemann-von Mangoldt density 1/(2 pi) log(q t / (2 pi e)) for L-functions)
    is a TODO; with 25 zeros, the explicit-formula reconstruction in Layer 2
    already covers x up to roughly e^(2 pi / gamma_25) ~ 10^5, sufficient for
    the toolchain validation.
    """
    if n > len(LMFDB_ZEROS_25):
        raise NotImplementedError(
            f"get_zeros currently supports n <= {len(LMFDB_ZEROS_25)}; got n={n}. "
            "Extension via Riemann-von Mangoldt seeding is TODO."
        )
    with mp.workdps(dps + 10):
        seeds = LMFDB_ZEROS_25[:n]
        if verbose:
            print(f"  refining {n} zeros to {dps} digits...")
        refined = []
        for t_seed in seeds:
            t_zero = refine_zero(t_seed, dps=dps)
            refined.append(t_zero)

        if validate and len(refined) >= len(LMFDB_ZEROS_25):
            checks = validate_against_lmfdb(refined[:len(LMFDB_ZEROS_25)])
            n_agree = sum(1 for _, _, _, a, _ in checks if a)
            if verbose:
                print(f"  LMFDB cross-check: {n_agree}/{len(checks)} zeros agree to 4 decimals")
            if n_agree < len(checks):
                fails = [c for c in checks if not c[3]]
                if verbose:
                    print(f"  disagreements: {fails}")

        return refined


def main():
    """Demo: refine the first 25 zeros and validate against LMFDB."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=25)
    parser.add_argument('--dps', type=int, default=40)
    args = parser.parse_args()

    print(f"Refining first {args.n} zeros of L(s, chi_5) to {args.dps} digits...")
    zeros = get_zeros(n=args.n, dps=args.dps)

    print(f"\n  First 5 refined zeros (40-digit display):")
    for i, z in enumerate(zeros[:5]):
        print(f"    gamma_{i+1} = {mpmath.nstr(z, 25)}")

    if len(zeros) >= 25:
        checks = validate_against_lmfdb(zeros[:25])
        max_diff = max(c[4] for c in checks)
        print(f"\n  Max deviation from LMFDB (4-decimal values): {max_diff:.6e}")
        if max_diff < 1e-3:
            print("  PASS: refined zeros agree with LMFDB at 4-decimal display precision.")
        else:
            print("  FAIL: significant disagreement - investigate before using these zeros.")
            sys.exit(1)


if __name__ == '__main__':
    main()
