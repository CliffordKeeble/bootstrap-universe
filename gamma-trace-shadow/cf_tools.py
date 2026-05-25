"""
Continued-fraction utilities for Paper 190 Phase 1 v0.2.

Pre-reg: Mr_Code_Brief_Paper_190_Phase_1_v0_2.md (git commit 971938b).

Provides:
- continued_fraction(x, N): bare CF computation.
- cf_with_precision_check(x_high, x_low, N): primary/cross precision CFs,
  returning the stable common prefix plus disagreement diagnostics per the
  refinement adopted from CinC: when a precision-stability failure occurs
  at step n_fail, we record BOTH the disagreement AND the precision-stable
  value at n_fail - 1.
- convergents(cf, m): first m convergents as exact Fractions.
"""

from fractions import Fraction
from mpmath import mp, mpf, floor


def continued_fraction(x, N=300):
    """Compute up to N partial quotients of x.

    Stops early if the fractional part becomes exactly zero (rare for our
    candidates - they are all irrational). Does NOT do its own precision
    check; for that, use cf_with_precision_check.
    """
    cf = []
    y = mpf(x)
    for _ in range(N):
        a = int(floor(y))
        cf.append(a)
        frac = y - a
        if frac == 0:
            break
        y = mpf(1) / frac
    return cf


def cf_with_precision_check(x_high, x_low, N=300):
    """Compute CFs of x_high and x_low in lockstep, return stable common
    prefix and disagreement diagnostics.

    Both iterations run at the ambient mp.dps. x_high should be the value
    computed at the higher input precision (e.g. 2000 dps), x_low at the
    cross-check precision (e.g. 1500 dps). Both are passed in at their
    natively-computed precision; mpmath padding at the higher iteration
    precision is harmless.

    Returns:
        cf_stable:   list of partial quotients where both inputs agree
        n_fail:      index of first disagreement (== N if all stable)
        info:        dict with disagreement details, including the value
                     at n_fail - 1 per the CinC refinement.
    """
    cf_stable = []
    yh = mpf(x_high)
    yl = mpf(x_low)

    for k in range(N):
        a_h = int(floor(yh))
        a_l = int(floor(yl))

        if a_h != a_l:
            info = {
                "k_fail": k,
                "a_high": a_h,
                "a_low": a_l,
                "frac_high_at_k_fail": float(yh - a_h),
                "frac_low_at_k_fail": float(yl - a_l),
                "a_stable_at_k_fail_minus_1": cf_stable[-1] if cf_stable else None,
                "stable_length": len(cf_stable),
            }
            return cf_stable, k, info

        cf_stable.append(a_h)

        frac_h = yh - a_h
        frac_l = yl - a_l
        if frac_h == 0 or frac_l == 0:
            info = {
                "k_fail": k + 1,
                "reason": "exact rational reached (shouldn't happen for irrationals)",
                "a_stable_at_k_fail_minus_1": cf_stable[-1],
                "stable_length": len(cf_stable),
            }
            return cf_stable, k + 1, info

        yh = mpf(1) / frac_h
        yl = mpf(1) / frac_l

    return cf_stable, N, None  # all N partial quotients stable


def convergents(cf, m=10):
    """Return the first m convergents (p_k/q_k) from a CF list as Fractions."""
    p_prev, p_curr = 0, 1
    q_prev, q_curr = 1, 0
    convs = []
    for a in cf[:m]:
        p_next = a * p_curr + p_prev
        q_next = a * q_curr + q_prev
        convs.append(Fraction(p_next, q_next))
        p_prev, p_curr = p_curr, p_next
        q_prev, q_curr = q_curr, q_next
    return convs


def cf_string(cf, max_display=300):
    """Pretty-print a CF as '[a_0; a_1, a_2, ...]'."""
    if not cf:
        return "[]"
    disp = cf[:max_display]
    head = str(disp[0])
    tail = ", ".join(str(a) for a in disp[1:])
    suffix = "" if len(cf) <= max_display else f", ... (+{len(cf) - max_display} more)"
    if tail:
        return f"[{head}; {tail}{suffix}]"
    return f"[{head}{suffix}]"
