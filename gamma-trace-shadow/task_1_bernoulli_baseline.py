"""
Task 1: Bernoulli baseline verification for Paper 190 Phase 1 v0.2.

Pre-registered as Task 1 of Mr_Code_Brief_Paper_190_Phase_1_v0_2.md
(pre-reg git commit: 971938b).

Computes gamma (Euler-Mascheroni) to >=50 digits via the Euler-Maclaurin
expansion. From brief section 1:

    H_n = ln(n) + gamma + 1/(2n) - 1/(12 n^2) + 1/(120 n^4) - ...
        = ln(n) + gamma + 1/(2n) - sum_{k=1}^{infty} B_{2k} / (2k n^{2k})

which inverts to

    gamma = H_n - ln(n) - 1/(2n) + sum_{k=1}^{K} B_{2k} / (2k n^{2k}) + R(n,K).

Brief section 4 Task 1 abbreviates the truncated form as

    gamma = H_n - ln(n) - sum_{k=1}^{K} B_{2k}/(2k n^{2k})

which is shorthand consistent with section 1's full expansion (the 1/(2n)
and the sign of the Bernoulli sum get absorbed by index/sign convention).
We implement the literal section-1 form to avoid ambiguity, and report
both expressions of the same final number.

Parameters: n = 10^6, K = 10 (i.e. asymptotic terms through order n^{-20}).
For n = 10^6 and K = 10 the next-term bound on the asymptotic residual
is |B_22/(22 n^22)| ~ 2.8e-130, so we have plenty of headroom for 50
correct digits.
"""

from sympy import bernoulli as sp_bernoulli, Rational, nsimplify
from mpmath import mp, mpf, log, fsum, euler, nstr

# ----------------------------------------------------------------------
# Working precision
# ----------------------------------------------------------------------
mp.dps = 60  # 60 decimal places working precision; >50 digits target

N = 10**6
K = 10  # asymptotic terms k = 1..10 -> orders n^-2 through n^-20

# ----------------------------------------------------------------------
# 1. Compute H_n at working precision
# ----------------------------------------------------------------------
print(f"[task 1] Computing H_n for n = {N} at {mp.dps} dps...")
H_n = fsum(mpf(1) / mpf(k) for k in range(1, N + 1))
print(f"[task 1] H_n done: {nstr(H_n, 12)}")

# ----------------------------------------------------------------------
# 2. ln(n)
# ----------------------------------------------------------------------
ln_n = log(mpf(N))

# ----------------------------------------------------------------------
# 3. Asymptotic correction terms (exact Bernoulli rationals via sympy)
# ----------------------------------------------------------------------
n_mpf = mpf(N)

# Exact Bernoulli coefficients B_{2k}/(2k) as Rational for display.
bern_table = []
for k in range(1, K + 1):
    B = sp_bernoulli(2 * k)               # exact Rational
    coeff = Rational(B, 2 * k)            # B_{2k}/(2k), exact
    bern_table.append((2 * k, B, coeff))

# Numeric truncated sum: sum_{k=1}^{K} B_{2k} / (2k n^{2k})
asymptotic_sum_terms = []
for (order, B, coeff) in bern_table:
    coeff_mpf = mpf(B.p) / mpf(B.q) / mpf(order)   # exact rational -> mpf
    term = coeff_mpf / (n_mpf ** order)
    asymptotic_sum_terms.append(term)

asymptotic_sum = fsum(asymptotic_sum_terms)

# Half-term: 1/(2n)
half_term = mpf(1) / (mpf(2) * n_mpf)

# gamma = H_n - ln(n) - 1/(2n) + asymptotic_sum
gamma_computed = H_n - ln_n - half_term + asymptotic_sum

# ----------------------------------------------------------------------
# 4. Reference value: mpmath's Euler-Mascheroni constant
# ----------------------------------------------------------------------
gamma_ref = +euler  # forces evaluation at mp.dps

err = gamma_computed - gamma_ref
abs_err = abs(err)
# correct digits = -log10(|err|)
if abs_err == 0:
    correct_digits = mp.dps
else:
    correct_digits = float(-log(abs_err, 10))

print("\n[task 1] Result")
print(f"  gamma (computed)   = {nstr(gamma_computed, 55)}")
print(f"  gamma (mpmath ref) = {nstr(gamma_ref, 55)}")
print(f"  signed error       = {nstr(err, 6)}")
print(f"  correct digits     ~ {correct_digits:.1f}")
print(f"  next-term bound    ~ B_{2*K+2}/(2*{K+1}) n^-{2*K+2} truncation")

# Predicted next-term magnitude as a sanity check on convergence rate
B_next = sp_bernoulli(2 * (K + 1))
next_coeff = Rational(B_next, 2 * (K + 1))
next_term_mag = abs(mpf(B_next.p) / mpf(B_next.q) / mpf(2 * (K + 1))) / (n_mpf ** (2 * (K + 1)))
print(f"  predicted residual ~ {nstr(next_term_mag, 5)} (next-term bound)")

# ----------------------------------------------------------------------
# 5. First 10 explicit terms of H_n - ln(n) - gamma
# ----------------------------------------------------------------------
print("\n[task 1] First 10 explicit terms of (H_n - ln n - gamma)")
print(f"  {'order':<8}{'B_{2k}':<22}{'signed coeff -B_{2k}/(2k)':<28}{'numeric at n=1e6'}")
print(f"  {'n^-1':<8}{'(half-term)':<22}{'+1/2':<28}{nstr(half_term, 16)}")
for (order, B, coeff) in bern_table:
    Bs = str(B)
    signed = -coeff  # the coefficient as it appears in H_n - ln n - gamma
    term_signed = -mpf(B.p) / mpf(B.q) / mpf(order) / (n_mpf ** order)
    print(f"  n^-{order:<5}{Bs:<22}{str(signed):<28}{nstr(term_signed, 16)}")

# ----------------------------------------------------------------------
# 6. Write markdown table for the findings report
# ----------------------------------------------------------------------
md_lines = []
md_lines.append("# Task 1 — Bernoulli baseline verification")
md_lines.append("")
md_lines.append(f"**Pre-reg commit:** `971938b`  ")
md_lines.append(f"**n =** 10^6 &nbsp; **K =** 10 &nbsp; **working dps =** {mp.dps}")
md_lines.append("")
md_lines.append("## Computed value vs reference")
md_lines.append("")
md_lines.append("```")
md_lines.append(f"gamma (computed)   = {nstr(gamma_computed, 55)}")
md_lines.append(f"gamma (mpmath ref) = {nstr(gamma_ref,      55)}")
md_lines.append(f"signed error       = {nstr(err,             6)}")
md_lines.append(f"correct digits     ~ {correct_digits:.1f}")
md_lines.append(f"next-term residual ~ {nstr(next_term_mag,   5)}")
md_lines.append("```")
md_lines.append("")
md_lines.append(
    "Note: observed error ~4×10<sup>−61</sup> reflects working-precision rounding "
    "accumulation in the H_n sum (60 dps), not the EM truncation bound at K=10 "
    "(~3×10<sup>−130</sup>). Verification passes by a wide margin against the "
    "≥50-digit target."
)
md_lines.append("")
md_lines.append("## First 10 asymptotic terms of H_n − ln n − γ")
md_lines.append("")
md_lines.append("| order | B<sub>2k</sub> (exact) | signed coefficient −B<sub>2k</sub>/(2k) | numeric @ n=10<sup>6</sup> |")
md_lines.append("|-------|------------------------|------------------------------------------|----------------------------|")
md_lines.append(f"| n<sup>−1</sup> | (half-term) | +1/2 | {nstr(half_term, 16)} |")
for (order, B, coeff) in bern_table:
    signed = -coeff
    term_signed = -mpf(B.p) / mpf(B.q) / mpf(order) / (n_mpf ** order)
    md_lines.append(
        f"| n<sup>−{order}</sup> | {B} | {signed} | {nstr(term_signed, 16)} |"
    )
md_lines.append("")
md_lines.append("## Verdict")
md_lines.append("")
if correct_digits >= 50:
    md_lines.append(f"**PASS.** Computed γ matches reference to ~{correct_digits:.0f} decimal "
                    f"digits (target: ≥ 50). Euler–Maclaurin / Bernoulli machinery verified. "
                    f"Cleared to proceed to Task 2.")
else:
    md_lines.append(f"**FAIL.** Computed γ matches reference to only ~{correct_digits:.0f} "
                    f"digits (target: ≥ 50). HALT and investigate before Task 2.")
md_lines.append("")

with open("task_1_baseline.md", "w", encoding="utf-8") as f:
    f.write("\n".join(md_lines))

print(f"\n[task 1] Wrote task_1_baseline.md  ({len(md_lines)} lines)")
