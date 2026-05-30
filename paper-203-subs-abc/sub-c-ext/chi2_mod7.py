"""
Cubic character chi_2 mod 7 (per brief Sub C-ext-2).

Definition: chi_2(a) = omega^k where k is the discrete log of a base 3
(primitive root) mod 7.
Specifically:
  a=1: chi_2(1) = 1 = omega^0
  a=2: chi_2(2) = omega^1
  a=3: chi_2(3) = omega^2   (3 is the primitive root)
  a=4: chi_2(4) = omega^2
  a=5: chi_2(5) = omega^1
  a=6: chi_2(6) = omega^0 = 1
  a=0: chi_2(0) = 0

with omega = exp(2*pi*i/3).
"""

from __future__ import annotations

import math
import cmath

# Cube root of unity
OMEGA = cmath.exp(2j * math.pi / 3)
OMEGA_2 = cmath.exp(4j * math.pi / 3)

CHI_2_MOD_7 = {
    0: 0+0j,
    1: 1+0j,
    2: OMEGA,
    3: OMEGA_2,
    4: OMEGA_2,
    5: OMEGA,
    6: 1+0j,
}


def chi2(a: int) -> complex:
    return CHI_2_MOD_7[a % 7]


def verify():
    out = {}
    # Sum over period
    total = sum(CHI_2_MOD_7.values())
    out['sum_over_period'] = total
    out['balanced'] = abs(total) < 1e-12

    # Multiplicativity
    mult_ok = True
    for a in range(7):
        for b in range(7):
            lhs = chi2(a * b)
            rhs = chi2(a) * chi2(b)
            if abs(lhs - rhs) > 1e-12:
                mult_ok = False
                break
        if not mult_ok:
            break
    out['multiplicative'] = mult_ok

    # Parity (even vs odd character): chi(-1) = chi(6) = 1 → EVEN
    out['chi_minus_1'] = chi2(-1)
    out['parity'] = 'even' if abs(chi2(-1) - 1) < 1e-12 else 'odd'

    # Order (smallest n s.t. chi^n is trivial)
    # chi(3) = omega -> chi^3(3) = 1. So order is 3 (cubic).
    chi_cubed = chi2(3) ** 3
    out['chi(3)^3'] = chi_cubed
    out['order'] = 'cubic (3)' if abs(chi_cubed - 1) < 1e-12 else '?'

    # Primitive: not factorable through smaller modulus
    # For conductor 7 (prime), the only non-trivial characters are mod 7, all primitive.
    out['conductor'] = 7
    return out


if __name__ == '__main__':
    print("=" * 60)
    print("chi_2 mod 7 (cubic character) verification")
    print("=" * 60)
    print("\nTable:")
    for a in range(7):
        c = CHI_2_MOD_7[a]
        print(f"  chi_2({a}) = {c.real:+.4f} + {c.imag:+.4f}i")
    print("\nVerification:")
    r = verify()
    for k, v in r.items():
        print(f"  {k}: {v}")
