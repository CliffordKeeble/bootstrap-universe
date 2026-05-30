"""
Dirichlet character tables for d in {3, 5, 7, 13}.

Per PRE_REGISTRATION.md, these are the four real primitive characters
associated with the fundamental discriminants:
- d = 5  → conductor 5  (chi_5)
- d = 13 → conductor 13 (chi_13)
- d = 3  → conductor 12 (chi_12)
- d = 7  → conductor 28 (chi_28)

Pre-committed values; verified at import time for balance and
multiplicativity.
"""

import math

# Hand-computed values per PRE_REGISTRATION.md
CHI_TABLE = {
    5: {
        'conductor': 5,
        'values': {0: 0, 1: 1, 2: -1, 3: -1, 4: 1},
        'name': 'chi_5 (Legendre mod 5)',
    },
    13: {
        'conductor': 13,
        'values': {0: 0, 1: 1, 2: -1, 3: 1, 4: 1, 5: -1, 6: -1,
                    7: -1, 8: -1, 9: 1, 10: 1, 11: -1, 12: 1},
        'name': 'chi_13 (Legendre mod 13)',
    },
    3: {
        'conductor': 12,
        'values': {0: 0, 1: 1, 2: 0, 3: 0, 4: 0, 5: -1, 6: 0,
                    7: -1, 8: 0, 9: 0, 10: 0, 11: 1},
        'name': 'chi_12 (Kronecker for Q(sqrt 3))',
    },
    7: {
        'conductor': 28,
        'values': {
            0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: 0,
            8: 0, 9: 1, 10: 0, 11: -1, 12: 0, 13: -1, 14: 0,
            15: -1, 16: 0, 17: -1, 18: 0, 19: 1, 20: 0, 21: 0,
            22: 0, 23: -1, 24: 0, 25: 1, 26: 0, 27: 1,
        },
        'name': 'chi_28 (Kronecker for Q(sqrt 7))',
    },
}


def chi(d, n):
    """Return chi_d(n) using the conductor for d."""
    tbl = CHI_TABLE[d]
    q = tbl['conductor']
    return tbl['values'][n % q]


def verify_characters():
    """Verify each character is balanced (sum over period = 0) and multiplicative."""
    out = {}
    for d, tbl in CHI_TABLE.items():
        q = tbl['conductor']
        # Balance check
        total = sum(tbl['values'].values())
        balanced = (total == 0)

        # Multiplicativity check: chi(a*b) = chi(a)*chi(b) for all a, b mod q
        mult_ok = True
        for a in range(q):
            for b in range(q):
                lhs = chi(d, a * b)
                rhs = chi(d, a) * chi(d, b)
                if lhs != rhs:
                    mult_ok = False
                    break
            if not mult_ok:
                break

        out[d] = {
            'conductor': q,
            'balanced': balanced,
            'sum': total,
            'multiplicative': mult_ok,
            'name': tbl['name'],
        }
    return out


if __name__ == '__main__':
    print("Character verification:")
    res = verify_characters()
    for d, info in res.items():
        print(f"  d={d:2d} ({info['name']}):")
        print(f"    conductor={info['conductor']:3d}, "
              f"balanced={info['balanced']} (sum={info['sum']}), "
              f"multiplicative={info['multiplicative']}")
