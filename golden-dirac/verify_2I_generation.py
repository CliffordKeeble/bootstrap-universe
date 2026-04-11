# verify_2I_generation.py
# Verify whether Gamma_seed generates the binary icosahedral group 2I
# Gap closure for Paper 191 Result 2, identified by Mr Adversary
# 2I Universe Programme / Bootstrap Universe Programme
# Dr Clifford Keeble, Woodbridge UK, ORCID 0009-0003-6828-2155
# Written by Mr Code, April 11 2026
#
# Strategy: exact arithmetic in Q(sqrt(5)).
# Every matrix entry is a + b*sqrt(5) with a, b rational.
# Represent as (a, b) tuples. No sympy needed for core computation.

from fractions import Fraction
from collections import Counter

# ── Q(sqrt(5)) arithmetic ────────────────────────────────────────────────────

def qs_add(x, y):  return (x[0] + y[0], x[1] + y[1])
def qs_sub(x, y):  return (x[0] - y[0], x[1] - y[1])
def qs_neg(x):     return (-x[0], -x[1])

def qs_mul(x, y):
    a, b = x; c, d = y
    return (a*c + 5*b*d, a*d + b*c)

def qs_inv(x):
    a, b = x
    denom = a*a - 5*b*b
    return (a / denom, -b / denom)

def qs_to_float(x):
    return float(x[0]) + float(x[1]) * 2.2360679774997896

ZERO  = (Fraction(0), Fraction(0))
ONE   = (Fraction(1), Fraction(0))
ROOT5 = (Fraction(0), Fraction(1))
NEG1  = (Fraction(-1), Fraction(0))

# ── 2x2 matrix over Q(sqrt(5)) ───────────────────────────────────────────────

def mat(a, b, c, d):  return (a, b, c, d)

def mat_mul(M, N):
    a, b, c, d = M; e, f, g, h = N
    return mat(
        qs_add(qs_mul(a,e), qs_mul(b,g)), qs_add(qs_mul(a,f), qs_mul(b,h)),
        qs_add(qs_mul(c,e), qs_mul(d,g)), qs_add(qs_mul(c,f), qs_mul(d,h))
    )

def mat_inv(M):
    a, b, c, d = M
    det = qs_sub(qs_mul(a,d), qs_mul(b,c))
    di = qs_inv(det)
    return mat(qs_mul(d,di), qs_mul(qs_neg(b),di),
               qs_mul(qs_neg(c),di), qs_mul(a,di))

def mat_det(M):
    return qs_sub(qs_mul(M[0],M[3]), qs_mul(M[1],M[2]))

def mat_key(M):
    return tuple(v for entry in M for v in entry)

EYE = mat(ONE, ZERO, ZERO, ONE)

# ── Generators ────────────────────────────────────────────────────────────────

Gamma_seed = mat(ZERO, ONE, ROOT5, ZERO)
Gamma_adj  = mat(ZERO, qs_neg(ONE), ROOT5, ZERO)

# ══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("2I GENERATION VERIFICATION")
print("=" * 60)

# ── 1. Determinant analysis ──────────────────────────────────────────────────

det_seed = mat_det(Gamma_seed)
det_adj  = mat_det(Gamma_adj)

print(f"\ndet(Gamma_seed) = {det_seed[0]} + {det_seed[1]}*sqrt(5)")
print(f"det(Gamma_adj)  = {det_adj[0]} + {det_adj[1]}*sqrt(5)")
print(f"Numerical: det(seed) = {qs_to_float(det_seed):.6f}")
print(f"Numerical: det(adj)  = {qs_to_float(det_adj):.6f}")

print("\n2I requires unit quaternions: det = 1 in SU(2).")
print("det(Gamma_seed) = -sqrt(5) != +/-1.")
print("=> Gamma_seed^n has det = (-sqrt(5))^n, never returns to 1.")
print("=> <Gamma_seed> is INFINITE under matrix multiplication.")

# ── 2. Verify: compute seed^n for small n ────────────────────────────────────

print("\n-- Powers of Gamma_seed --")
power = EYE
for n in range(1, 13):
    power = mat_mul(power, Gamma_seed)
    d = mat_det(power)
    print(f"  seed^{n:2d}: det = {qs_to_float(d):>12.4f}  "
          f"= I? {mat_key(power) == mat_key(EYE)}")

# ── 3. Seed^2 = sqrt(5)*I (confirming Clifford, not group) ──────────────────

seed_sq = mat_mul(Gamma_seed, Gamma_seed)
print(f"\nseed^2 = [[{qs_to_float(seed_sq[0]):.4f}, {qs_to_float(seed_sq[1]):.4f}],"
      f" [{qs_to_float(seed_sq[2]):.4f}, {qs_to_float(seed_sq[3]):.4f}]]")
print(f"       = sqrt(5) * I_2  (Clifford relation, not group closure)")

# ── 4. What DOES close? Check seed * adj and adj * seed ──────────────────────

print("\n-- Products of seed and adj --")
sa = mat_mul(Gamma_seed, Gamma_adj)
as_ = mat_mul(Gamma_adj, Gamma_seed)
print(f"seed * adj = [[{qs_to_float(sa[0]):.4f}, {qs_to_float(sa[1]):.4f}],"
      f" [{qs_to_float(sa[2]):.4f}, {qs_to_float(sa[3]):.4f}]]")
print(f"adj * seed = [[{qs_to_float(as_[0]):.4f}, {qs_to_float(as_[1]):.4f}],"
      f" [{qs_to_float(as_[2]):.4f}, {qs_to_float(as_[3]):.4f}]]")
print(f"det(seed*adj) = {qs_to_float(mat_det(sa)):.4f}")

# ── 5. Normalised generators (det = 1) ───────────────────────────────────────
# To get 2I, we need elements of SL(2, Q(sqrt(5))[i]) — det exactly 1.
# Gamma_seed / 5^(1/4) would have det = 1, but 5^(1/4) is NOT in Q(sqrt(5)).
# This is the fundamental obstruction.

print("\n-- Can we normalise to det = 1? --")
print("Need: seed_normalised = Gamma_seed / 5^(1/4)")
print("det(seed/5^(1/4)) = (-sqrt(5)) / sqrt(sqrt(5))^2 = (-sqrt(5))/sqrt(5) = -1")
print("But 5^(1/4) is NOT in Q(sqrt(5)). It lives in Q(5^(1/4)).")
print("=> Cannot normalise within Q(sqrt(5)).")
print("=> 2I does not embed in GL(2, Q(sqrt(5))) as a FINITE subgroup")
print("   with the Gamma_seed structure.")

# ── 6. The correct statement about 2I and the golden Dirac ───────────────────

print("\n" + "=" * 60)
print("DIAGNOSIS")
print("=" * 60)

print("""
Gamma_seed generates an INFINITE group under matrix multiplication
because det(Gamma_seed) = -sqrt(5) is not a root of unity.

The binary icosahedral group 2I (order 120) consists of UNIT
quaternions with entries in Z[phi]. As 2x2 matrices, these live
in SU(2) — complex 2x2 matrices with det = 1. They cannot be
faithfully represented as 2x2 REAL matrices.

2I can be represented as:
  (a) 2x2 complex matrices in SU(2) with entries in Q(sqrt(5))[i]
  (b) 4x4 real matrices via the quaternion -> matrix embedding
  (c) A quotient: in PGL(2, Q(sqrt(5))), modding out scalars

What Gamma_seed IS:
  - An element of the Clifford algebra Cl(1,0) over Q(sqrt(5))
  - seed^2 = sqrt(5)*I is the CLIFFORD relation, not group closure
  - The golden Clifford algebra {G1, G2, G3} encodes the STRUCTURE
    of 2I (the symmetries), but the matrices themselves are not
    elements of 2I

Paper 191 Result 2 options:
  (A) Restate: the golden Clifford algebra is the matrix algebra
      in which 2I acts, not a group generated by seed directly
  (B) Prove via 4x4: use the quaternion representation of the
      actual 2I generators (s, t) as 4x4 real matrices over
      Q(sqrt(5)), verify |group| = 120 in GL(4, Q(sqrt(5)))
  (C) Prove via SU(2): represent 2I generators as 2x2 complex
      matrices with entries in Q(sqrt(5))[i], verify closure
""")

# ── 7. Quick check: do the ACTUAL 2I generators close to 120? ────────────────
# Use 4x4 real quaternion representation, exact arithmetic in Q(sqrt(5))

print("=" * 60)
print("BONUS: ACTUAL 2I GENERATORS (4x4 quaternion representation)")
print("=" * 60)

# 4x4 matrices over Q(sqrt(5))
# Stored as 16-tuple of Q(sqrt(5)) elements

def mat4(rows):
    """Create 4x4 matrix from list of 16 Q(sqrt(5)) elements."""
    return tuple(rows)

def mat4_mul(A, B):
    """Multiply two 4x4 matrices."""
    result = []
    for i in range(4):
        for j in range(4):
            s = ZERO
            for k in range(4):
                s = qs_add(s, qs_mul(A[4*i+k], B[4*k+j]))
            result.append(s)
    return tuple(result)

def mat4_key(M):
    return tuple(v for entry in M for v in entry)

EYE4 = tuple(ONE if i == j else ZERO for i in range(4) for j in range(4))

# Quaternion q = a + bi + cj + dk -> 4x4 real matrix
def quat_to_mat4(a, b, c, d):
    return mat4([
        a, qs_neg(b), qs_neg(c), qs_neg(d),
        b, a, qs_neg(d), c,
        c, d, a, qs_neg(b),
        d, qs_neg(c), b, a
    ])

# phi = (1+sqrt(5))/2 = 1/2 + (1/2)*sqrt(5)
PHI = (Fraction(1,2), Fraction(1,2))
# psi = 1 - phi = 1/2 - (1/2)*sqrt(5)
PSI = (Fraction(1,2), Fraction(-1,2))
HALF = (Fraction(1,2), Fraction(0))

# 2I generator s = (1/2)(1 + phi*i + psi*j)
# a=1/2, b=phi/2, c=psi/2, d=0
s_a = HALF
s_b = qs_mul(PHI, HALF)
s_c = qs_mul(PSI, HALF)
s_d = ZERO

S = quat_to_mat4(s_a, s_b, s_c, s_d)

# 2I generator t = (1/2)(phi + 1*i + psi*k)
# a=phi/2, b=1/2, c=0, d=psi/2
t_a = qs_mul(PHI, HALF)
t_b = HALF
t_c = ZERO
t_d = qs_mul(PSI, HALF)

T = quat_to_mat4(t_a, t_b, t_c, t_d)

# Verify unit quaternion: |s|^2 = a^2 + b^2 + c^2 + d^2 should = 1
s_norm = qs_add(qs_add(qs_mul(s_a,s_a), qs_mul(s_b,s_b)),
                qs_add(qs_mul(s_c,s_c), qs_mul(s_d,s_d)))
t_norm = qs_add(qs_add(qs_mul(t_a,t_a), qs_mul(t_b,t_b)),
                qs_add(qs_mul(t_c,t_c), qs_mul(t_d,t_d)))

print(f"\n|s|^2 = {qs_to_float(s_norm):.6f}  (should be 1)")
print(f"|t|^2 = {qs_to_float(t_norm):.6f}  (should be 1)")

# Generate 2I from {S, T} as 4x4 matrices
print("\nGenerating 2I from quaternion generators s, t (4x4 real)...")

def mat4_inv(M):
    """Invert 4x4 matrix. For unit quaternion, inv = transpose."""
    # For orthogonal matrices (unit quaternion rep), M^-1 = M^T
    result = []
    for i in range(4):
        for j in range(4):
            result.append(M[4*j+i])
    return tuple(result)

all_gens = [S, mat4_inv(S), T, mat4_inv(T)]

seen = set()
elements = []

k0 = mat4_key(EYE4)
seen.add(k0)
elements.append(EYE4)
queue = [EYE4]

while queue:
    elem = queue.pop(0)
    for gen in all_gens:
        new = mat4_mul(elem, gen)
        k = mat4_key(new)
        if k not in seen:
            seen.add(k)
            elements.append(new)
            queue.append(new)

            if len(elements) % 20 == 0:
                print(f"  [2I] {len(elements)} elements, queue {len(queue)}")

            if len(elements) > 200:
                print(f"  [2I] SAFETY: > 200 elements, stopping")
                queue = []
                break

count_2I = len(elements)
print(f"\n|<s, t>| = {count_2I}")

if count_2I == 120:
    print("PASS: The standard 2I generators produce exactly 120 elements")
    print("as 4x4 real matrices over Q(sqrt(5)).")
    print("The binary icosahedral group 2I is confirmed in GL(4, Q(sqrt(5))).")

    # Check element orders
    print("\n-- Element orders in 2I --")
    order_counts = Counter()
    for elem in elements:
        p = EYE4
        for n in range(1, 200):
            p = mat4_mul(p, elem)
            if mat4_key(p) == mat4_key(EYE4):
                order_counts[n] += 1
                break
        else:
            order_counts['inf'] += 1

    print("  Order : Count")
    for o, c in sorted((k,v) for k,v in order_counts.items() if isinstance(k, int)):
        print(f"  {o:5d} : {c}")
    if 'inf' in order_counts:
        print(f"   >200 : {order_counts['inf']}")

    # Expected for 2I:
    print("\n  Expected for 2I:")
    print("  Order 1:  1  (identity)")
    print("  Order 2:  1  (-I)")
    print("  Order 3: 20")
    print("  Order 4: 30")
    print("  Order 5: 24")
    print("  Order 6: 20")
    print("  Order 10: 24")
    print("  Total:  120")
else:
    print(f"Count = {count_2I}, expected 120.")

print("\n" + "=" * 60)
print("COMPLETE")
print("=" * 60)
