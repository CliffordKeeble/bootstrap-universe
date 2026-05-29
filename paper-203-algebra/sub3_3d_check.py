"""
Sub 3 extra check: compute (Gamma_5)^2 for all 3-gamma signatures,
where Gamma_5 = Gamma_1 * Gamma_2 * Gamma_3. This is the Paper 191
chirality at D=3 — used to verify Mr A's circularity catch.
"""

import sympy as sp
from sympy import sqrt, Matrix, eye, zeros, kronecker_product, simplify
import itertools

phi = (1 + sqrt(5)) / 2
root5 = 2*phi - 1
I2 = eye(2)
I4 = eye(4)

K = Matrix([[0, 1], [1, 0]])
J = Matrix([[0, -1], [1, 0]])
Z = K * J
Gamma_seed = Matrix([[0, 1], [root5, 0]])
Gamma_adj  = Matrix([[0, -1], [root5, 0]])


def anticomm(A, B):
    return sp.simplify(A * B + B * A)


def is_zero(M):
    return simplify(M) == sp.zeros(*M.shape)


# Build pool of all A x B with squares ±sqrt(5)*I_4
basis = {'I': I2, 'K': K, 'J': J, 'Z': Z, 'Gseed': Gamma_seed, 'Gadj': Gamma_adj}
pool = []
for nA, A in basis.items():
    for nB, B in basis.items():
        G = kronecker_product(A, B)
        sq = simplify(G * G)
        if is_zero(sq - root5 * I4):
            pool.append((f"{nA}x{nB}", G, +1))
        elif is_zero(sq + root5 * I4):
            pool.append((f"{nA}x{nB}", G, -1))

# Find all 3-cliques (mutually anticommuting triples)
n = len(pool)
ac = [[False]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if is_zero(anticomm(pool[i][1], pool[j][1])):
            ac[i][j] = True
            ac[j][i] = True

triples = []
for combo in itertools.combinations(range(n), 3):
    if all(ac[a][b] for a, b in itertools.combinations(combo, 2)):
        triples.append(combo)

print(f"Found {len(triples)} 3-cliques of mutually anticommuting gammas")

sig_groups = {}
for t in triples:
    signs = [pool[i][2] for i in t]
    p = signs.count(+1)
    q = signs.count(-1)
    sig_groups.setdefault((p, q), []).append(t)

for sig in sorted(sig_groups.keys()):
    cliques = sig_groups[sig]
    print(f"\nSignature ({sig[0]}, {sig[1]}): {len(cliques)} cliques")
    # Compute (G1 G2 G3)^2 for the first clique in each signature
    c = cliques[0]
    G1, G2, G3 = pool[c[0]][1], pool[c[1]][1], pool[c[2]][1]
    G5 = G1 * G2 * G3
    G5_sq = simplify(G5 * G5)
    if is_zero(G5_sq):
        scalar = 0
    else:
        scalar = sp.simplify(G5_sq[0, 0])
        if not is_zero(G5_sq - scalar * I4):
            scalar = 'non-scalar'
    names = [pool[i][0] for i in c]
    print(f"  example: {names}")
    print(f"  (G5)^2 = {scalar} * I_4 = {sp.simplify(scalar)} * I_4")
