"""
Paper 203 v0.2 algebraic verification — Mr A's three deepest catches.

Sub 1: §5 hinge — is i_QM identifiable with Gamma_adj structurally?
Sub 2: §4/§6/§8 time tension — is iteration of Gamma_seed = fourth-gamma orbit?
Sub 3: §8 Constraint 3 — does (Gamma_5)^2 = -5^(3/2)·I uniquely select (1,3)?

All exact symbolic Sympy over Z[phi].
"""

from __future__ import annotations

import sympy as sp
from sympy import sqrt, Matrix, eye, zeros, kronecker_product, simplify, Rational, Symbol, I as sp_I

# ---------------------------------------------------------------------------
# Foundation (matches golden-dirac/golden_dirac.py)
# ---------------------------------------------------------------------------

phi = (1 + sqrt(5)) / 2
psi = 1 - phi
root5 = 2*phi - 1            # = sqrt(5) ∈ Z[phi]

I2 = eye(2)
I4 = eye(4)
Z2 = zeros(2)
Z4 = zeros(4)

# 2x2 building blocks
K = Matrix([[0, 1], [1, 0]])    # K^2 = +I, swap (sigma_x like)
J = Matrix([[0, -1], [1, 0]])   # J^2 = -I (sigma_y like times i)
Z = K * J                        # = [[1, 0], [0, -1]]  Z^2 = +I

Gamma_seed = Matrix([[0, 1], [root5, 0]])     # Gamma_seed^2 = +sqrt(5) I
Gamma_adj  = Matrix([[0, -1], [root5, 0]])    # Gamma_adj^2  = -sqrt(5) I

# 4x4 spatial gammas (Paper 191)
G1 = kronecker_product(K, Gamma_seed)         # +sqrt(5) I_4
G2 = kronecker_product(J, Gamma_seed)         # -sqrt(5) I_4
G3 = kronecker_product(I2, Gamma_adj)         # -sqrt(5) I_4


def anticomm(A, B):
    return sp.simplify(A * B + B * A)


def comm(A, B):
    return sp.simplify(A * B - B * A)


def is_zero(M):
    return simplify(M) == sp.zeros(*M.shape)


def is_scalar_times_I(M, scalar):
    """Check M = scalar * I_n."""
    n = M.shape[0]
    return is_zero(M - scalar * eye(n))


def matrix_subs_root5_flip(M):
    """Apply Galois conjugation sigma: sqrt(5) -> -sqrt(5) entrywise."""
    return M.applyfunc(lambda x: sp.expand(x.subs(sqrt(5), -sqrt(5))))


# ---------------------------------------------------------------------------
# Sub 1: §5 hinge
# ---------------------------------------------------------------------------

def sub1_a():
    """i_golden = Gamma_adj / 5^(1/4). Verify (i_golden)^2 = -I_2."""
    fifth_root = root5 ** sp.Rational(1, 2)   # 5^(1/4)
    i_golden = Gamma_adj / fifth_root
    sq = simplify(i_golden * i_golden)
    target = -I2
    pass_ = is_zero(sq - target)
    return dict(
        i_golden=i_golden,
        square=sq,
        passed=pass_,
    )


def sub1_b():
    """
    Galois conjugation sigma: sqrt(5) -> -sqrt(5). Apply to Gamma_adj.
    Check whether sigma(Gamma_adj) = -Gamma_adj.
    Also check what relationship actually holds.
    """
    sigma_Gamma_adj = matrix_subs_root5_flip(Gamma_adj)
    minus_Gamma_adj = -Gamma_adj
    literal_pass = is_zero(sigma_Gamma_adj - minus_Gamma_adj)

    # Also check alternative relationships
    minus_Gamma_seed = -Gamma_seed
    matches_minus_seed = is_zero(sigma_Gamma_adj - minus_Gamma_seed)

    plus_Gamma_seed = Gamma_seed
    matches_plus_seed = is_zero(sigma_Gamma_adj - plus_Gamma_seed)

    # Transposed
    minus_Gamma_adj_T = -Gamma_adj.T
    matches_minus_adj_T = is_zero(sigma_Gamma_adj - minus_Gamma_adj_T)

    return dict(
        sigma_Gamma_adj=sigma_Gamma_adj,
        literal_pass=literal_pass,
        matches_minus_seed=matches_minus_seed,
        matches_plus_seed=matches_plus_seed,
        matches_minus_adj_T=matches_minus_adj_T,
    )


def sub1_c():
    """
    R(theta) = exp(theta * Gamma_adj).
    Using Gamma_adj^2 = -sqrt(5) I:
      R(theta) = cos(5^(1/4) theta) I + sin(5^(1/4) theta) (Gamma_adj / 5^(1/4))
    Find conserved quadratic form M: Gamma_adj^T M + M Gamma_adj = 0.
    """
    # Find M such that Gamma_adj^T * M + M * Gamma_adj = 0
    a, b, c, d = sp.symbols('a b c d', real=True)
    M = Matrix([[a, b], [c, d]])
    cond = Gamma_adj.T * M + M * Gamma_adj
    cond_simplified = simplify(cond)
    eqs = [cond_simplified[i, j] for i in range(2) for j in range(2)]
    sol = sp.solve(eqs, (a, b, c, d), dict=True)

    # Filter trivial (all zero)
    nontrivial = []
    for s in sol:
        if not all(simplify(v) == 0 for v in s.values()):
            nontrivial.append(s)

    M_solutions = []
    for s in sol:
        M_sol = M.subs(s)
        M_solutions.append(M_sol)

    # Verify group composition R(t1) R(t2) = R(t1+t2) using known closed form
    t1, t2 = sp.symbols('t1 t2', real=True)
    fifth_root = root5 ** sp.Rational(1, 2)
    R1 = sp.cos(fifth_root * t1) * I2 + sp.sin(fifth_root * t1) * (Gamma_adj / fifth_root)
    R2 = sp.cos(fifth_root * t2) * I2 + sp.sin(fifth_root * t2) * (Gamma_adj / fifth_root)
    R12 = simplify(R1 * R2)
    R_sum = sp.cos(fifth_root * (t1+t2)) * I2 + sp.sin(fifth_root * (t1+t2)) * (Gamma_adj / fifth_root)
    diff = simplify(R12 - R_sum)
    group_passes = is_zero(diff)

    return dict(
        M_solutions=M_solutions,
        nontrivial_M=len(nontrivial),
        group_composition_passes=group_passes,
    )


def sub1_d():
    """
    Canonical commutator analogue.
    Search building-block basis for X (Hermitian over Q(sqrt 5)) and P
    such that [X, P] is proportional to Gamma_adj (or i_golden).
    """
    fifth_root = root5 ** sp.Rational(1, 2)
    i_golden = Gamma_adj / fifth_root

    # 2x2 candidates
    basis = {
        'I': I2, 'K': K, 'J': J, 'Z': Z,
        'Gseed': Gamma_seed, 'Gadj': Gamma_adj,
    }

    results = []
    for name_X, X in basis.items():
        # X^dagger = X.T (real entries; conjugate-transpose on Q(sqrt 5))
        if not is_zero(X - X.T):
            continue  # X not Hermitian-symmetric
        for name_P, P in basis.items():
            comm_XP = comm(X, P)
            if is_zero(comm_XP):
                continue
            # Is [X, P] proportional to Gamma_adj?
            # Try ratios: assume [X, P] = c * Gamma_adj for scalar c
            # Solve: comm_XP[0,1] = c * Gamma_adj[0,1] = c * (-1)
            #        comm_XP[1,0] = c * Gamma_adj[1,0] = c * sqrt(5)
            if Gamma_adj[0, 1] != 0:
                c_candidate = comm_XP[0, 1] / Gamma_adj[0, 1]
            elif Gamma_adj[1, 0] != 0:
                c_candidate = comm_XP[1, 0] / Gamma_adj[1, 0]
            else:
                continue
            diff = simplify(comm_XP - c_candidate * Gamma_adj)
            if is_zero(diff):
                results.append(dict(
                    X=name_X, P=name_P, c=simplify(c_candidate),
                    commutator=comm_XP,
                ))

    return dict(solutions=results, n_solutions=len(results))


def sub1_e(H_label="K"):
    """
    Schrodinger analogue. Hermitian H; evolve U = exp(-i_golden * H * dt).
    Test whether U preserves the golden inner product M (from sub1_c).
    """
    fifth_root = root5 ** sp.Rational(1, 2)
    i_golden = Gamma_adj / fifth_root

    H_map = {'K': K, 'Z': Z, 'I': I2}
    H = H_map[H_label]
    if not is_zero(H - H.T):
        return dict(error="H not Hermitian")

    # Try the closed form: since (i_golden)^2 = -I,
    # exp(-i_golden * H * dt) = cos(H dt) I - sin(H dt) i_golden ... but H, i_golden may not commute
    # For diagonal H = Z (eigenvalues +-1), check:
    # If [H, i_golden] = 0: clean exponential. Otherwise: Baker-Campbell.

    commute = is_zero(comm(H, i_golden))

    # Without commuting, we cannot simply factorize. Check norm preservation
    # numerically-symbolically for small dt via Taylor:
    dt = sp.Symbol('dt', real=True)
    A = -i_golden * H * dt  # 2x2
    # Truncated exponential to order 4:
    U4 = I2 + A + A*A/2 + A*A*A/6 + A*A*A*A/24

    # Conserved form M from sub1_c: pick one of the M solutions
    sub1c = sub1_c()
    if not sub1c['M_solutions']:
        return dict(error="no M to test")
    # The first solution: pick the one with non-trivial parameterisation.
    # In sub1_c, M = [[a, b], [c, d]] with Gamma_adj^T M + M Gamma_adj = 0.
    # Solve and choose one (we'll do this fresh here).
    a, b, c_var, d = sp.symbols('a b c d', real=True)
    M_param = Matrix([[a, b], [c_var, d]])
    cond = simplify(Gamma_adj.T * M_param + M_param * Gamma_adj)
    eqs = [cond[i, j] for i in range(2) for j in range(2)]
    M_sols = sp.solve(eqs, (a, b, c_var, d), dict=True)
    if not M_sols:
        return dict(error="no M solutions")
    M_sol = M_param.subs(M_sols[0])
    # Pick a specific non-zero parametrisation
    M_specific = M_sol.subs(list(M_sol.free_symbols)[0], 1) if M_sol.free_symbols else M_sol
    # If still has free symbols, pick another
    for s in list(M_specific.free_symbols):
        M_specific = M_specific.subs(s, 0)
    # Check norm preservation at order dt^2:
    norm_diff = simplify(U4.T * M_specific * U4 - M_specific)
    # Coefficient of dt^0: should be 0 (M_specific identity check)
    # Coefficient of dt^1: should be 0
    # Coefficient of dt^2: should be 0 if M is conserved

    # Easier: check Lie-algebra condition directly: Gamma_adj^T M + M Gamma_adj = 0 (already known)
    # AND H^T M = M H (the condition for H to preserve M)
    H_preserves_M = is_zero(H.T * M_specific - M_specific * H)

    return dict(
        H=H_label,
        H_commutes_with_i_golden=commute,
        M_specific=M_specific,
        H_preserves_M=H_preserves_M,
    )


def sub1_null():
    """
    Test alternative i-candidates (must square to -I_2).
    """
    fifth_root = root5 ** sp.Rational(1, 2)
    candidates = {
        'J (Pauli i)': J,
        'i_golden = Gamma_adj / 5^(1/4)': Gamma_adj / fifth_root,
        'Gamma_seed (NOT a candidate; sq=+sqrt5)': Gamma_seed,
    }
    out = {}
    for name, M in candidates.items():
        sq = simplify(M * M)
        is_minus_I = is_zero(sq + I2)
        out[name] = dict(square=sq, is_minus_I=is_minus_I)
    return out


# ---------------------------------------------------------------------------
# Sub 2: fourth gamma
# ---------------------------------------------------------------------------

def sub2_a():
    """
    Enumerate Gamma_0 = A x B candidates with A, B in building blocks.
    Check anticommutation with G1, G2, G3 and (Gamma_0)^2.
    """
    basis = {
        'I': I2, 'K': K, 'J': J, 'Z': Z,
        'Gseed': Gamma_seed, 'Gadj': Gamma_adj,
    }

    valid_candidates = []
    for nA, A in basis.items():
        for nB, B in basis.items():
            G0 = kronecker_product(A, B)
            sq = simplify(G0 * G0)
            # Categorise the square
            sq_class = None
            if is_zero(sq - root5 * I4):
                sq_class = '+sqrt(5)*I'
            elif is_zero(sq + root5 * I4):
                sq_class = '-sqrt(5)*I'
            elif is_zero(sq - I4):
                sq_class = '+I'
            elif is_zero(sq + I4):
                sq_class = '-I'
            else:
                sq_class = 'other'

            ac1 = is_zero(anticomm(G0, G1))
            ac2 = is_zero(anticomm(G0, G2))
            ac3 = is_zero(anticomm(G0, G3))

            if ac1 and ac2 and ac3:
                valid_candidates.append(dict(
                    A=nA, B=nB,
                    sq_class=sq_class,
                    sq=sq,
                ))
    return dict(candidates=valid_candidates)


def sub2_b():
    """Iteration tower Gamma_seed^n for n=1..6; check norm cascade."""
    rows = []
    M = I2
    for n in range(1, 7):
        M = simplify(M * Gamma_seed)
        # Operator norm via singular values: largest sqrt(eigenvalue of M^T M)
        MtM = simplify(M.T * M)
        evals = MtM.eigenvals()
        max_eig = sp.simplify(sp.Max(*[sp.sqrt(e) for e in evals.keys()]))
        # Frobenius norm
        frob = sp.sqrt(simplify(sum(M[i, j]**2 for i in range(2) for j in range(2))))
        rows.append(dict(n=n, M=M, op_norm=max_eig, frob_norm=simplify(frob)))
    return rows


def sub2_c(G0):
    """Compute exp(t * G0) symbolically, given G0^2 = +sqrt(5) I."""
    t = sp.Symbol('t', real=True)
    fifth_root = root5 ** sp.Rational(1, 2)
    # exp(t G0) = cosh(5^(1/4) t) I + sinh(5^(1/4) t) (G0 / 5^(1/4))
    closed_form = sp.cosh(fifth_root * t) * I4 + sp.sinh(fifth_root * t) * (G0 / fifth_root)
    # Sanity: Taylor first three terms
    # I + t G0 + (t G0)^2/2 + ...
    A = t * G0
    taylor = I4 + A + (A*A) / 2 + (A*A*A)/6
    # Reduce both to same form
    diff = simplify(closed_form.series(t, 0, 4).removeO() - taylor) if False else None
    return dict(closed_form=closed_form)


def sub2_d(G0):
    """Identification test: relate Gamma_seed^n iteration to exp(t G0)."""
    # Gamma_seed^n acts on C^2; G0 acts on C^4. They're in different algebras.
    # Question: is there an embedding of Gamma_seed into 4x4 that matches G0?
    # The natural embedding: G1 = K x Gamma_seed has Gamma_seed in it.
    # G1^2 = +sqrt(5) I_4, same as G0. So G1 and G0 are both "timelike" gammas.
    # Are G0 and G1 the same up to similarity?
    rows_2b = sub2_b()
    # Comparison: G0 vs G1
    G1_sq = simplify(G1 * G1)
    G0_sq = simplify(G0 * G0)
    same_square = is_zero(G1_sq - G0_sq)
    # Do G0 and G1 commute or anticommute?
    G0_G1_comm = is_zero(comm(G0, G1))
    G0_G1_anticomm = is_zero(anticomm(G0, G1))
    return dict(
        G1_sq_equals_G0_sq=same_square,
        G0_G1_commute=G0_G1_comm,
        G0_G1_anticommute=G0_G1_anticomm,
    )


# ---------------------------------------------------------------------------
# Sub 3: signature enumeration
# ---------------------------------------------------------------------------

def sub3_enumerate_4tuples():
    """
    Find all 4-tuples of mutually anticommuting 4x4 gammas A x B,
    grouped by signature (p, q).
    """
    basis = {
        'I': I2, 'K': K, 'J': J, 'Z': Z,
        'Gseed': Gamma_seed, 'Gadj': Gamma_adj,
    }
    # Build the full pool of candidate gammas (A x B for each basis pair),
    # keep only those with squares = +sqrt(5)*I or -sqrt(5)*I
    pool = []
    for nA, A in basis.items():
        for nB, B in basis.items():
            G = kronecker_product(A, B)
            sq = simplify(G * G)
            if is_zero(sq - root5 * I4):
                pool.append(dict(name=f"{nA}x{nB}", G=G, sign=+1))
            elif is_zero(sq + root5 * I4):
                pool.append(dict(name=f"{nA}x{nB}", G=G, sign=-1))

    # Quick check: known three from Paper 191
    print(f"  Pool size (gammas with sqrt(5)·I squares): {len(pool)}")

    # Now find 4-tuples that pairwise anticommute
    # For efficiency, precompute anticommutation graph
    n_pool = len(pool)
    ac_graph = [[False] * n_pool for _ in range(n_pool)]
    for i in range(n_pool):
        for j in range(i+1, n_pool):
            if is_zero(anticomm(pool[i]['G'], pool[j]['G'])):
                ac_graph[i][j] = True
                ac_graph[j][i] = True

    # Find all 4-cliques in this graph
    import itertools
    cliques = []
    indices = range(n_pool)
    for combo in itertools.combinations(indices, 4):
        ok = True
        for a, b in itertools.combinations(combo, 2):
            if not ac_graph[a][b]:
                ok = False
                break
        if ok:
            cliques.append(combo)

    # Group by signature (p, q) where p = number of +1 squares
    sig_groups: dict[tuple, list] = {}
    for clique in cliques:
        signs = [pool[i]['sign'] for i in clique]
        p = signs.count(+1)
        q = signs.count(-1)
        sig = (p, q)
        sig_groups.setdefault(sig, []).append(clique)

    return dict(pool=pool, cliques=cliques, sig_groups=sig_groups)


def sub3_chirality_per_signature(enum_result):
    """For each signature with at least one clique, compute (Gamma^5)^2."""
    pool = enum_result['pool']
    sig_groups = enum_result['sig_groups']
    out = {}
    for sig, cliques in sig_groups.items():
        # Take the first clique; compute Gamma_0 * Gamma_1 * Gamma_2 * Gamma_3
        # and its square
        c = cliques[0]
        gammas = [pool[i]['G'] for i in c]
        G5 = gammas[0] * gammas[1] * gammas[2] * gammas[3]
        G5_sq = simplify(G5 * G5)
        # Try to express as scalar * I_4
        if is_zero(G5_sq):
            scalar = 0
        else:
            # extract scalar from diagonal
            scalar = sp.simplify(G5_sq[0, 0])
            # verify scalar * I
            if not is_zero(G5_sq - scalar * I4):
                scalar = 'non-scalar'
        out[sig] = dict(
            n_cliques=len(cliques),
            example_names=[pool[i]['name'] for i in c],
            G5_sq_scalar=scalar,
        )
    return out


if __name__ == '__main__':
    print("=" * 60)
    print("Paper 203 v0.2 algebraic verification")
    print("=" * 60)

    # ----- Sub 1 -----
    print("\n=== Sub 1 — §5 hinge ===")

    print("\n(1.a) (i_golden)^2 = -I_2?")
    r = sub1_a()
    print(f"  PASS = {r['passed']}")
    print(f"  square = {r['square'].tolist()}")

    print("\n(1.b) sigma(Gamma_adj) = -Gamma_adj?")
    r = sub1_b()
    print(f"  sigma(Gamma_adj) = {r['sigma_Gamma_adj'].tolist()}")
    print(f"  literal pass (= -Gamma_adj): {r['literal_pass']}")
    print(f"  alternative: = -Gamma_seed: {r['matches_minus_seed']}")
    print(f"  alternative: = +Gamma_seed: {r['matches_plus_seed']}")
    print(f"  alternative: = -Gamma_adj^T: {r['matches_minus_adj_T']}")

    print("\n(1.c) R(theta) = exp(theta Gamma_adj): group composition and conserved form")
    r = sub1_c()
    print(f"  group composition R(t1)R(t2) = R(t1+t2): {r['group_composition_passes']}")
    print(f"  number of M solutions to Gamma_adj^T M + M Gamma_adj = 0: {len(r['M_solutions'])}")
    for i, M in enumerate(r['M_solutions']):
        print(f"  M[{i}] = {M.tolist()}")

    print("\n(1.d) Canonical commutator [X, P] = c Gamma_adj")
    r = sub1_d()
    print(f"  Number of non-trivial solutions found: {r['n_solutions']}")
    for s in r['solutions'][:10]:
        print(f"  [{s['X']}, {s['P']}] = ({s['c']}) Gamma_adj")

    print("\n(1.e) Schrodinger analogue")
    for H_lbl in ('K', 'Z', 'I'):
        r = sub1_e(H_lbl)
        if 'error' in r:
            print(f"  H={H_lbl}: {r['error']}")
        else:
            print(f"  H={H_lbl}: [H, i_golden] = 0: {r['H_commutes_with_i_golden']}, H preserves M: {r['H_preserves_M']}")

    print("\n(Sub 1 null) Candidates for i (must square to -I_2)")
    r = sub1_null()
    for name, info in r.items():
        print(f"  {name}: square = {info['square'].tolist()}, is -I_2 = {info['is_minus_I']}")

    # ----- Sub 2 -----
    print("\n\n=== Sub 2 — Time tension ===")

    print("\n(2.a) Fourth-gamma candidates")
    r = sub2_a()
    print(f"  Total Gamma_0 candidates anticommuting with G1, G2, G3: {len(r['candidates'])}")
    for c in r['candidates']:
        print(f"    {c['A']} x {c['B']}: G0^2 = {c['sq_class']}")

    # Find the canonical +sqrt(5)*I candidate
    positive_candidates = [c for c in r['candidates'] if c['sq_class'] == '+sqrt(5)*I']
    print(f"\n  +sqrt(5)*I candidates: {len(positive_candidates)}")
    if positive_candidates:
        first = positive_candidates[0]
        # Reconstruct G0
        basis_map = {'I': I2, 'K': K, 'J': J, 'Z': Z, 'Gseed': Gamma_seed, 'Gadj': Gamma_adj}
        G0_canonical = kronecker_product(basis_map[first['A']], basis_map[first['B']])
        print(f"  Canonical G0 = {first['A']} x {first['B']}")

        print("\n(2.b) Iteration tower of Gamma_seed")
        tower = sub2_b()
        for row in tower:
            print(f"    n={row['n']}: op_norm = {row['op_norm']}, frob = {row['frob_norm']}")

        print("\n(2.c) exp(t G0): closed form computed")
        r = sub2_c(G0_canonical)
        # don't print giant matrix; show shape
        print(f"    Closed form has shape {r['closed_form'].shape}")

        print("\n(2.d) Identification test")
        r = sub2_d(G0_canonical)
        print(f"    G0^2 = G1^2: {r['G1_sq_equals_G0_sq']}")
        print(f"    [G0, G1] = 0: {r['G0_G1_commute']}")
        print(f"    {{G0, G1}} = 0: {r['G0_G1_anticommute']}")

    # ----- Sub 3 -----
    print("\n\n=== Sub 3 — Signature enumeration ===")
    print("\n(3.a) Enumerating 4-tuples of mutually anticommuting gammas...")
    enum = sub3_enumerate_4tuples()
    print(f"  Total 4-cliques found: {len(enum['cliques'])}")
    print(f"  Signature distribution:")
    for sig in sorted(enum['sig_groups'].keys()):
        n = len(enum['sig_groups'][sig])
        print(f"    ({sig[0]}, {sig[1]}): {n} cliques")

    print("\n(3.b)/(3.c) (Gamma^5)^2 per signature:")
    chir = sub3_chirality_per_signature(enum)
    for sig in sorted(chir.keys()):
        info = chir[sig]
        print(f"  ({sig[0]}, {sig[1]}): {info['n_cliques']} cliques; "
              f"example {info['example_names']}; (G5)^2 = {info['G5_sq_scalar']} * I_4")
