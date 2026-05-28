"""
Continuous S^3 quaternion dynamics with lag-coupled icosahedral drift.

Each unit i carries q_i in S^3 (unit quaternion in R^4). At each step:

    q_drift  = c · q_i(t)                         # 2I drift (5-fold rotation)
    q_neigh  = mean_{j != i} K_ij · q_j(t - tau)  # lagged neighbour mean in R^4
    q_i(t+1) = normalize( (1 - eps) * q_drift + eps * q_neigh )

Per PRE_REGISTRATION.md.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

PHI = (1.0 + np.sqrt(5.0)) / 2.0

# 2I drift element: 5-fold rotation, quaternion (phi/2, 1/(2phi), 1/2, 0).
# This is a vertex of the 600-cell; represents rotation by 2pi/5 about the
# axis (1/(2 phi), 1/2, 0) / ||·||.
C_2I_DRIFT = np.array([PHI / 2.0, 1.0 / (2.0 * PHI), 0.5, 0.0], dtype=np.float64)
assert abs(np.linalg.norm(C_2I_DRIFT) - 1.0) < 1e-12, "c must be unit norm"


def quat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton quaternion product (w, x, y, z) · (w', x', y', z')."""
    aw, ax, ay, az = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    bw, bx, by, bz = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ], axis=-1)


def normalize_q(q: np.ndarray, eps_norm: float = 1e-12) -> np.ndarray:
    n = np.linalg.norm(q, axis=-1, keepdims=True)
    return q / np.maximum(n, eps_norm)


def random_unit_quaternions(N: int, rng: np.random.Generator) -> np.ndarray:
    q = rng.standard_normal(size=(N, 4))
    return normalize_q(q)


@dataclass
class RunConfig:
    N: int = 240
    T: int = 8000
    transient: int = 2000
    eps: float = 0.1
    tau: int = 5
    drift_q: np.ndarray | None = None  # if None: use C_2I_DRIFT
    K_mode: str = "all-to-all"         # "all-to-all" | "shuffled"
    K_matrix: np.ndarray | None = None  # only used when K_mode="shuffled"
    seed: int = 20260528


def run_dynamics(cfg: RunConfig) -> dict:
    """Run the dynamics. Returns dict with keys 'states', 'corr', 'meta'."""
    rng = np.random.default_rng(cfg.seed)
    drift = cfg.drift_q if cfg.drift_q is not None else C_2I_DRIFT
    drift = np.asarray(drift, dtype=np.float64)
    assert drift.shape == (4,)
    assert abs(np.linalg.norm(drift) - 1.0) < 1e-9

    # Initial states: uniform on S^3
    q = random_unit_quaternions(cfg.N, rng)

    # Ring buffer for lag
    buf_len = max(cfg.tau + 1, 1)
    buf = np.tile(q, (buf_len, 1, 1))  # shape (buf_len, N, 4)
    buf_idx = 0

    # Coupling matrix
    if cfg.K_mode == "all-to-all":
        K = np.full((cfg.N, cfg.N), 1.0 / max(cfg.N - 1, 1), dtype=np.float64)
        np.fill_diagonal(K, 0.0)
    elif cfg.K_mode == "shuffled":
        assert cfg.K_matrix is not None
        K = cfg.K_matrix.copy()
        # Row-normalise so each unit's neighbour mean is well-defined
        rowsum = K.sum(axis=1, keepdims=True)
        K = np.where(rowsum > 0, K / np.maximum(rowsum, 1e-12), K)
    else:
        raise ValueError(cfg.K_mode)

    # Storage for the recording window
    record_T = cfg.T - cfg.transient
    record = np.zeros((record_T, cfg.N, 4), dtype=np.float64)

    # Time loop
    for t in range(cfg.T):
        # lagged states from ring buffer
        lag_idx = (buf_idx - cfg.tau) % buf_len
        q_lag = buf[lag_idx]                          # (N, 4)
        q_neigh = K @ q_lag                           # (N, 4)
        q_drift = quat_mul(np.broadcast_to(drift, (cfg.N, 4)), q)
        q_new = (1.0 - cfg.eps) * q_drift + cfg.eps * q_neigh
        q = normalize_q(q_new)

        # advance buffer
        buf_idx = (buf_idx + 1) % buf_len
        buf[buf_idx] = q

        if t >= cfg.transient:
            record[t - cfg.transient] = q

    # Correlation matrix C_ij = <q_i(t) · q_j(t)>_t
    # Use stable formulation: C = (1/T) sum_t q(t).T @ q(t) - but we want
    # inner products per pair. Vectorised:
    #   C_ij = (1/T) sum_t sum_d q[t,i,d] * q[t,j,d]
    # This equals (1/T) (Q.T @ Q) where Q is (T, N, 4) reshaped to (T*4, N) ...
    # easier: just average per-step inner products.
    C = np.zeros((cfg.N, cfg.N), dtype=np.float64)
    for t in range(record_T):
        C += record[t] @ record[t].T
    C /= record_T
    # numerical clip
    C = np.clip(C, -1.0, 1.0)

    # Sanity: variance of correlation distribution
    iu = np.triu_indices(cfg.N, k=1)
    off = C[iu]
    meta = dict(
        N=cfg.N, T=cfg.T, transient=cfg.transient,
        eps=cfg.eps, tau=cfg.tau, K_mode=cfg.K_mode,
        mean_corr=float(np.mean(off)),
        std_corr=float(np.std(off)),
        min_corr=float(np.min(off)),
        max_corr=float(np.max(off)),
        drift_q=drift.tolist(),
    )
    return dict(record=record, corr=C, meta=meta)


def correlation_to_distance(C: np.ndarray) -> np.ndarray:
    """d_ij = arccos(C_ij). Symmetric in [0, pi]."""
    return np.arccos(np.clip(C, -1.0, 1.0))


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------

def weighted_graph_from_corr(C: np.ndarray, drop_negative: bool = True):
    """Return W (similarity weights) and L (graph Laplacian) as numpy arrays."""
    W = C.copy()
    np.fill_diagonal(W, 0.0)
    if drop_negative:
        W = np.maximum(W, 0.0)
    D = np.diag(W.sum(axis=1))
    L = D - W
    return W, L


def spectral_dimension(L: np.ndarray, t_range: tuple[float, float] = (0.5, 20.0),
                       n_t: int = 30) -> tuple[float, np.ndarray, np.ndarray]:
    """
    Fit p(t) ~ t^{-d_s/2} from heat-kernel return probability.
    """
    n = L.shape[0]
    # Use all eigenvalues (dense N up to 240 is fine).
    evals = np.linalg.eigvalsh(L)
    evals = np.maximum(evals, 0.0)
    ts = np.geomspace(t_range[0], t_range[1], n_t)
    p = np.zeros_like(ts)
    for i, t in enumerate(ts):
        p[i] = np.mean(np.exp(-t * evals))
    # subtract floor 1/n
    floor = 1.0 / n
    p_signal = p - floor
    valid = p_signal > 1e-12
    if valid.sum() < 5:
        return float('nan'), ts, p
    # fit middle 60% of valid range
    valid_idx = np.where(valid)[0]
    lo = int(0.2 * len(valid_idx))
    hi = int(0.8 * len(valid_idx))
    sl = valid_idx[lo:hi] if hi - lo >= 5 else valid_idx
    slope, _ = np.polyfit(np.log(ts[sl]), np.log(p_signal[sl]), 1)
    d_s = -2.0 * slope
    return float(d_s), ts, p


def volume_growth(D_mat: np.ndarray, n_r: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """
    For each unit u, count N_u(r) = number of units j with d_uj <= r.
    Average over u. Return (r_grid, N_mean_of_r).
    """
    n = D_mat.shape[0]
    # symmetric distances; consider only off-diagonal
    D = D_mat.copy()
    np.fill_diagonal(D, 0.0)
    rmax = float(np.max(D))
    r_grid = np.linspace(0.0, rmax, n_r)
    Nr = np.zeros(n_r)
    for k, r in enumerate(r_grid):
        Nr[k] = np.mean(np.sum(D <= r, axis=1))
    return r_grid, Nr


def fit_volume_dim_and_saturation(r_grid: np.ndarray, Nr: np.ndarray,
                                   N_total: int) -> dict:
    """
    Fit small-r power-law slope; check saturation (Nr -> N_total at large r).
    Return dict.
    """
    # avoid r=0
    valid = (r_grid > 0) & (Nr > 1)
    rv = r_grid[valid]
    Nv = Nr[valid]
    if len(rv) < 5:
        return dict(slope=float('nan'), saturated=False, sat_frac=float('nan'))
    # use first 40% of range for slope
    cut = max(int(0.4 * len(rv)), 5)
    slope, intercept = np.polyfit(np.log(rv[:cut]), np.log(Nv[:cut]), 1)
    # saturation: ratio of final Nr to N_total
    sat_frac = float(Nv[-1]) / N_total
    saturated = sat_frac > 0.95  # 95% of units within max range
    return dict(slope=float(slope), saturated=bool(saturated), sat_frac=sat_frac,
                Nr_final=float(Nv[-1]))


def mean_triangle_defect(D_mat: np.ndarray, n_samples: int = 5000,
                         seed: int = 0) -> dict:
    """
    Sample triples (i, j, k); compute defect = d_ij + d_jk + d_ki − π.
    For a spherical triangle on S^3 (positive curvature), defect can be
    positive or negative depending on triangle size; we report the mean
    and the fraction positive.
    """
    n = D_mat.shape[0]
    rng = np.random.default_rng(seed)
    defects = []
    sums = []
    for _ in range(n_samples):
        idx = rng.choice(n, size=3, replace=False)
        i, j, k = idx
        s = D_mat[i, j] + D_mat[j, k] + D_mat[k, i]
        sums.append(s)
        defects.append(s - np.pi)
    defects = np.array(defects)
    sums = np.array(sums)
    return dict(
        mean_defect=float(np.mean(defects)),
        median_defect=float(np.median(defects)),
        frac_positive=float(np.mean(defects > 0)),
        mean_perimeter=float(np.mean(sums)),
    )
