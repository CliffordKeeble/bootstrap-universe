"""
v0.2 sweep driver: Scheme A (edge subdivision) and Scheme B (barycentric
cell subdivision), with random-graph null matched on degree sequence.

Stop-on-fail at Scheme A level 1 was identified in the implementation
smoke test (deviation -70.76% > 30% threshold). Per pre-registered
protocol, Scheme A levels 2-4 are NOT executed. Scheme B levels 0-1 are
still executed because Scheme B is a separate refinement scheme.
"""

from __future__ import annotations

import csv
import time
import numpy as np
import networkx as nx

from build_refinement import (
    scheme_A_levels, edge_subdivide, scheme_B_level1, analyse_level,
    null_distribution, smallest_eigs, gap, TARGET_LAMBDA1_X_N23, SEED,
)
from reference_600cell import build_600cell_graph


def fmt_pct(x: float) -> str:
    return f"{x:+.2f}%"


def run_with_null(G: nx.Graph, scheme: str, level: int, n_null: int = 20):
    t0 = time.time()
    res = analyse_level(G, scheme, level)
    t_eig = time.time() - t0
    print(f"  {scheme} lvl{level}: N={res.N} E={res.E} avg_deg={res.avg_degree:.2f} "
          f"comp={res.n_components} lam1={res.lambda1_raw:.6f} "
          f"lam1*N^(2/3)={res.lambda1_x_N23:.3f} dev={fmt_pct(res.deviation_pct)} "
          f"(eig {t_eig:.1f}s)")
    t0 = time.time()
    null = null_distribution(G, n_samples=n_null, seed=SEED + 1000 * level)
    t_null = time.time() - t0
    null_mean = float(np.mean(null))
    null_std = float(np.std(null))
    z = (res.lambda1_x_N23 - null_mean) / (null_std + 1e-12)
    print(f"  {scheme} lvl{level} null: mean={null_mean:.3f} std={null_std:.3f} "
          f"z={z:+.2f} ({t_null:.1f}s)")
    return res, null, null_mean, null_std, z


def main():
    rows = []

    print(f"=== Target lambda_1 * N^(2/3) = {TARGET_LAMBDA1_X_N23:.4f} ===\n")

    # ---- Scheme A ----
    print("=== Scheme A: iterated edge subdivision ===")
    G = build_600cell_graph()
    G = nx.convert_node_labels_to_integers(G)
    # Level 0
    res, null, nm, ns, z = run_with_null(G, 'A', 0)
    rows.append({
        'scheme': 'A', 'level': 0, 'N': res.N, 'E': res.E,
        'avg_degree': res.avg_degree, 'n_components': res.n_components,
        'lambda1_raw': res.lambda1_raw, 'lambda1_x_N23': res.lambda1_x_N23,
        'deviation_pct': res.deviation_pct,
        'null_mean': nm, 'null_std': ns, 'z_score': z,
        'protocol_status': 'reported',
    })
    # Level 1 — pre-registered evaluation point
    G = edge_subdivide(G)
    res, null, nm, ns, z = run_with_null(G, 'A', 1)
    failed_30 = abs(res.deviation_pct) > 30.0
    rows.append({
        'scheme': 'A', 'level': 1, 'N': res.N, 'E': res.E,
        'avg_degree': res.avg_degree, 'n_components': res.n_components,
        'lambda1_raw': res.lambda1_raw, 'lambda1_x_N23': res.lambda1_x_N23,
        'deviation_pct': res.deviation_pct,
        'null_mean': nm, 'null_std': ns, 'z_score': z,
        'protocol_status': 'stop-on-fail-triggered' if failed_30 else 'reported',
    })
    if failed_30:
        print(f"\n*** STOP-ON-FAIL TRIGGERED at Scheme A level 1: "
              f"|dev| = {abs(res.deviation_pct):.2f}% > 30% threshold. "
              f"Scheme A levels 2-4 NOT run per pre-registered protocol. ***\n")
    else:
        # If somehow level 1 passes, continue
        for lvl in range(2, 5):
            G = edge_subdivide(G)
            res, null, nm, ns, z = run_with_null(G, 'A', lvl)
            rows.append({
                'scheme': 'A', 'level': lvl, 'N': res.N, 'E': res.E,
                'avg_degree': res.avg_degree, 'n_components': res.n_components,
                'lambda1_raw': res.lambda1_raw, 'lambda1_x_N23': res.lambda1_x_N23,
                'deviation_pct': res.deviation_pct,
                'null_mean': nm, 'null_std': ns, 'z_score': z,
                'protocol_status': 'reported',
            })
            if abs(res.deviation_pct) > 30.0:
                print(f"\n*** STOP-ON-FAIL at level {lvl}. ***\n")
                break

    # ---- Scheme B ----
    print("\n=== Scheme B: barycentric cell subdivision ===")
    G_B0 = build_600cell_graph()
    G_B0 = nx.convert_node_labels_to_integers(G_B0)
    res, null, nm, ns, z = run_with_null(G_B0, 'B', 0)
    rows.append({
        'scheme': 'B', 'level': 0, 'N': res.N, 'E': res.E,
        'avg_degree': res.avg_degree, 'n_components': res.n_components,
        'lambda1_raw': res.lambda1_raw, 'lambda1_x_N23': res.lambda1_x_N23,
        'deviation_pct': res.deviation_pct,
        'null_mean': nm, 'null_std': ns, 'z_score': z,
        'protocol_status': 'reported',
    })
    G_B1 = scheme_B_level1()
    res, null, nm, ns, z = run_with_null(G_B1, 'B', 1)
    rows.append({
        'scheme': 'B', 'level': 1, 'N': res.N, 'E': res.E,
        'avg_degree': res.avg_degree, 'n_components': res.n_components,
        'lambda1_raw': res.lambda1_raw, 'lambda1_x_N23': res.lambda1_x_N23,
        'deviation_pct': res.deviation_pct,
        'null_mean': nm, 'null_std': ns, 'z_score': z,
        'protocol_status': 'reported',
    })

    # Scheme A vs B agreement at level 0 (common N=120 by construction)
    a0 = next(r for r in rows if r['scheme'] == 'A' and r['level'] == 0)
    b0 = next(r for r in rows if r['scheme'] == 'B' and r['level'] == 0)
    diff = (a0['lambda1_x_N23'] - b0['lambda1_x_N23']) / b0['lambda1_x_N23'] * 100.0
    print(f"\nScheme A vs B at level 0 (N=120, same graph): diff = {fmt_pct(diff)} "
          f"(should be 0%; sanity check)")

    # Write CSV
    out_path = 'observables_v0_2.csv'
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
