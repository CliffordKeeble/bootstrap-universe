"""
Drive the full N=12, 60, 120 sweep for all three constructions, plus the
600-cell reference, with null at each N. Write observables CSV.

Per stop-on-fail at N=120 the run halts at N=120; no N=600/1200/6000 runs
are performed for face/axis (primary protocol triggered).
"""

import csv
import time
import numpy as np
import networkx as nx
from build_network import (
    build_face_shared, build_axis_shared, build_vertex_shared,
    network_to_graph, analyse_graph, null_distribution, SEED,
)
from reference_600cell import build_600cell_graph

N_SIZES = [12, 60, 120]


def main():
    rows = []
    for cname, builder in [('face', build_face_shared),
                            ('axis', build_axis_shared),
                            ('vertex', build_vertex_shared)]:
        t0 = time.time()
        icos, pool, edges = builder(max(N_SIZES))
        G_full = network_to_graph(icos, edges)
        t_build = time.time() - t0
        print(f"\n=== {cname} ===")
        print(f"  built in {t_build:.2f}s: {len(icos)} icosahedra, {len(pool.positions)} vertices, "
              f"{G_full.number_of_edges()} edges")
        for N in N_SIZES:
            if N > len(icos):
                continue
            Gn = G_full.subgraph(range(N)).copy()
            Gn = nx.convert_node_labels_to_integers(Gn)
            res = analyse_graph(Gn, N, want_heat=False)
            lam_null, d_null = null_distribution(N, res.avg_degree, n_samples=20)
            z_lam = (res.lambda1_raw - lam_null.mean()) / (lam_null.std() + 1e-12)
            row = {
                'construction': cname,
                'N': N,
                'edges': res.n_edges,
                'avg_degree': res.avg_degree,
                'n_components': res.n_components,
                'lambda1_raw': res.lambda1_raw,
                'lambda1_per_deg': res.lambda1_per_deg,
                'lambda1_x_N23': res.lambda1_x_N23,
                'dim_d': res.dim_d,
                'dim_residual': res.dim_residual,
                'null_lambda1_mean': float(lam_null.mean()),
                'null_lambda1_std': float(lam_null.std()),
                'null_d_mean': float(np.nanmean(d_null)),
                'null_d_std': float(np.nanstd(d_null)),
                'z_lambda1': float(z_lam),
            }
            rows.append(row)
            print(f"  N={N}: lam1={res.lambda1_raw:.4f} "
                  f"lam1*N^2/3={res.lambda1_x_N23:.3f} "
                  f"d={res.dim_d:.3f} z={z_lam:.2f}")

    # 600-cell reference at N=120
    print(f"\n=== 600-cell reference ===")
    G600 = build_600cell_graph()
    res = analyse_graph(G600, 120, want_heat=False)
    lam_null, d_null = null_distribution(120, res.avg_degree, n_samples=20)
    z_lam = (res.lambda1_raw - lam_null.mean()) / (lam_null.std() + 1e-12)
    print(f"  N=120 (built exactly): edges={res.n_edges} avg_deg={res.avg_degree:.2f} "
          f"lam1={res.lambda1_raw:.4f} lam1*N^2/3={res.lambda1_x_N23:.3f} d={res.dim_d:.3f} z={z_lam:.2f}")
    rows.append({
        'construction': '600cell',
        'N': 120,
        'edges': res.n_edges,
        'avg_degree': res.avg_degree,
        'n_components': res.n_components,
        'lambda1_raw': res.lambda1_raw,
        'lambda1_per_deg': res.lambda1_per_deg,
        'lambda1_x_N23': res.lambda1_x_N23,
        'dim_d': res.dim_d,
        'dim_residual': res.dim_residual,
        'null_lambda1_mean': float(lam_null.mean()),
        'null_lambda1_std': float(lam_null.std()),
        'null_d_mean': float(np.nanmean(d_null)),
        'null_d_std': float(np.nanstd(d_null)),
        'z_lambda1': float(z_lam),
    })

    # Write CSV
    out_path = 'observables.csv'
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
