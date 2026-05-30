"""
Precompute L-function zeros for d in {5, 13, 3, 7} using the fast
truncated-Dirichlet method (validated against mpmath chi_5 to max delta
0.02, well under the W=1.0 matching window).
"""

import time
import argparse
from l_fast import get_or_compute_zeros_fast


# Target: ~1000 zeros per character to match Paper 150 v2.0's N=1000 spec.
T_MAX_PER_D = {
    5:  1000.0,   # ~ 900 zeros
    13: 600.0,    # ~ 800 zeros
    3:  700.0,    # ~ 700 zeros (cond 12)
    7:  500.0,    # ~ 700 zeros (cond 28)
}

M_TRUNC = 1000   # truncation; validated chi_5 to delta <= 0.02
DT_SCAN = 0.05


def main(force=False):
    for d, tmax in T_MAX_PER_D.items():
        cond = {3: 12, 5: 5, 7: 28, 13: 13}[d]
        print(f"\n=== d = {d}, conductor = {cond} ===")
        t0 = time.time()
        zeros = get_or_compute_zeros_fast(d, tmax, dt=DT_SCAN, M=M_TRUNC,
                                          force=force)
        print(f"  {len(zeros)} zeros up to t = {tmax} ({time.time()-t0:.1f}s)")
        if zeros:
            print(f"  first 3: {[f'{z:.3f}' for z in zeros[:3]]}")
            print(f"  last 3:  {[f'{z:.3f}' for z in zeros[-3:]]}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    main(force=args.force)
