# Golden Zeta Investigation

Re-examination of Paper 150 §3.1.1: does ζ_φ's minima-to-Riemann-zero
matching degrade at high t?

Part of the Bootstrap Universe Programme.  
Dr Clifford Keeble, Woodbridge UK  
ORCID: 0009-0003-6828-2155

## How to reproduce

Each script runs standalone with Python 3. Requires `numpy` and `mpmath`.

```bash
pip install numpy mpmath

# Step 1: Baseline reproduction
python golden-zeta/zeta_phi.py

# Steps 2-3: Window sensitivity + null hypothesis
python golden-zeta/null_test.py

# Step 4: Scaling in N
python golden-zeta/scaling_test.py

# Step 5: Per-zero diagnostic (writes per_zero.csv)
python golden-zeta/diagnostic.py
```

## Files

| File | Description |
|------|-------------|
| `zeta_phi.py` | ζ_φ construction, minima detection, zero matching (importable) |
| `null_test.py` | Window sensitivity + 1000-trial Monte Carlo null test |
| `scaling_test.py` | Matching at N=5000, 10000, 20000 |
| `diagnostic.py` | Per-zero forensic table, produces `per_zero.csv` |
| `per_zero.csv` | 100 rows: each zero, its match, Δ, spacing, null comparison |
| `findings.md` | Narrative write-up of results and recommendation |

## Key finding

The matching does not degrade at high t. It is uniformly modest
(z = -1.4 to -2.2 vs null) across [14, 237]. ζ_φ minima are
~30% closer to Riemann zeros than random, in every range bin.

See `findings.md` for full analysis and recommended correction
to Paper 150 §3.1.1.

## Licence

MIT
