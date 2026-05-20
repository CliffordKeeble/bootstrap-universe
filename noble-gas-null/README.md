# noble-gas-null — Paper 17 v2.0 null forensic

Companion code for Paper 17 v2.0's §3 noble-gas retirement.

## Provenance

Written for **Paper 17 v2.0 Notice of Supersession** of the icosahedral proton periodic table (concept DOI [10.5281/zenodo.18145312](https://doi.org/10.5281/zenodo.18145312), v2.0 issued May 2026). The script computes the load-bearing figure cited in §3: how cheap it is to hit the six noble-gas atomic numbers {2, 10, 18, 36, 54, 86} from the icosahedral atom inventory under simple arithmetic.

## Question

Does the icosahedral integer-combination space at expression depth *k* cover the six noble-gas atomic numbers?

- **Atom inventory:** `{V=12, F=20, E=30, χ=2, D=3, D!=6}` — icosahedral vertex/face/edge counts, Euler characteristic, dimensional count, dimensional factorial.
- **Operations:** `{+, −, ×, ·²}` — binary add/subtract/multiply, plus unary squaring.
- **Depth convention:** parse-tree depth, not operator count. Every operand is an inventory member; no free integers are introduced. The script's `Level k` corresponds to the paper's `Depth k`.

If the expressible space already covers most integers ≤ 100 at low depth, then hitting all six noble gases is not evidence of structure — it is what any dense expression space would produce. That is the diagnosis Paper 17 v2.0 cites.

## Expected output (for re-run verification)

| Depth | Integers ≤ 100 expressible | Noble gases hit |
|:----:|:----:|:----|
| 1 | 31 / 100 | {2, 10, 18, 36} — 4 of 6 |
| 2 | **98 / 100** | **{2, 10, 18, 36, 54, 86} — all 6** |
| 3 | 100 / 100 | all 6 |

At depth 2, Xe = 54 is reached via D! × D² and Rn = 86 via F + E + (D!)². The "98/100 at depth 2, 100/100 at depth 3" pair is the figure Paper 17 v2.0 cites as the supersession's empirical foundation.

A canonical run transcript is saved alongside the script as `run_output.txt`.

## How to run

```
python noble_gas_null.py
```

Stdlib only — no external dependencies. Exits cleanly with the table above printed to stdout.

**Windows note:** the script prints the character `≤`. On a default Windows console (cp1252) Python's `print` will raise `UnicodeEncodeError`. Set UTF-8 stdout first:

```
set PYTHONIOENCODING=utf-8 && python noble_gas_null.py     # cmd
$env:PYTHONIOENCODING='utf-8'; python noble_gas_null.py     # PowerShell
PYTHONIOENCODING=utf-8 python noble_gas_null.py             # bash / WSL / macOS / Linux
```

Linux and macOS typically default to UTF-8 and run cleanly without the prefix.

## Honest note — retrospective archival

This is a **retrospective code commit**, not a pre-registered computation. The figure (98/100 at depth 2, 100/100 at depth 3) was used in Paper 17 v2.0 prior to this commit; the script existed at the time the paper was written but was not in the repo. Standard Mr Code pre-registration discipline (commit protocol → commit code → run → commit results) does not apply here. The commit exists for **reproducibility of a load-bearing figure**, so the supersession claim in §3 of Paper 17 v2.0 remains re-runnable on demand.

This framing is recorded explicitly because the programme distinguishes pre-registered numerical work from retrospective archival, and the audit trail should be honest about which this is.
