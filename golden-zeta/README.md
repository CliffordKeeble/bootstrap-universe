# golden-zeta

Computational companion to **Paper 150 v2.0 — *The Golden Norm Hears the Dedekind Zeta***
DOI: [10.5281/zenodo.19162130](https://doi.org/10.5281/zenodo.19162130) (all-versions)
Author: Clifford Keeble · ORCID [0009-0003-6828-2155](https://orcid.org/0009-0003-6828-2155)

---

## Summary

The Dedekind factorisation ζ<sub>ℚ(√5)</sub>(s) = ζ(s) · L(s, χ<sub>5</sub>) (Dedekind, 1871) makes the Riemann zeros one factor of the Dedekind zeta function of the smallest non-trivial real quadratic field. This directory tests the hypothesis that those zeros are detectable as statistical minima of a Fejér-weighted Dirichlet sum under the golden-field norm N(a + bφ) = a² + ab − b² on ℤ[φ].

At N = 10,000 Riemann zeros, the test statistic reaches **z = −22.91** against a same-density uniform-random null. Discovery threshold is |z| ≥ 5. The investigation was executed under a pre-registered four-test falsification protocol on 12 April 2026; all four tests passed, with mechanism isolated to the field norm rather than to φ specifically.

## The four tests and what they showed

| Test | Question | Result |
|---|---|---|
| 1. Primary | Does the golden-norm detector separate Riemann zeros from random? | z = −22.91 at N = 10,000 |
| 2. φ-specificity | Is φ the essential phase? | No. √2 phase gives \|z\| = 9.38; e and π both exceed \|z\| > 6. Any equidistributed irrational works. |
| 3. Norm decoupling | Is the golden field norm essential? | Yes. φ-phase with circular norm fails (\|z\| = 2.63, sub-threshold). √2-phase with golden norm passes (\|z\| = 9.38). The norm carries the effect. |
| 4. Dedekind partition | Are both factors of ζ<sub>ℚ(√5)</sub> independently detected? | Yes. Riemann zeros at \|z\| = 8.14, L(χ<sub>5</sub>) zeros at \|z\| = 6.18. |

The headline is not that φ is special. The headline is that the **indefinite form p² − 5q² of ℚ(√5) is an empirical Dedekind probe**. Any sufficiently dense equidistributed phase sequence excites it.

## Pre-registration audit trail

Two commits establish the protocol before any analysis was run:

- `21f74ed` — primary protocol ([view commit](https://github.com/CliffordKeeble/bootstrap-universe/commit/21f74ed))
  Committed 2026-04-12 08:24:28 UTC

- `c2b9d0d` — compound addendum (Tests 2, 3, 4) ([view commit](https://github.com/CliffordKeeble/bootstrap-universe/commit/c2b9d0d))
  Committed 2026-04-12 08:49:09 UTC

Every dependent analysis commit is timestamped 12+ minutes after its pre-registration parent. Linear history, no rewrites. The chronology is inspectable.

This matters because the headline test statistic is only meaningful against a protocol that was fixed *before* the analyses were run. That protocol is in `pre_registration_N1000.md` and `pre_registration_compound.md` at the SHAs above.

## Directory contents

### Pre-registration
- `pre_registration_N1000.md` — primary protocol (commit `21f74ed`)
- `pre_registration_compound.md` — compound addendum with four falsification tests (commit `c2b9d0d`)

### Analysis
- `run_N1000.py` — primary test at N = 1000
- `run_N10000.py` — primary test at N = 10000 (headline result)
- `run_alt_alpha.py` — Test 2 (φ-specificity)
- `run_norm_decoupling.py` — Test 3 (norm essentiality)
- `run_dedekind_partition.py` — Test 4 (both factors of Dedekind product)

### Forensic data (per-zero output tables)
- `per_zero_1000.csv`, `per_zero_N10000.csv` — primary-test outputs
- `per_zero_sqrt2.csv`, `per_zero_e.csv`, `per_zero_pi.csv` — alternative-phase probes
- `per_zero_circular_norm.csv`, `per_zero_sqrt2_goldenNorm.csv` — norm-decoupling probes
- `per_zero_Lchi5.csv` — L(χ<sub>5</sub>) detection table

### Narrative findings
- `findings_N1000.md`
- `findings_compound.md`

## Reproduction

```bash
git clone https://github.com/CliffordKeeble/bootstrap-universe.git
cd bootstrap-universe/golden-zeta
pip install -r requirements.txt
python run_N1000.py              # primary test, smaller N, faster turnaround
python run_N10000.py             # headline test — longer runtime; see note
python run_alt_alpha.py          # Test 2
python run_norm_decoupling.py    # Test 3
python run_dedekind_partition.py # Test 4
```

**Runtime note.** The primary result at N = 1000 completes in minutes on a modern laptop. The headline result at N = 10,000 is substantially longer because `mpmath.zetazero` computes each Riemann zero from scratch; expect this to run into several hours depending on hardware. For verification of the protocol and mechanism rather than the headline scale, `run_N1000.py` is the recommended entry point.

Riemann zero data is computed on the fly via `mpmath`; no external dataset download is required. LMFDB zeros are cited in the paper for cross-reference but are not a runtime dependency.

## Licence

MIT. See [`../LICENSE`](../LICENSE).

## Correspondence

Dr. Clifford Keeble, PhD — ORCID [0009-0003-6828-2155](https://orcid.org/0009-0003-6828-2155)
[academia.edu/CliffordKeeble](https://independent.academia.edu/CliffordKeeble)
Woodbridge, UK

🐕☕⬡
