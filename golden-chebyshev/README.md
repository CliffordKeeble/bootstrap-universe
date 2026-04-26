# golden-chebyshev/

Companion to **Paper 196 (candidate)** — *Chebyshev Bias and Golden Phase Structure in L(χ₅)*.

Three-layer numerical investigation of the Dedekind factorisation
ζ_{ℚ(√5)}(s) = ζ(s) · L(s, χ₅) viewed through the chirality channel L(s, χ₅).
Layer 1 reproduces the classical Chebyshev bias for q = 5; Layer 2 validates the
explicit-formula toolchain by reconstructing the bias from L(χ₅) zeros; Layer 3
(held — pending CinC redesign) tests whether the arctan(1/φ) phase lock from
Paper 125 inherits to L(χ₅) zeros via the Hurwitz N = 120 mod-5 character
decomposition.

## Layout

| File | Purpose | Layer | Status flag |
|---|---|---|---|
| `chebyshev_bias.py` | Direct count of E(x) = #stubborns − #splitters mod 5 | 1 | OBSERVED |
| `lchi5_zeros.py` | L(χ₅) zero data: LMFDB-anchored + mpmath-refined | 2 setup | DERIVED |
| `explicit_formula_lchi5.py` | Spectral reproduction of E(x) via L(χ₅) zeros | 2 | DERIVED |
| `density_qch5.py` | R-S logarithmic density at q = 5, three paths | 1+2 | mixed |
| `test_chebyshev.py` | Sieve correctness + sign convention + sanity vs Granville-Martin Table 4 | 1 | — |
| `test_explicit.py` | Spectral reconstruction agreement + truncation sensitivity | 2 | — |
| `paper_196_results.md` | Findings markdown with status-flagged claims | all | — |

## Sign convention

E(x) = #{p ≤ x : p ≡ ±2 mod 5} − #{p ≤ x : p ≡ ±1 mod 5}

Stubborns minus splitters (non-residues minus residues). Positive E(x) = stubborn
lead. Confirmed via Granville-Martin 2006 §3 (equation 3 framing).

## Layer 3 status

**Held.** Layer 3 redesigned post-Paper-125 read (brief addendum 2). When run, it
will need a new module `hurwitz_120_mod5.py` and the full Hurwitz N=120 mod-5
character decomposition machinery. Not started here.

## References

- Granville & Martin, *Prime Number Races*, Amer Math Monthly 113 (2006) 1–33
- Rubinstein & Sarnak, *Chebyshev's bias*, Experiment. Math. 3 (1994) 173–197
- Stark, Acta Arith. 68 (1971) 311–320 — mod 5 sign-change
- Keeble, *The Golden Phase Lock at Zeta Zeros*, Bootstrap Universe Programme
  Paper 125. DOI: 10.5281/zenodo.19022277
- Keeble, *The Golden Norm Hears the Dedekind Zeta*, Bootstrap Universe Programme
  Paper 150 v2.0. DOI: 10.5281/zenodo.19162130 (sibling investigation, golden-zeta/)

## Reuse from sibling directories

- `golden-zeta/run_dedekind_partition.py:45` — `chi5(n)` Kronecker (n/5)
  implementation. Reused via in-module reimplementation (5 lines, not worth a
  cross-directory import for a Legendre symbol).
- `golden-zeta/zeta_phi.py` — `find_minima`, `match_zeros` — **not used here**
  (Layer 2 uses LMFDB high-precision zeros, not minima detection).

## Requirements

See `requirements.txt`. Standard scientific Python — numpy, mpmath, sympy,
matplotlib. No exotic dependencies.
