# Pre-registration — Sub C-extension (positive-definite + cubic L-function)

**Date**: 30 May 2026
**Author**: Mr Code
**Brief**: [../BRIEF.md](../BRIEF.md) (CinC, 30 May 2026, Sub C-extension section)
**Binding**: this document is committed to git **before any positive-definite
probe is run and before cubic L-function zeros are pulled**.

## Sub C-ext-1 — positive-definite control

**Probe**: `probe_pos(t) = sqrt(Re²(Z_φ(t)) + 5 · Im²(Z_φ(t)))`
where `Z_φ(t)` is the standard golden-angle Fejér Dirichlet sum
from Paper 150 v2.0 (`zeta_gen.compute_zeta_gen` with α = φ).

This is positive-definite: it has no null cone (no t at which the
form vanishes for general Re, Im); zero only at the origin Re = Im = 0
which is not generically reached in the Dirichlet sum.

**Target**: zeros of L(s, χ_5) cached at
`paper-203-sub-c/zeros_chi5_tmax1000_fast.csv` (904 zeros, t ∈ [6.6, 999.8]).
Reuses Sub C-1's cache verbatim.

**Run parameters** (matching Sub C-1's (5,5) cell exactly):
- N_terms = 5000
- dt = 0.008
- t_pad = 5.0
- W = 1.0
- N_null = 1000
- seed = 42

**Output**: z-score against Poisson null; effect size ε.

## Sub C-ext-2 — cubic-character L-function

**Cubic character χ₂ mod 7** (from brief §5):
- χ₂(0) = 0
- χ₂(1) = 1 = ω⁰
- χ₂(2) = ω = e^(2πi/3)
- χ₂(3) = ω²
- χ₂(4) = ω²
- χ₂(5) = ω
- χ₂(6) = 1

Verification checks before pulling zeros:
- Multiplicativity: χ₂(ab) = χ₂(a)·χ₂(b) for all a, b ∈ {0..6}.
- χ₂(6) = χ₂(−1) = 1 → EVEN character (Γ-factor is Γ(s/2)).
- Sum over residues = 0 (primitive).

**L-function zero finder**: extend `paper-203-sub-c/l_fast.py` to handle
complex characters. The Hardy phase for an EVEN primitive complex
character of conductor q is:
  θ(t, χ) = arg Γ(1/4 + it/2) + (t/2) log(q/π) + (1/2) arg(ε(χ))
where ε(χ) = τ(χ)/√q is the root number. For complex χ, arg(ε(χ))
contributes a constant offset that doesn't affect sign-change locations.

**Implementation choice**: I will rotate the L-function by ε* (the conjugate
root number) so that the rotated function is REAL on the critical line.
Equivalently, compute the completed L-function Λ(s, χ) and use its sign
changes on Re(s)=1/2 as zeros.

The rotation: define `Z_cubic(t) := exp(i·θ_phase(t))·L(½+it, χ₂)` with
appropriate θ_phase including ε(χ) phase, so Z is real-valued on Re(s)=1/2.
Sign changes of Z are zeros of L.

**Target range**: t ∈ (1, 500], matching Sub C-1's d=7 target range.
**Expected zero count**: density ~ (1/2π) log(7t/(2π)) at large t. For
t ≤ 500: about 400 zeros expected.

**Validation**: cross-check the first 3 zeros against published LMFDB
values for L(s, χ₂ mod 7) if available. If not available, sanity-check
that the first zero is at approximately the right t (low-t scaling
log(7·t/2π) crosses 1 at t ≈ π/7·e ≈ 1.21 — first zero typically a few
times that).

**Probe**: ζ_φ,5(t) = sqrt(|Re²(Z_φ(t)) - 5·Im²(Z_φ(t))|), Paper 150 v2.0
canonical (matches the d_probe = 5 row of Sub C-1's matrix).

**Run parameters** (same as Sub C-ext-1; matching Sub C-1's framework).

**Output**: z-score and effect size against the cubic L-function zeros.

## Pre-registered thresholds (re-stated)

**Sub C-ext-1**:
- CONFIRMS-INDEFINITENESS-ESSENTIAL: |z| < 3
- AMBIGUOUS: 3 ≤ |z| < 5
- REFUTES-INDEFINITENESS-CLAIM: |z| ≥ 5

**Sub C-ext-2**:
- NOT DETECTED: |z| < 3
- AMBIGUOUS: 3 ≤ |z| < 5
- DETECTED: |z| ≥ 5

## What I will NOT do

- Switch probe forms after seeing data.
- Adjust the cubic character table after computing L-zeros.
- Reframe the verdict thresholds.

## Files

- `sub_c_ext_1.py` — positive-definite probe driver.
- `sub_c_ext_2.py` — cubic character L-zero finder + probe driver.
- `chi2_mod7.py` — cubic character table + verification.
- `findings_paper_203_sub_c_extension.md` — long-form report.
