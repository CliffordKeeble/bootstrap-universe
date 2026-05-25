"""
Task 2 orchestrator: continued-fraction table for 8 candidates.

Pre-reg: Mr_Code_Brief_Paper_190_Phase_1_v0_2.md (git commit 971938b).

Candidates (per brief sec 4 Task 2):
  a) gamma                             - Euler-Mascheroni
  b) 1/sqrt(3)                          - the v0.1 1/sqrt 3 quantity
  c) log(phi)                           - natural icosahedral log
  d) pi/(3 sqrt 3) = L(1, chi_{-3})     - Eisenstein L-value
  e) 2 log(phi)/sqrt(5) = L(1, chi_5)   - Dedekind L-value for Q(sqrt 5)
  f) gamma_0(Q(sqrt 5))                 - Dedekind Stieltjes, real-quadratic
  g) gamma_0(Q(sqrt -3))                - Dedekind Stieltjes, imag-quadratic
  h) gamma_{2I-spectral}                - Stieltjes-like for icosahedral spectrum

Method
  Rows (a)-(g): compute value at primary precision (2000 dps) and cross-check
  precision (1500 dps); run CF in lockstep, record stable prefix.
  Row (h): two Richardson fits at different N-range (sanity-check at full N,
  and reduced N) used as primary/cross precision pair. CF length capped by
  N^{-1} convergence to ~15-25 terms.

Outputs
  task_2_cf_data.json          : machine-readable raw CFs + values
  task_2_cf_table.md           : the unified markdown table
"""

import json
import time
from pathlib import Path

from mpmath import mp, mpf, log, sqrt, pi, euler, nstr

from cf_tools import continued_fraction, cf_with_precision_check, convergents
from dedekind_stieltjes import chi_5, chi_minus_3, gamma_0_method_B
from spectral_2I import (
    compute_molien_coefficients, compute_partial_sum_S_N_float,
    fit_log_slope_richardson,
)


DPS_PRIMARY = 2000
DPS_CROSS = 1500
DPS_ITER = 2200
N_CF = 300
M_CONV = 10

ROOT = Path(__file__).parent
JSON_OUT = ROOT / "task_2_cf_data.json"
MD_OUT = ROOT / "task_2_cf_table.md"


def compute_with_two_precisions(name, factory, dps_high=DPS_PRIMARY, dps_low=DPS_CROSS):
    """Compute factory() at two precisions; return (x_high, x_low) as mpf
    at iteration precision DPS_ITER."""
    print(f"[{name}] primary at {dps_high} dps...", flush=True)
    t0 = time.time()
    mp.dps = dps_high
    x_high_str = mp.nstr(factory(), dps_high + 5, strip_zeros=False)
    t1 = time.time()
    print(f"[{name}]   {t1-t0:.1f}s, first 30 chars: {x_high_str[:30]}...", flush=True)

    print(f"[{name}] cross at {dps_low} dps...", flush=True)
    t0 = time.time()
    mp.dps = dps_low
    x_low_str = mp.nstr(factory(), dps_low + 5, strip_zeros=False)
    t1 = time.time()
    print(f"[{name}]   {t1-t0:.1f}s, first 30 chars: {x_low_str[:30]}...", flush=True)

    mp.dps = DPS_ITER
    x_high = mpf(x_high_str)
    x_low = mpf(x_low_str)
    return x_high, x_low, x_high_str, x_low_str


def cf_block(name, factory):
    x_h, x_l, sh, sl = compute_with_two_precisions(name, factory)
    print(f"[{name}] CF iteration (N <= {N_CF}, iter dps = {DPS_ITER})...", flush=True)
    t0 = time.time()
    cf_stable, n_fail, info = cf_with_precision_check(x_h, x_l, N=N_CF)
    t1 = time.time()
    print(f"[{name}]   {t1-t0:.1f}s, stable length = {len(cf_stable)}", flush=True)
    if info is not None:
        print(f"[{name}]   precision-stability info: {info}", flush=True)
    convs = convergents(cf_stable, m=M_CONV)
    return {
        "name": name,
        "value_high_dps": DPS_PRIMARY,
        "value_low_dps": DPS_CROSS,
        "value_high_str": sh,
        "value_low_str": sl,
        "cf": list(cf_stable),
        "cf_length": len(cf_stable),
        "n_fail": n_fail,
        "precision_info": info,
        "convergents": [[c.numerator, c.denominator] for c in convs],
    }


def spectral_cf_block():
    """Row (h): Richardson at two N-ranges as primary/cross precision pair."""
    name = "h) gamma_{2I-spectral}"
    print(f"\n[{name}] Computing Molien coefficients (N_max = 1e7)...", flush=True)
    t0 = time.time()
    m = compute_molien_coefficients(10**7)
    t1 = time.time()
    print(f"[{name}]   {t1-t0:.1f}s", flush=True)

    print(f"[{name}] Partial sums at multiple N...", flush=True)
    N_vals_full = [10**4, 10**5, 10**6, 10**7]
    S_vals_full = [compute_partial_sum_S_N_float(m, N) for N in N_vals_full]
    for N, S in zip(N_vals_full, S_vals_full):
        print(f"[{name}]   S_{N} = {S:.15f}", flush=True)

    fit_full = fit_log_slope_richardson(N_vals_full, S_vals_full, mp_dps=30)
    fit_short = fit_log_slope_richardson(N_vals_full[:-1], S_vals_full[:-1], mp_dps=30)
    g_high = fit_full["gamma"]
    g_low = fit_short["gamma"]
    print(f"[{name}]   gamma (full fit, N up to 1e7) = {nstr(g_high, 12)}", flush=True)
    print(f"[{name}]   gamma (short fit, N up to 1e6) = {nstr(g_low, 12)}", flush=True)
    print(f"[{name}]   |full - short| = {nstr(abs(g_high - g_low), 6)}", flush=True)

    mp.dps = DPS_ITER
    x_h = mpf(g_high)
    x_l = mpf(g_low)
    cf_stable, n_fail, info = cf_with_precision_check(x_h, x_l, N=N_CF)
    print(f"[{name}]   CF stable length = {len(cf_stable)} "
          f"(N^{{-1}} convergence limits precision)", flush=True)
    if info is not None:
        print(f"[{name}]   info: {info}", flush=True)
    convs = convergents(cf_stable, m=M_CONV)

    return {
        "name": name,
        "value_high_dps": "Richardson N=1e7",
        "value_low_dps": "Richardson N=1e6",
        "value_high_str": nstr(g_high, 20),
        "value_low_str": nstr(g_low, 20),
        "fit_c_empirical": nstr(fit_full["c"], 15),
        "fit_c_predicted": nstr(mpf(1)/120, 15),
        "cf": list(cf_stable),
        "cf_length": len(cf_stable),
        "n_fail": n_fail,
        "precision_info": info,
        "convergents": [[c.numerator, c.denominator] for c in convs],
    }


def main():
    print("=" * 70)
    print("Task 2 CF table: 8 candidates")
    print(f"  Primary dps = {DPS_PRIMARY}, Cross dps = {DPS_CROSS}, Iter dps = {DPS_ITER}")
    print(f"  N_cf target = {N_CF}, first {M_CONV} convergents reported per row.")
    print("=" * 70)

    results = {}

    # (a) gamma
    results["a"] = cf_block("a) gamma", lambda: +euler)
    # (b) 1/sqrt(3)
    results["b"] = cf_block("b) 1/sqrt(3)", lambda: mpf(1) / sqrt(mpf(3)))
    # (c) log(phi)
    results["c"] = cf_block("c) log(phi)",
                            lambda: log((mpf(1) + sqrt(mpf(5))) / mpf(2)))
    # (d) pi/(3 sqrt 3) = L(1, chi_{-3})
    results["d"] = cf_block("d) pi/(3 sqrt 3)",
                            lambda: pi / (mpf(3) * sqrt(mpf(3))))
    # (e) 2 log(phi)/sqrt(5) = L(1, chi_5)
    results["e"] = cf_block("e) 2 log(phi)/sqrt(5)",
                            lambda: mpf(2) * log((mpf(1) + sqrt(mpf(5))) / mpf(2))
                                    / sqrt(mpf(5)))
    # (f) gamma_0(Q(sqrt 5)) via Method B
    def g_f():
        R = mpf(2) * log((mpf(1) + sqrt(mpf(5))) / mpf(2)) / sqrt(mpf(5))
        return gamma_0_method_B(chi_5, 5, L_value=R, name="Q(sqrt 5)")["gamma_0_K"]
    results["f"] = cf_block("f) gamma_0(Q(sqrt 5))", g_f)

    # (g) gamma_0(Q(sqrt -3)) via Method B
    def g_g():
        R = pi / (mpf(3) * sqrt(mpf(3)))
        return gamma_0_method_B(chi_minus_3, 3, L_value=R, name="Q(sqrt -3)")["gamma_0_K"]
    results["g"] = cf_block("g) gamma_0(Q(sqrt -3))", g_g)

    # (h) gamma_{2I-spectral}
    results["h"] = spectral_cf_block()

    # ---- Save JSON ----
    out = {
        "metadata": {
            "pre_reg_commit": "971938b",
            "dps_primary": DPS_PRIMARY,
            "dps_cross": DPS_CROSS,
            "dps_iter": DPS_ITER,
            "N_cf_target": N_CF,
            "M_convergents_reported": M_CONV,
            "task": "Task 2 of Mr_Code_Brief_Paper_190_Phase_1_v0_2.md",
        },
        "candidates": results,
    }
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {JSON_OUT}", flush=True)

    # ---- Build markdown table ----
    md = []
    md.append("# Task 2 - Continued-fraction table for 8 candidates")
    md.append("")
    md.append(f"**Pre-reg commit:** `971938b`")
    md.append(f"**Primary precision:** {DPS_PRIMARY} dps  ")
    md.append(f"**Cross-check precision:** {DPS_CROSS} dps  ")
    md.append(f"**Iteration precision:** {DPS_ITER} dps  ")
    md.append(f"**Target N:** {N_CF} partial quotients per row")
    md.append("")
    md.append("## Summary")
    md.append("")
    md.append("| row | candidate | value (30 sig fig) | CF length | n_fail |")
    md.append("|-----|-----------|---------------------|-----------|--------|")
    for key in ["a", "b", "c", "d", "e", "f", "g", "h"]:
        r = results[key]
        val30 = r["value_high_str"][:32]
        md.append(f"| {key} | {r['name']} | `{val30}...` | {r['cf_length']} | {r['n_fail']} |")
    md.append("")

    md.append("## Continued fractions (first 30 partial quotients per row)")
    md.append("")
    md.append("| row | CF prefix |")
    md.append("|-----|-----------|")
    for key in ["a", "b", "c", "d", "e", "f", "g", "h"]:
        r = results[key]
        cf30 = r["cf"][:30]
        if cf30:
            head = str(cf30[0])
            tail = ", ".join(str(a) for a in cf30[1:])
            cfstr = f"[{head}; {tail}]" if tail else f"[{head}]"
        else:
            cfstr = "(empty)"
        md.append(f"| {key} | `{cfstr}` |")
    md.append("")

    md.append("## First 10 convergents per row")
    md.append("")
    for key in ["a", "b", "c", "d", "e", "f", "g", "h"]:
        r = results[key]
        md.append(f"### Row ({key}) {r['name']}")
        md.append("")
        md.append("| k | p_k | q_k | p_k/q_k |")
        md.append("|---|-----|-----|---------|")
        for k, (p, q) in enumerate(r["convergents"]):
            md.append(f"| {k} | {p} | {q} | {p}/{q} |")
        md.append("")

    md.append("## Full CFs (machine-readable in task_2_cf_data.json)")
    md.append("")
    for key in ["a", "b", "c", "d", "e", "f", "g", "h"]:
        r = results[key]
        md.append(f"### Row ({key}) {r['name']} - {r['cf_length']} partial quotients")
        md.append("")
        md.append("```")
        # Wrap at 80 chars
        line = ""
        for i, a in enumerate(r["cf"]):
            sep = "; " if i == 1 else ", " if i > 0 else ""
            piece = f"{sep}{a}"
            if len(line) + len(piece) > 78:
                md.append(line)
                line = piece.lstrip(",; ")
            else:
                line += piece
        if i == 0:
            line = f"[{r['cf'][0]}]"
        else:
            line = "[" + line if not line.startswith("[") else line
            if not line.endswith("]"):
                line = line + "]"
        md.append(line)
        md.append("```")
        md.append("")

    md.append("## Precision-stability notes")
    md.append("")
    md.append("Rows (a)-(g): primary 2000 dps vs cross 1500 dps; n_fail records "
              "the first partial quotient where the two precisions disagree (== N_CF "
              "if all stable). Per CinC refinement, when precision-stability fails, "
              "both the disagreement AND the stable value at n_fail-1 are recorded "
              "in `precision_info`.")
    md.append("")
    md.append("Row (h) gamma_{2I-spectral}: no analytic continuation available; "
              "Richardson at N up to 1e7 (primary) vs N up to 1e6 (cross) used as "
              "precision-stability pair. CF length is capped by the precision of "
              "the fit, which is bounded by the 1/N convergence rate of the partial "
              "sum (faster than the brief's N^{-1/2} assumption).")
    md.append("")

    with open(MD_OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"Wrote {MD_OUT}", flush=True)


if __name__ == "__main__":
    main()
