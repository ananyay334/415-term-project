# piv_pressure/analysis/summary.py
"""Write results/SUMMARY.md automatically from the study CSV files.

Every number and every observation in the summary is computed here from the
CSVs; nothing is typed in by hand.  Observations are simple, checkable
comparisons (e.g. "Lagrangian error < Eulerian error at the largest dt").

Usage:  python analysis/summary.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from analysis.metrics import convergence_slope  # noqa: E402
from analysis.plots import RESULTS_DIR  # noqa: E402
from studies.common import CASE_LABELS, mean_std, read_csv  # noqa: E402

LABEL = {"poisson": "Poisson", "direct": "Direct (avg.)", "omni": "Omni",
         "direct_path1": "Direct path 1", "direct_path2": "Direct path 2",
         "poisson_velocity": "Poisson (velocity source)", "eulerian": "Eulerian",
         "lagrangian": "Lagrangian", "neumann": "All-Neumann, pin mean",
         "neumann_node": "All-Neumann, pin node", "bernoulli_1": "Bernoulli, 1 edge",
         "exact_1": "Exact p, 1 edge", "bernoulli_2": "Bernoulli, 2 edges",
         "exact_2": "Exact p, 2 edges", "none": "none", "gaussian": "Gaussian",
         "savgol": "Savitzky-Golay", "jacobi": "Jacobi (omega 0.95)", "gs": "Gauss-Seidel",
         "sor": "SOR", "cg": "CG", "jacobi_classic": "Jacobi, classic (omega 1)"}


def _vals(rows: list[dict], key: str, **match) -> list[str]:
    """Sorted distinct values of ``key`` among matching rows."""
    out = {r[key] for r in rows if all(str(r.get(k)) == str(v) for k, v in match.items())}
    try:
        return sorted(out, key=float)
    except ValueError:
        return sorted(out)


def _cases(rows: list[dict]) -> list[str]:
    return sorted({r["case"] for r in rows}, key=lambda c: list(CASE_LABELS).index(c))


def fmt(m: float, s: float | None = None) -> str:
    """'1.23e-02' or '1.23e-02 ± 4.5e-03'."""
    if not np.isfinite(m):
        return "-"
    return f"{m:.3g}" if s is None or s == 0 or not np.isfinite(s) else f"{m:.3g} ± {s:.2g}"


def section_convergence(rows: list[dict]) -> list[str]:
    """Grid-spacing study: tables, fitted slopes, optimum h."""
    out = ["## 1. Grid spacing h (run_convergence.py)", "",
           "Clean = exact point samples; noisy = W = 1 h, sigma = 1 % (mean ± std over seeds).", ""]
    for case in _cases(rows):
        methods = _vals(rows, "method", case=case)
        hs = _vals(rows, "h", case=case)
        out += [f"### {CASE_LABELS[case]}", "",
                "| h | n | " + " | ".join(f"{LABEL[m]} clean" for m in methods) + " | "
                + " | ".join(f"{LABEL[m]} noisy" for m in methods) + " |",
                "|" + "---|" * (2 + 2 * len(methods))]
        for h in sorted(hs, key=float, reverse=True):
            n = _vals(rows, "n_points", case=case, h=h)[0]
            clean = [fmt(mean_std(rows, "eps_l2", case=case, condition="clean", method=m, h=h)[0])
                     for m in methods]
            noisy = [fmt(*mean_std(rows, "eps_l2", case=case, condition="noisy", method=m, h=h)[:2])
                     for m in methods]
            out.append(f"| {float(h):.4g} | {n} | " + " | ".join(clean + noisy) + " |")
        out.append("")
        for m in methods:
            h_c = [float(h) for h in hs if np.isfinite(mean_std(rows, "eps_l2", case=case,
                                                                condition="clean", method=m, h=h)[0])]
            e_c = [mean_std(rows, "eps_l2", case=case, condition="clean", method=m, h=h)[0] for h in h_c]
            slope = convergence_slope(h_c, e_c)
            h_n = [float(h) for h in hs if mean_std(rows, "eps_l2", case=case, condition="noisy",
                                                    method=m, h=h)[2] > 0]
            e_n = [mean_std(rows, "eps_l2", case=case, condition="noisy", method=m, h=h)[0] for h in h_n]
            line = f"- {LABEL[m]}: fitted clean slope = **{slope:.2f}**"
            if e_n:
                k = int(np.argmin(e_n))            # h_n is sorted from smallest to largest h
                line += (f"; noisy optimum h = **{h_n[k]:.4g}** (eps_L2 = {e_n[k]:.3g}); "
                         + ("the noisy error rises again for smaller h." if k > 0
                            else "no rise within the tested range (optimum at the smallest h)."))
            out.append(line)
        out.append("")
    return out


def section_noise(rows: list[dict]) -> list[str]:
    """Noise study: mean ± std vs sigma, growth exponent, timing."""
    out = ["## 2. Noise level sigma (run_noise.py)", "",
           "n = 64, W = 1 h, dt = 0.05; eps_L2 mean ± std over seeds.", ""]
    for case in _cases(rows):
        methods = _vals(rows, "method", case=case)
        sig = _vals(rows, "sigma", case=case)
        out += [f"### {CASE_LABELS[case]}", "", "| sigma (%) | " + " | ".join(LABEL[m] for m in methods) + " |",
                "|" + "---|" * (1 + len(methods))]
        for s in sig:
            out.append(f"| {100 * float(s):g} | " + " | ".join(
                fmt(*mean_std(rows, "eps_l2", case=case, method=m, sigma=s)[:2]) for m in methods) + " |")
        out.append("")
        for m in methods:
            ss = [float(s) for s in sig if float(s) > 0]
            mm = [mean_std(rows, "eps_l2", case=case, method=m, sigma=s)[0] for s in ss]
            t = mean_std(rows, "time_s", case=case, method=m)[0]
            out.append(f"- {LABEL[m]}: log-log slope of eps_L2 vs sigma (sigma > 0) = "
                       f"**{convergence_slope(ss, mm):.2f}** (1 = linear growth); "
                       f"mean wall-clock time {1e3 * t:.1f} ms per reconstruction.")
        best = {m: mean_std(rows, "eps_l2", case=case, method=m, sigma=sig[-1])[0] for m in methods}
        out += [f"- Lowest error at sigma = {100 * float(sig[-1]):g} %: **{LABEL[min(best, key=best.get)]}**.", ""]
    return out


def section_dt(rows: list[dict]) -> list[str]:
    """Time-step study: Eulerian vs Lagrangian."""
    out = ["## 3. Time step dt (run_dt.py)", ""]
    for case in _cases(rows):
        dts = _vals(rows, "dt", case=case)
        out += [f"### {CASE_LABELS[case]}", "",
                "| dt | Eulerian clean | Lagrangian clean | Eulerian noisy | Lagrangian noisy |",
                "|---|---|---|---|---|"]
        for d in dts:
            cells = [fmt(*mean_std(rows, "eps_l2", case=case, condition=c, scheme=s, dt=d)[:2])
                     for c in ("clean", "noisy") for s in ("eulerian", "lagrangian")]
            out.append(f"| {float(d):.3g} | " + " | ".join(cells) + " |")
        dmax = dts[-1]
        e = mean_std(rows, "eps_l2", case=case, condition="clean", scheme="eulerian", dt=dmax)[0]
        lg = mean_std(rows, "eps_l2", case=case, condition="clean", scheme="lagrangian", dt=dmax)[0]
        e0 = mean_std(rows, "eps_l2", case=case, condition="clean", scheme="eulerian", dt=dts[0])[0]
        l0 = mean_std(rows, "eps_l2", case=case, condition="clean", scheme="lagrangian", dt=dts[0])[0]
        out += ["", f"- Clean data, largest dt = {float(dmax):.3g}: Eulerian {e:.3g}, Lagrangian {lg:.3g} -> "
                f"**{'Lagrangian' if lg < e else 'Eulerian'}** is lower; growth from smallest to largest dt: "
                f"Eulerian x{e / e0:.3g}, Lagrangian x{lg / l0:.3g}."]
        for s in ("eulerian", "lagrangian"):
            nm = [mean_std(rows, "eps_l2", case=case, condition="noisy", scheme=s, dt=d)[0] for d in dts]
            if any(np.isfinite(nm)):
                k = int(np.nanargmin(nm))
                out.append(f"- Noisy data, {LABEL[s]}: lowest mean error at dt = {float(dts[k]):.3g} "
                           f"({nm[k]:.3g}).")
        out.append("")
    return out


def section_variants(rows: list[dict], title: str, key: str) -> list[str]:
    """Generic clean/noisy table for the BC and method studies."""
    out = [title, ""]
    for case in _cases(rows):
        variants = list(dict.fromkeys(r[key] for r in rows if r["case"] == case))
        out += [f"### {CASE_LABELS[case]}", "", "| variant | eps_L2 clean | eps_L2 noisy | eps_max noisy | time (ms) |",
                "|---|---|---|---|---|"]
        noisy = {}
        for v in variants:
            c = mean_std(rows, "eps_l2", case=case, condition="clean", **{key: v})[0]
            n, sd, _ = mean_std(rows, "eps_l2", case=case, condition="noisy", **{key: v})
            mx = mean_std(rows, "eps_max", case=case, condition="noisy", **{key: v})[0]
            t = mean_std(rows, "time_s", case=case, **{key: v})[0]
            noisy[v] = n
            out.append(f"| {LABEL.get(v, v)} | {fmt(c)} | {fmt(n, sd)} | {fmt(mx)} | {1e3 * t:.1f} |")
        best = min(noisy, key=lambda v: noisy[v] if np.isfinite(noisy[v]) else np.inf)
        out += ["", f"- Lowest noisy error: **{LABEL.get(best, best)}** ({noisy[best]:.3g}).", ""]
    return out


def section_method_extra(rows: list[dict]) -> list[str]:
    """Path-dependence indicators from the method study."""
    out = []
    for case in _cases(rows):
        pd = mean_std(rows, "path_rms_difference", case=case, condition="noisy", variant="direct")[0]
        sp = mean_std(rows, "direction_spread_rms", case=case, condition="noisy", variant="omni")[0]
        out.append(f"- {CASE_LABELS[case]} (noisy): RMS difference between direct paths 1 and 2 = {pd:.3g}; "
                   f"RMS spread of single-direction omni estimates = {sp:.3g} (pressure units).")
    return out + [""]


def section_window(rows: list[dict]) -> list[str]:
    """Window study: eps and core suction ratio vs W."""
    out = ["## 6. Interrogation window W (run_window.py)", "",
           "Core ratio = reconstructed / reference vortex-core pressure (after mean removal); "
           "< 1 means the suction peak is under-predicted.", ""]
    for case in _cases(rows):
        ws = _vals(rows, "window_h", case=case)
        out += [f"### {CASE_LABELS[case]}", "", "| W / h | eps_L2 clean | core ratio clean | eps_L2 noisy | core ratio noisy |",
                "|---|---|---|---|---|"]
        ratios = []
        for w in ws:
            ec = mean_std(rows, "eps_l2", case=case, condition="clean", window_h=w)[0]
            rc = mean_std(rows, "core_ratio", case=case, condition="clean", window_h=w)[0]
            en, sn, _ = mean_std(rows, "eps_l2", case=case, condition="noisy", window_h=w)
            rn, srn, _ = mean_std(rows, "core_ratio", case=case, condition="noisy", window_h=w)
            ratios.append(rc)
            out.append(f"| {float(w):g} | {fmt(ec)} | {fmt(rc)} | {fmt(en, sn)} | {fmt(rn, srn)} |")
        if all(np.isfinite(ratios)):
            mono = all(b <= a + 1e-12 for a, b in zip(ratios, ratios[1:]))
            out += ["", f"- Clean core ratio goes from {ratios[0]:.4g} (W = {float(ws[0]):g} h) to "
                    f"{ratios[-1]:.4g} (W = {float(ws[-1]):g} h); "
                    f"{'decreases monotonically with W' if mono else 'not monotonic in W'}."]
        out.append("")
    return out


def section_denoise(rows: list[dict]) -> list[str]:
    """Denoising study: eps per sigma and denoiser; where it starts to help."""
    out = ["## 7. Denoising (run_denoise.py)", ""]
    for case in _cases(rows):
        sig = _vals(rows, "sigma", case=case)
        dn = ["none", "gaussian", "savgol"]
        out += [f"### {CASE_LABELS[case]}", "", "| sigma (%) | " + " | ".join(LABEL[d] for d in dn) + " |",
                "|" + "---|" * 4]
        for s in sig:
            out.append(f"| {100 * float(s):g} | " + " | ".join(
                fmt(*mean_std(rows, "eps_l2", case=case, denoise=d, sigma=s)[:2]) for d in dn) + " |")
        out.append("")
        for d in ("gaussian", "savgol"):
            helps = [float(s) for s in sig if mean_std(rows, "eps_l2", case=case, denoise=d, sigma=s)[0]
                     < mean_std(rows, "eps_l2", case=case, denoise="none", sigma=s)[0]]
            hurts = [float(s) for s in sig if float(s) not in helps]
            out.append(f"- {LABEL[d]}: lower error than no denoising at sigma = "
                       f"{', '.join(f'{100 * s:g} %' for s in helps) or 'none'}; higher at "
                       f"{', '.join(f'{100 * s:g} %' for s in hurts) or 'none'}.")
        out.append("")
    return out


def section_residuals(results: Path) -> list[str]:
    """Iterative solvers: iterations, times, growth with N, SOR stability, norms."""
    out = ["## 9. Iterative Poisson solvers and residual norms (run_residuals.py)", ""]
    hist_path = results / "residual_history.csv"
    if hist_path.is_file():
        rows = read_csv(hist_path)
        solvers = list(dict.fromkeys(r["solver"] for r in rows))
        n = rows[0]["n"]
        out += [f"Taylor-Green {n}x{n}, all-Neumann (zero mean), tol on the scaled L2_n residual.", "",
                "| solver | omega | iterations | wall-clock (s) | converged | final scaled residual |",
                "|---|---|---:|---:|---|---:|"]
        last = {}
        for sv in solvers:
            r = [x for x in rows if x["solver"] == sv][-1]
            last[sv] = r
            out.append(f"| {LABEL.get(sv, sv)} | {float(r['omega']):.4g} | {r['iteration']} | "
                       f"{float(r['time_s']):.3g} | {r['converged']} | {float(r['residual_scaled']):.2e} |")
        conv = {k: int(v["iteration"]) for k, v in last.items() if v["converged"] == "True"}
        if "jacobi" in conv:
            for sv in ("gs", "sor", "cg"):
                if sv in conv:
                    out.append(f"- {LABEL[sv]} needed {conv['jacobi'] / conv[sv]:.3g}x fewer iterations than "
                               f"Jacobi (omega = 0.95).")
        if "jacobi_classic" in last and last["jacobi_classic"]["converged"] != "True":
            out.append(f"- Classic Jacobi (omega = 1) did not reach the tolerance: its residual levelled off "
                       f"at {float(last['jacobi_classic']['residual_scaled']):.2e}.")
        out.append("")
    grid_path = results / "iterations_vs_grid.csv"
    if grid_path.is_file():
        rows = read_csv(grid_path)
        ns = sorted({int(r["n"]) for r in rows})
        methods = list(dict.fromkeys(r["method"] for r in rows))
        out += ["| n | N | " + " | ".join(f"{LABEL.get(m, m)} it. (s)" for m in methods) + " |",
                "|" + "---|" * (2 + len(methods))]
        for n in ns:
            cells = []
            for m in methods:
                r = [x for x in rows if int(x["n"]) == n and x["method"] == m]
                cells.append(f"{r[0]['iterations']} ({float(r[0]['time_s']):.2g})" if r else "-")
            out.append(f"| {n} | {n * n} | " + " | ".join(cells) + " |")
        out.append("")
        for m in methods:
            r = [x for x in rows if x["method"] == m and x["converged"] == "True"]
            if len(r) > 1:
                a = convergence_slope([float(x["N"]) for x in r], [float(x["iterations"]) for x in r])
                out.append(f"- {LABEL.get(m, m)}: iterations grow like N^{a:.2f} "
                           f"(fitted over n = {', '.join(x['n'] for x in r)}).")
        out.append("")
    sor_path = results / "sor_omega.csv"
    if sor_path.is_file():
        rows = [r for r in read_csv(sor_path) if r["converged"] == "True"]
        if rows:
            best = min(rows, key=lambda r: int(r["iterations"]))
            out.append(f"- SOR: fewest iterations ({best['iterations']}) at omega = {float(best['omega']):.4g} "
                       f"among the tested values.")
        div = results / "sor_divergence.csv"
        if div.is_file():
            drows = read_csv(div)
            for w in dict.fromkeys(r["omega"] for r in drows):
                rr = [r for r in drows if r["omega"] == w]
                r0, r1 = float(rr[0]["residual_scaled"]), float(rr[-1]["residual_scaled"])
                state = "diverged" if r1 > 1e3 * r0 else "no convergence"
                out.append(f"- SOR with omega = {float(w):g}: scaled residual went from {r0:.3g} to "
                           f"{r1:.3g} after {rr[-1]['iteration']} iterations ({state}).")
        out.append("")
    ns_path = results / "norm_scaling.csv"
    if ns_path.is_file():
        rows = read_csv(ns_path)
        N = [float(r["N"]) for r in rows]
        out.append("- Norms of the same initial residual (interior rows) on grids " +
                   ", ".join(r["n"] for r in rows) + ": growth exponents with N - " + ", ".join(
                       f"{k.upper()} {convergence_slope(N, [float(r[k]) for r in rows]):.2f}"
                       for k in ("l1", "l2", "l1n", "l2n")) + ".")
        out.append("- Including the ghost-point boundary rows (flux terms 2G/h) the L2_n exponent is "
                   f"{convergence_slope(N, [float(r['l2n_with_boundary_rows']) for r in rows]):.2f}.")
        out.append("")
    return out


def section_order(results: Path) -> list[str]:
    """First- vs second-order grid convergence: slopes and the 1 % resolution."""
    path = results / "order_summary.csv"
    if not path.is_file():
        return []
    rows = read_csv(path)
    out = ["## 10. First- vs second-order scheme (run_order.py)", "",
           "Slopes of error (%) vs dx; 'fine' = fit over the finest grids (asymptotic range).", "",
           "| case | order | slope (fine) | slope (all grids) | dx where error first < 1 % | "
           "1 % crossing (interpolated) | noisy: min error % at dx |", "|---|---|---:|---:|---:|---:|---|"]
    for r in rows:
        cross = float(r["dx_cross"])
        note = " (extrapolated)" if r["cross_extrapolated"] == "True" else ""
        first = float(r["dx_first_below"])
        out.append(f"| {CASE_LABELS.get(r['case'], r['case'])} | {r['order']} | {float(r['slope_fine']):.2f} | "
                   f"{float(r['slope_all']):.2f} | {first:.4g} | {cross:.4g}{note} | "
                   f"{float(r['noisy_min_err_pct']):.3g} % at {float(r['noisy_optimum_dx']):.3g} |")
    return out + [""]


def section_cfd(results: Path) -> list[str]:
    """CFD grid independence and validation (if present)."""
    out = ["## 8. CFD cylinder wake", ""]
    gi = results / "cfd" / "grid_independence.md"
    out += (["### Grid independence (cfd/grid_independence.py)", "", gi.read_text(), ""]
            if gi.is_file() else ["CFD grid-independence results not found (run cfd/run_cases.sh).", ""])
    tl = results / "cfd" / "time_loop_summary.json"
    summ = json.loads(tl.read_text()) if tl.is_file() else {}
    if summ.get("periodicity") or "window_change" in summ or "steady" in summ:
        out += ["### Time-loop criteria (cfd/change_norm.py)", ""]
        for lev, d in summ.get("periodicity", {}).items():
            out.append(f"- {lev}: {d['n_cycles']} shedding cycles; cycle-to-cycle change of Cl amplitude and "
                       f"period < {summ['tol']:g} for 3 cycles from t = {d['t_periodic']:.4g} "
                       f"(satisfied: {d['satisfied']}).")
        if "window_change" in summ:
            w = summ["window_change"]
            out.append(f"- Re = 100 window: ||u^(n+1) - u^n|| stays between {w['min']:.3g} and {w['max']:.3g} "
                       "(periodic flow - no steady state).")
        if "steady" in summ:
            st = summ["steady"]
            out.append(f"- Steady Re = 20 case: change norm first below eps = {st['eps']:g} at t = "
                       f"{st['t_below_eps']:.4g}; the OpenFOAM run stopped at t = {st['t_stop']:.4g}.")
        out.append("")
    val = results / "cfd_validation.csv"
    if val.is_file():
        rows = read_csv(val)
        out += ["### Reconstruction of one snapshot (run_cfd_validation.py)", "",
                "| method | eps_L2 | eps_max | core ratio |", "|---|---|---|---|"]
        for r in rows:
            if r.get("kind") == "snapshot":
                out.append(f"| {r['method']} | {float(r['eps_l2']):.3g} | {float(r['eps_max']):.3g} | "
                           f"{float(r['core_ratio']):.3g} |")
        for s in ("eulerian", "lagrangian"):
            m, sd, k = mean_std(rows, "eps_l2", kind="history", scheme=s)
            if k:
                out.append(f"\n- Over one shedding period ({k} snapshots), {LABEL[s]}: eps_L2 = {fmt(m, sd)}.")
        out.append("")
    return out


def main(results: Path = RESULTS_DIR) -> Path:
    """Build SUMMARY.md from whatever CSVs exist."""
    lines = ["# Results summary (generated automatically)", "",
             f"Generated {datetime.now():%Y-%m-%d %H:%M} by `analysis/summary.py` from the CSV files in "
             "`results/`. All numbers below come from those files; re-running the studies regenerates them.", ""]
    parts = [
        ("convergence.csv", section_convergence), ("noise.csv", section_noise), ("dt.csv", section_dt),
        ("bc.csv", lambda r: section_variants(r, "## 4. Boundary conditions (run_bc.py)", "variant")),
        ("method.csv", lambda r: section_variants(r, "## 5. Integration methods (run_method.py)", "variant")
         + section_method_extra(r)),
        ("window.csv", section_window), ("denoise.csv", section_denoise),
    ]
    for fname, fn in parts:
        path = results / fname
        lines += fn(read_csv(path)) if path.is_file() else [f"_{fname} not found - study not run._", ""]
    lines += section_cfd(results)
    lines += section_residuals(results)
    lines += section_order(results)
    out = results / "SUMMARY.md"
    out.write_text("\n".join(lines))
    print(f"[summary] -> {out}")
    return out


if __name__ == "__main__":
    main()
