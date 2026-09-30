# piv_pressure/studies/run_convergence.py
"""Grid-spacing study: error vs PIV spacing h (Figures 1 and 2).

* clean data  = exact point samples (no window averaging, no noise):
  the error is pure discretisation error and should fall as h^2;
* noisy data  = baseline PIV model (W = 1 h, sigma = 1 %, 10 seeds):
  noise amplified by differentiation makes the error rise again at small h,
  so there is an optimum spacing, which is reported.

Usage:  python studies/run_convergence.py [--quick] [--cases tg kov20 kov40 cfd]
Output: results/convergence.csv, figures fig01_tg_convergence,
        fig02_kovasznay_fields, fig01b_convergence_all_cases.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analysis.metrics import convergence_slope, error_map  # noqa: E402
from analysis.plots import (add_slope_triangle, error_panel, field_row, plot_series,  # noqa: E402
                            save_figure, set_style, style_of)
from core.pipeline import reconstruct  # noqa: E402
from studies.common import (BASE_SIGMA, BASE_WINDOW, CASE_LABELS, available_cases,  # noqa: E402
                            base_parser, evaluate, make_case, mean_std, seeds_for, write_csv)

CFD_STRIDES = (4, 2, 1)   # PIV spacing h = 0.2, 0.1, 0.05 D


def run(args) -> list[dict]:
    """Run the sweep and return all result rows."""
    ns = [16, 32, 64] if args.quick else args.ns
    seeds = seeds_for(args)
    rows: list[dict] = []
    for case in available_cases(args.cases):
        grid = [(s, None) for s in CFD_STRIDES] if case == "cfd" else [(1, n) for n in ns]
        for stride, n in grid:
            for condition in ("clean", "noisy"):
                window, sigma = (0.0, 0.0) if condition == "clean" else (BASE_WINDOW, args.sigma)
                for seed in ([0] if condition == "clean" else seeds):
                    fld = (make_case("cfd", window_h=window, sigma=sigma, seed=seed, stride=stride)
                           if case == "cfd" else make_case(case, n, window, sigma, seed))
                    omni_limit = args.omni_max_n if condition == "clean" else args.omni_max_n_noisy
                    for method in args.methods:
                        if method == "omni" and fld.nx > omni_limit:
                            continue
                        row = evaluate(fld, case, method=method)
                        row.update(case=case, condition=condition, method=method, seed=seed,
                                   sigma=sigma, window_h=window)
                        rows.append(row)
                print(f"  {case:6s} nx={fld.nx:4d} {condition} done")
    return rows


def optimum_h(rows: list[dict], case: str, method: str = "poisson") -> tuple[float, float]:
    """(h, mean eps_L2) minimising the noisy error for one case/method."""
    hs = sorted({float(r["h"]) for r in rows if r["case"] == case and r["condition"] == "noisy"
                 and r["method"] == method})
    means = [mean_std(rows, "eps_l2", case=case, condition="noisy", method=method, h=h)[0] for h in hs]
    if not hs:
        return float("nan"), float("nan")
    k = int(np.nanargmin(means))
    return hs[k], means[k]


def plot_convergence(rows: list[dict], cases: list[str], name: str, methods) -> None:
    """Log-log error vs h: clean (left) and noisy mean +/- std (right) per case."""
    fig, axes = plt.subplots(len(cases), 2, figsize=(7.2, 2.9 * len(cases)), squeeze=False)
    for r_i, case in enumerate(cases):
        for c_i, condition in enumerate(("clean", "noisy")):
            ax = axes[r_i, c_i]
            ref_line = None
            for method in methods:
                hs = sorted({float(r["h"]) for r in rows if r["case"] == case
                             and r["condition"] == condition and r["method"] == method})
                if not hs:
                    continue
                stats = [mean_std(rows, "eps_l2", case=case, condition=condition, method=method, h=h)
                         for h in hs]
                m = np.array([s[0] for s in stats])
                sd = np.array([s[1] for s in stats])
                slope = convergence_slope(hs, m)
                lab = None
                if condition == "clean":
                    lab = f"{style_of(method)['label']} (slope {slope:.2f})"
                plot_series(ax, hs, m, method, yerr=sd if condition == "noisy" else None, label=lab)
                if condition == "clean" and method == "poisson":
                    ref_line = (hs, m)
            if ref_line is not None:
                add_slope_triangle(ax, *ref_line, 2.0, "reference slope 2")
            if condition == "noisy":
                h_opt, e_opt = optimum_h(rows, case)
                if np.isfinite(h_opt):
                    ax.plot([h_opt], [e_opt], marker="*", ms=11, color="#0b0b0b", ls="",
                            label=f"optimum h = {h_opt:.3g}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("PIV grid spacing h (non-dimensional)")
            ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
            title = "clean (exact samples)" if condition == "clean" else \
                f"noisy (W = h, $\\sigma$ = {100 * BASE_SIGMA:g}% $U$), mean $\\pm$ std"
            ax.set_title(f"{CASE_LABELS[case]}: {title}")
            ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, name)


def plot_kovasznay_fields(n: int = 64) -> None:
    """Figure 2: exact vs reconstructed Kovasznay pressure + error map."""
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4))
    for r_i, case in enumerate(("kov20", "kov40")):
        fld = make_case(case, n, window_h=0.0, sigma=0.0)
        rec = reconstruct(fld, method="poisson")
        p_ref = rec.p_ref - rec.p_ref.mean()
        p = rec.p - rec.p.mean()
        vmax = max(np.abs(p_ref).max(), np.abs(p).max())
        field_row(axes[r_i, :2], fld.x, fld.y, [p_ref, p],
                  [f"{CASE_LABELS[case]}: exact $p$", f"reconstructed $p$ (Poisson, {n}$^2$)"],
                  clim=(-vmax, vmax), cbar_label=r"$p - \bar p$ ($\rho U^2$)")
        error_panel(axes[r_i, 2], fld.x, fld.y, error_map(p, p_ref), "error map $|p - p_{ref}|$",
                    label=r"$|p - p_{ref}|$ ($\rho U^2$)")
    save_figure(fig, "fig02_kovasznay_fields")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg", "kov20", "kov40", "cfd"))
    ap.add_argument("--ns", type=int, nargs="+", default=[16, 32, 64, 128, 256])
    ap.add_argument("--sigma", type=float, default=BASE_SIGMA)
    ap.add_argument("--methods", nargs="+", default=["poisson", "direct", "omni"])
    ap.add_argument("--omni-max-n", type=int, default=256, help="largest grid for omni (clean)")
    ap.add_argument("--omni-max-n-noisy", type=int, default=64, help="largest grid for omni (noisy)")
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "convergence.csv")
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    if "tg" in cases:
        plot_convergence(rows, ["tg"], "fig01_tg_convergence", args.methods)
    plot_convergence(rows, cases, "fig01b_convergence_all_cases", args.methods)
    plot_kovasznay_fields()
    for case in cases:
        for method in args.methods:
            hs = sorted({float(r["h"]) for r in rows if r["case"] == case
                         and r["condition"] == "clean" and r["method"] == method})
            errs = [mean_std(rows, "eps_l2", case=case, condition="clean", method=method, h=h)[0]
                    for h in hs]
            print(f"  slope {case:6s} {method:8s} clean: {convergence_slope(hs, errs):.2f}")
        print(f"  optimum h {case}: {optimum_h(rows, case)}")
    return rows


if __name__ == "__main__":
    main()
