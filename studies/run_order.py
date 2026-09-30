# piv_pressure/studies/run_order.py
"""Grid convergence, lecture style: first- vs second-order scheme (Figure 16).

For each case the full reconstruction (G -> Poisson) is repeated with
``order=1`` (forward differences, first-order Neumann BCs) and ``order=2``
(the default scheme) on grids from 8^2 to 256^2 (refinement factor ~sqrt 2):

* error (%) = 100 x relative L2 error, plotted on log axes against dx with
  the dx axis reversed (finer grids to the right), slope triangles for 1
  and 2 and the fitted slopes in the legend;
* the 1 % target line with guide lines down to the dx at which each scheme
  first drops below 1 % (log-log interpolation between the two bracketing
  grids; if a scheme never reaches 1 % the crossing is extrapolated with its
  asymptotic slope and marked as such);
* an overlay for noisy data (W = h, sigma = 1 %, mean +/- std over seeds)
  showing that the error stops falling (or rises) at small dx.

Slopes are fitted over the finest ``--fit-points`` grids (the asymptotic
range); the fit over all grids is also written to the summary CSV.

Usage:  python studies/run_order.py [--quick] [--cases tg kov20]
Output: results/order_convergence.csv, results/order_summary.csv,
        figures fig16_grid_convergence_tg, fig16b_grid_convergence_kov20.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analysis.metrics import convergence_slope  # noqa: E402
from analysis.plots import INK_2, draw_slope_triangle, save_figure, set_style, style_of, threshold_guides  # noqa: E402
from studies.common import (BASE_SIGMA, BASE_WINDOW, CASE_LABELS, available_cases,  # noqa: E402
                            base_parser, evaluate, make_case, mean_std, seeds_for, write_csv)

ORDERS = (1, 2)
TARGET_PCT = 1.0
FIG_NAMES = {"tg": "fig16_grid_convergence_tg", "kov20": "fig16b_grid_convergence_kov20",
             "kov40": "fig16c_grid_convergence_kov40"}


def run(args) -> list[dict]:
    """All reconstructions (one row each)."""
    ns = [8, 11, 16, 23, 32, 45, 64] if args.quick else args.ns
    rows: list[dict] = []
    for case in available_cases([c for c in args.cases if c != "cfd"]):
        for n in ns:
            for condition in ("clean", "noisy"):
                window, sigma = (0.0, 0.0) if condition == "clean" else (BASE_WINDOW, BASE_SIGMA)
                for seed in ([0] if condition == "clean" else seeds_for(args)):
                    fld = make_case(case, n, window, sigma, seed)
                    for order in ORDERS:
                        row = evaluate(fld, case, method="poisson", order=order)
                        row.update(case=case, n=n, dx=fld.dx, order=order, condition=condition,
                                   seed=seed, sigma=sigma)
                        rows.append(row)
            print(f"  {case:6s} n={n:4d} done")
    return rows


def crossing(dx: list[float], err: list[float], slope: float, target: float = TARGET_PCT) -> dict:
    """First dx (coarse -> fine) with err < target, and the interpolated crossing.

    ``dx`` is sorted from coarse to fine.  Returns dx_first_below (a tested
    grid, NaN if none), dx_cross (log-log interpolation, or extrapolation with
    ``slope`` beyond the finest grid) and whether it was extrapolated.
    """
    dx, err = np.asarray(dx, float), np.asarray(err, float)
    below = np.nonzero(err < target)[0]
    if below.size:
        k = int(below[0])
        if k == 0:
            return {"dx_first_below": dx[0], "dx_cross": np.nan, "cross_extrapolated": False,
                    "note": "below target already on the coarsest grid"}
        t = (np.log(target) - np.log(err[k - 1])) / (np.log(err[k]) - np.log(err[k - 1]))
        return {"dx_first_below": dx[k],
                "dx_cross": float(np.exp(np.log(dx[k - 1]) + t * (np.log(dx[k]) - np.log(dx[k - 1])))),
                "cross_extrapolated": False, "note": ""}
    if not np.isfinite(slope) or slope <= 0:
        return {"dx_first_below": np.nan, "dx_cross": np.nan, "cross_extrapolated": True,
                "note": "not reached; no positive slope to extrapolate"}
    return {"dx_first_below": np.nan,
            "dx_cross": float(dx[-1] * (target / err[-1]) ** (1.0 / slope)),
            "cross_extrapolated": True, "note": "not reached within the tested grids (extrapolated)"}


def summarise(rows: list[dict], fit_points: int) -> list[dict]:
    """Slopes, 1 % crossings and noisy optimum per case and order."""
    out = []
    for case in sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index):
        for order in ORDERS:
            dxs = sorted({float(r["dx"]) for r in rows if r["case"] == case}, reverse=True)
            clean = [mean_std(rows, "err_pct", case=case, order=order, condition="clean", dx=d)[0] for d in dxs]
            noisy = [mean_std(rows, "err_pct", case=case, order=order, condition="noisy", dx=d)[0] for d in dxs]
            slope_fine = convergence_slope(dxs[-fit_points:], clean[-fit_points:])
            k_opt = int(np.nanargmin(noisy))
            out.append({
                "case": case, "order": order, "slope_all": convergence_slope(dxs, clean),
                "slope_fine": slope_fine, "fit_points": min(fit_points, len(dxs)),
                **crossing(dxs, clean, slope_fine),
                "noisy_min_err_pct": noisy[k_opt], "noisy_optimum_dx": dxs[k_opt],
                "finest_dx": dxs[-1], "finest_err_pct": clean[-1],
            })
    return out


def plot_case(rows: list[dict], summary: list[dict], case: str) -> None:
    """Figure 16 for one case."""
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    dxs = sorted({float(r["dx"]) for r in rows if r["case"] == case}, reverse=True)
    for order in ORDERS:
        key = f"order{order}"
        st = style_of(key)
        s = next(x for x in summary if x["case"] == case and x["order"] == order)
        clean = [mean_std(rows, "err_pct", case=case, order=order, condition="clean", dx=d)[0] for d in dxs]
        stats = [mean_std(rows, "err_pct", case=case, order=order, condition="noisy", dx=d) for d in dxs]
        m = np.array([q[0] for q in stats])
        sd = np.array([q[1] for q in stats])
        ax.plot(dxs, clean, color=st["color"], marker=st["marker"], lw=1.6,
                label=f"order {order}, clean: fitted slope {s['slope_fine']:.2f}")
        ax.plot(dxs, m, color=st["color"], marker=st["marker"], mfc="white", ls="--", lw=1.1,
                label=rf"order {order}, noisy ($\sigma$ = 1 %, W = h), mean $\pm$ std")
        ax.fill_between(dxs, np.maximum(m - sd, 1e-6), m + sd, color=st["color"], alpha=0.12, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(TARGET_PCT, color=INK_2, ls="--", lw=1.0, label="1 % target")
    # slope triangles just below the clean curves, in the middle of the range
    for order in ORDERS:
        clean = [mean_std(rows, "err_pct", case=case, order=order, condition="clean", dx=d)[0] for d in dxs]
        i = len(dxs) // 2
        draw_slope_triangle(ax, dxs[i] / 2.0, clean[i] / 2.0 ** order * 0.35, order, factor=2.0,
                            color=style_of(f"order{order}")["color"])
    # guide lines down to the 1 % crossing
    for order in ORDERS:
        s = next(x for x in summary if x["case"] == case and x["order"] == order)
        if np.isfinite(s["dx_cross"]):
            lab = (f"order {order}: 1 % reached at dx = {s['dx_cross']:.3g}"
                   + (" (extrapolated)" if s["cross_extrapolated"] else ""))
            if s["cross_extrapolated"]:
                ax.plot([s["finest_dx"], s["dx_cross"]], [s["finest_err_pct"], TARGET_PCT],
                        color=style_of(f"order{order}")["color"], ls=(0, (1, 3)), lw=1.2, marker="")
                lo, hi = ax.get_xlim()                          # not yet reversed: lo = finest
                ax.set_xlim(min(lo, s["dx_cross"] / 1.3), hi)
            threshold_guides(ax, s["dx_cross"], TARGET_PCT, style_of(f"order{order}")["color"], lab,
                             s["cross_extrapolated"])
    ax.set_xlim(max(ax.get_xlim()), min(ax.get_xlim()))      # reverse: coarse left, fine right
    ax.set_xlabel(r"grid spacing $\Delta x$ (non-dimensional; finer $\rightarrow$)")
    ax.set_ylabel("pressure error (%)  = 100 x relative L2")
    ax.set_title(f"{CASE_LABELS[case]}: grid convergence of the reconstructed pressure")
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    save_figure(fig, FIG_NAMES.get(case, f"fig16_grid_convergence_{case}"))


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[0], ("tg", "kov20"))
    ap.add_argument("--ns", type=int, nargs="+", default=[8, 11, 16, 23, 32, 45, 64, 91, 128, 181, 256])
    ap.add_argument("--fit-points", type=int, default=5, help="finest grids used for the fitted slope")
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "order_convergence.csv")
    summary = summarise(rows, args.fit_points)
    write_csv(summary, args.outdir / "order_summary.csv")
    for case in sorted({r["case"] for r in rows}):
        plot_case(rows, summary, case)
    for s in summary:
        print(f"  {s['case']:6s} order {s['order']}: slope {s['slope_fine']:.2f} (all {s['slope_all']:.2f}), "
              f"1 % at dx = {s['dx_cross']:.3g}{' (extrapolated)' if s['cross_extrapolated'] else ''}")
    return rows


if __name__ == "__main__":
    main()
