# piv_pressure/studies/run_dt.py
"""Time-step study: Eulerian vs Lagrangian material derivative (Figure 7).

Truncation error of both schemes grows as dt^2, but the Eulerian central
difference in time has to follow structures that are *advected* past the
fixed node (du/dt ~ -U du/dx), whereas the Lagrangian scheme follows them.
Noise enters Du/Dt as ~ sigma / dt, so noisy data favour larger dt.

Leading truncation errors (see EXPLANATIONS.md): Eulerian ~ (dt^2/6) d3u/dt3,
which for a pattern carried at speed c scales like c^3; Lagrangian
~ (dt^2/6) (w.grad)^3 u, which scales like w^3 with w the velocity relative
to the pattern.  The Lagrangian scheme therefore wins when w < c.

Cases: ``tgc`` (Taylor-Green of amplitude A = 0.4 advected at c = 1 - an
analytic analogue of a wake vortex, w/c = 0.4), ``tg`` (stationary vortex:
du/dt = -2 nu u is tiny, so the Eulerian time step hardly matters - a
control case) and ``cfd`` (multiples k of the export step).

Usage:  python studies/run_dt.py [--quick] [--cases tgc tg cfd]
Output: results/dt.csv, figure fig07_dt.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402

from analysis.plots import plot_series, save_figure, set_style  # noqa: E402
from studies.common import (BASE_N, BASE_SIGMA, BASE_WINDOW, CASE_LABELS,  # noqa: E402
                            available_cases, base_parser, evaluate, get_cfd, make_case,
                            mean_std, seeds_for, write_csv)

SCHEMES = ("eulerian", "lagrangian")


def run(args) -> list[dict]:
    """Run the sweep and return all result rows."""
    rows: list[dict] = []
    for case in available_cases(args.cases):
        if case == "cfd":
            step = get_cfd().dt
            dts = [k * step for k in (1, 2, 3, 4)]
        else:
            dts = [0.02, 0.05, 0.2] if args.quick else args.dts
        for dt in dts:
            for condition in ("clean", "noisy"):
                window, sigma = (0.0, 0.0) if condition == "clean" else (BASE_WINDOW, BASE_SIGMA)
                for seed in ([0] if condition == "clean" else seeds_for(args)):
                    fld = make_case(case, BASE_N, window, sigma, seed, dt=dt)
                    for scheme in SCHEMES:
                        row = evaluate(fld, case, method="poisson", material_derivative=scheme)
                        row.update(case=case, dt=float(fld.dt), condition=condition, seed=seed,
                                   scheme=scheme, sigma=sigma)
                        rows.append(row)
            print(f"  {case:6s} dt={dt:.3f} done")
    return rows


def plot(rows: list[dict]) -> None:
    """Figure 7: eps_L2 vs dt; clean (top) and noisy (bottom) per case."""
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    fig, axes = plt.subplots(2, len(cases), figsize=(3.6 * len(cases), 6.0), squeeze=False)
    for c_i, case in enumerate(cases):
        for r_i, condition in enumerate(("clean", "noisy")):
            ax = axes[r_i, c_i]
            dts = sorted({float(r["dt"]) for r in rows if r["case"] == case})
            for scheme in SCHEMES:
                st = [mean_std(rows, "eps_l2", case=case, condition=condition, scheme=scheme, dt=d)
                      for d in dts]
                plot_series(ax, dts, [s[0] for s in st], scheme,
                            yerr=[s[1] for s in st] if condition == "noisy" else None)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"time step $\Delta t$ ($D/U$ or $L/U$)")
            ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
            label = "clean" if condition == "clean" else r"noisy ($\sigma$ = 1%, W = h)"
            ax.set_title(f"{CASE_LABELS[case]}\n{label}")
            ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, "fig07_dt")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tgc", "tg", "cfd"))
    ap.add_argument("--dts", type=float, nargs="+", default=[0.01, 0.02, 0.05, 0.1, 0.2])
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "dt.csv")
    plot(rows)
    for case in sorted({r["case"] for r in rows}):
        dmax = max(float(r["dt"]) for r in rows if r["case"] == case)
        e = {s: mean_std(rows, "eps_l2", case=case, condition="clean", scheme=s, dt=dmax)[0]
             for s in SCHEMES}
        print(f"  {case:6s} clean at dt={dmax:.2f}: " + ", ".join(f"{k} {v:.3e}" for k, v in e.items()))
    return rows


if __name__ == "__main__":
    main()
