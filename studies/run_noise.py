# piv_pressure/studies/run_noise.py
"""Noise study: error vs velocity-noise level sigma (Figure 3).

Baseline PIV (n = 64, W = 1 h, dt = 0.05); sigma = 0, 0.5, 1, 2, 5, 10 % of U;
>= 10 noise realisations per level; Poisson vs direct vs omni integration.
Also records the wall-clock time of every reconstruction.

Usage:  python studies/run_noise.py [--quick] [--cases tg kov20 cfd]
Output: results/noise.csv, figure fig03_noise.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analysis.metrics import convergence_slope  # noqa: E402
from analysis.plots import plot_series, save_figure, set_style  # noqa: E402
from studies.common import (BASE_N, BASE_WINDOW, CASE_LABELS, available_cases,  # noqa: E402
                            base_parser, evaluate, make_case, mean_std, seeds_for, write_csv)


def run(args) -> list[dict]:
    """Run the sweep and return all result rows."""
    sigmas = [0.0, 0.01, 0.05, 0.10] if args.quick else args.sigmas
    rows: list[dict] = []
    for case in available_cases(args.cases):
        for sigma in sigmas:
            for seed in ([0] if sigma == 0.0 else seeds_for(args)):
                fld = make_case(case, BASE_N, BASE_WINDOW, sigma, seed)
                for method in args.methods:
                    row = evaluate(fld, case, method=method)
                    row.update(case=case, sigma=sigma, seed=seed, method=method)
                    rows.append(row)
            print(f"  {case:6s} sigma={sigma:.3f} done")
    return rows


def noise_slope(rows: list[dict], case: str, method: str) -> float:
    """log-log slope of mean eps_L2 vs sigma over sigma > 0 (1 = linear growth)."""
    sig = sorted({float(r["sigma"]) for r in rows if r["case"] == case and float(r["sigma"]) > 0})
    m = [mean_std(rows, "eps_l2", case=case, method=method, sigma=s)[0] for s in sig]
    return convergence_slope(sig, m)


def plot(rows: list[dict], methods: list[str]) -> None:
    """Figure 3: mean +/- std of eps_L2 vs sigma, one panel per case."""
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    fig, axes = plt.subplots(1, len(cases), figsize=(3.6 * len(cases), 3.2), squeeze=False)
    for ax, case in zip(axes[0], cases):
        sig = sorted({float(r["sigma"]) for r in rows if r["case"] == case})
        for method in methods:
            st = [mean_std(rows, "eps_l2", case=case, method=method, sigma=s) for s in sig]
            m = np.array([s[0] for s in st])
            sd = np.array([s[1] for s in st])
            plot_series(ax, 100 * np.array(sig), m, method, yerr=sd)
        ax.set_xlabel(r"noise level $\sigma$ (% of $U_{ref}$)")
        ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
        ax.set_title(CASE_LABELS[case])
        ax.set_ylim(bottom=0)
    axes[0, 0].legend(loc="upper left")
    fig.suptitle(r"Error vs noise (n = 64, W = h, $\Delta t$ = 0.05), mean $\pm$ std over seeds", y=1.02)
    fig.tight_layout()
    save_figure(fig, "fig03_noise")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg", "kov20", "cfd"))
    ap.add_argument("--sigmas", type=float, nargs="+", default=[0.0, 0.005, 0.01, 0.02, 0.05, 0.10])
    ap.add_argument("--methods", nargs="+", default=["poisson", "direct", "omni"])
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "noise.csv")
    plot(rows, args.methods)
    for case in sorted({r["case"] for r in rows}):
        for method in args.methods:
            t = mean_std(rows, "time_s", case=case, method=method)[0]
            print(f"  {case:6s} {method:8s} noise slope {noise_slope(rows, case, method):.2f}"
                  f"  mean time {t * 1e3:.1f} ms")
    return rows


if __name__ == "__main__":
    main()
