# piv_pressure/studies/run_denoise.py
"""Denoising study: none vs Gaussian vs 2-D Savitzky-Golay (Figure 9).

Filtering the velocity before differentiation removes noise (good at high
sigma) but also smooths the true velocity peaks (bias, bad at low sigma).
The vortex-core suction is tracked because it is the quantity most hurt by
smoothing.

Usage:  python studies/run_denoise.py [--quick] [--cases tg_core kov20 cfd]
Output: results/denoise.csv, figures fig09b_denoise and (if
        results/window.csv exists) the combined fig09_window_denoise.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analysis.plots import plot_series, save_figure, set_style  # noqa: E402
from piv.denoise import denoise_field  # noqa: E402
from studies.common import (BASE_N, BASE_WINDOW, CASE_LABELS, available_cases,  # noqa: E402
                            base_parser, evaluate, has_core, make_case, mean_std, read_csv,
                            seeds_for, write_csv)
from studies.run_window import plot_window  # noqa: E402

DENOISERS = ("none", "gaussian", "savgol")


def run(args) -> list[dict]:
    """Run the sweep and return all result rows."""
    sigmas = [0.0, 0.01, 0.05, 0.10] if args.quick else args.sigmas
    sizes = {"none": 0, "gaussian": args.gauss_size, "savgol": args.sg_size}
    rows: list[dict] = []
    for case in available_cases(args.cases):
        for sigma in sigmas:
            for seed in ([0] if sigma == 0.0 else seeds_for(args)):
                raw = make_case(case, BASE_N, BASE_WINDOW, sigma, seed)
                for method in DENOISERS:
                    fld = denoise_field(raw, method, sizes[method], order=args.sg_order)
                    row = evaluate(fld, case, with_core=has_core(case), method="poisson")
                    row.update(case=case, sigma=sigma, seed=seed, denoise=method,
                               denoise_size=sizes[method])
                    rows.append(row)
            print(f"  {case:8s} sigma={sigma:.3f} done")
    return rows


def plot_denoise(ax_eps: plt.Axes, ax_core: plt.Axes | None, rows: list[dict], case: str) -> None:
    """eps_L2 (and core ratio) vs sigma for each denoiser, one case."""
    sig = sorted({float(r["sigma"]) for r in rows if r["case"] == case})
    for method in DENOISERS:
        st = [mean_std(rows, "eps_l2", case=case, denoise=method, sigma=s) for s in sig]
        plot_series(ax_eps, 100 * np.array(sig), [s[0] for s in st], method, yerr=[s[1] for s in st])
        if ax_core is not None and has_core(case):
            st = [mean_std(rows, "core_ratio", case=case, denoise=method, sigma=s) for s in sig]
            plot_series(ax_core, 100 * np.array(sig), [s[0] for s in st], method,
                        yerr=[s[1] for s in st])
    ax_eps.set_xlabel(r"noise level $\sigma$ (% of $U_{ref}$)")
    ax_eps.set_ylabel(r"$\epsilon_{L2}$ (-)")
    ax_eps.set_yscale("log")
    ax_eps.set_title(f"{CASE_LABELS[case]}: field error")
    ax_eps.legend(loc="best")
    if ax_core is not None:
        ax_core.axhline(1.0, color="#52514e", lw=1.0, ls="--", label="exact core pressure")
        ax_core.set_xlabel(r"noise level $\sigma$ (% of $U_{ref}$)")
        ax_core.set_ylabel(r"$p_{core,rec} / p_{core,ref}$ (-)")
        ax_core.set_title(f"{CASE_LABELS[case]}: core suction")
        ax_core.legend(loc="best")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg_core", "kov20", "cfd"))
    ap.add_argument("--sigmas", type=float, nargs="+", default=[0.0, 0.005, 0.01, 0.02, 0.05, 0.10])
    ap.add_argument("--gauss-size", type=float, default=1.0, help="Gaussian std (grid points)")
    ap.add_argument("--sg-size", type=int, default=5, help="Savitzky-Golay window (odd)")
    ap.add_argument("--sg-order", type=int, default=2, help="Savitzky-Golay polynomial order")
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "denoise.csv")
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    fig, axes = plt.subplots(2, len(cases), figsize=(3.7 * len(cases), 6.2), squeeze=False)
    for c_i, case in enumerate(cases):
        plot_denoise(axes[0, c_i], axes[1, c_i] if has_core(case) else None, rows, case)
        if not has_core(case):
            axes[1, c_i].axis("off")
    fig.tight_layout()
    save_figure(fig, "fig09b_denoise")
    win_csv = args.outdir / "window.csv"
    core_case = next((c for c in cases if has_core(c)), None)
    if win_csv.is_file() and core_case is not None:
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
        plot_window(axes[0], None, read_csv(win_csv))
        plot_denoise(axes[2], axes[1], rows, core_case)
        axes[0].set_title("(a) window averaging: core suction")
        axes[1].set_title(f"(b) denoising: core suction ({CASE_LABELS[core_case]})")
        axes[2].set_title(f"(c) denoising: field error ({CASE_LABELS[core_case]})")
        fig.tight_layout()
        save_figure(fig, "fig09_window_denoise")
    for case in cases:
        for s in sorted({float(r["sigma"]) for r in rows}):
            e = {m: mean_std(rows, "eps_l2", case=case, denoise=m, sigma=s)[0] for m in DENOISERS}
            print(f"  {case:8s} sigma={s:.3f} " + " ".join(f"{k} {v:.3e}" for k, v in e.items()))
    return rows


if __name__ == "__main__":
    main()
