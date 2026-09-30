# piv_pressure/studies/run_window.py
"""Interrogation-window study: effect of window averaging W (Figure 9a).

A box average of width W attenuates the velocity peaks around a vortex, and
since pressure depends quadratically on velocity, the suction peak in the
vortex core is under-predicted.  W = 0 means ideal point sampling
(reference); W = 1, 2, 4 h are the spec values.

Cases: ``tg_core`` (Taylor-Green window centred on the core at (pi, pi)),
``kov20`` (no vortex core - eps only) and ``cfd`` (PIV spacing 0.1 D, the
same 2-node margin cropped for every W so the domain is identical).

Usage:  python studies/run_window.py [--quick] [--cases tg_core kov20 cfd]
Output: results/window.csv, figure fig09a_window.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402

from analysis.plots import plot_series, save_figure, set_style  # noqa: E402
from studies.common import (BASE_N, BASE_SIGMA, CASE_LABELS, available_cases,  # noqa: E402
                            base_parser, cfd_piv, evaluate, has_core, make_case, mean_std,
                            seeds_for, write_csv)

WINDOWS = (0.0, 1.0, 2.0, 4.0)
CFD_STRIDE = 2           # PIV spacing 0.1 D so that W = 1 h already averages 2 export points
CFD_CROP = 2             # = max W / 2 in PIV spacings: filter never reaches past the data


def build(case: str, window_h: float, sigma: float, seed: int):
    """Synthetic PIV field for this study."""
    if case == "cfd":
        return cfd_piv(stride=CFD_STRIDE, window_h=window_h, sigma=sigma, seed=seed, crop=CFD_CROP)
    return make_case(case, BASE_N, window_h, sigma, seed)


def run(args) -> list[dict]:
    """Run the sweep and return all result rows."""
    rows: list[dict] = []
    for case in available_cases(args.cases):
        for window_h in WINDOWS:
            for condition in ("clean", "noisy"):
                sigma = 0.0 if condition == "clean" else BASE_SIGMA
                for seed in ([0] if condition == "clean" else seeds_for(args)):
                    fld = build(case, window_h, sigma, seed)
                    row = evaluate(fld, case, with_core=has_core(case), method="poisson")
                    row.update(case=case, window_h=window_h, condition=condition, seed=seed,
                               sigma=sigma)
                    rows.append(row)
            print(f"  {case:8s} W={window_h:g}h done")
    return rows


def plot_window(ax_core: plt.Axes, ax_eps: plt.Axes | None, rows: list[dict]) -> None:
    """Core-suction ratio (left) and eps_L2 (right, skipped if None) vs W/h; clean data."""
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    for case in cases:
        ws = sorted({float(r["window_h"]) for r in rows if r["case"] == case})
        eps = [mean_std(rows, "eps_l2", case=case, condition="clean", window_h=w)[0] for w in ws]
        if ax_eps is not None:
            plot_series(ax_eps, ws, eps, f"case:{case}")
        if has_core(case):
            ratio = [mean_std(rows, "core_ratio", case=case, condition="clean", window_h=w)[0]
                     for w in ws]
            plot_series(ax_core, ws, ratio, f"case:{case}")
    ax_core.axhline(1.0, color="#52514e", lw=1.0, ls="--", label="exact core pressure")
    ax_core.set_xlabel("interrogation window W / h (-)")
    ax_core.set_ylabel(r"$p_{core,rec} / p_{core,ref}$ (-)")
    ax_core.set_title("Vortex-core suction vs window size (clean)")
    ax_core.legend(loc="best")
    if ax_eps is None:
        return
    ax_eps.set_xlabel("interrogation window W / h (-)")
    ax_eps.set_ylabel(r"$\epsilon_{L2}$ (-)")
    ax_eps.set_yscale("log")
    ax_eps.set_title("Field error vs window size (clean)")
    ax_eps.legend(loc="best")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg_core", "kov20", "cfd"))
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "window.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.3))
    plot_window(a1, a2, rows)
    fig.tight_layout()
    save_figure(fig, "fig09a_window")
    for case in sorted({r["case"] for r in rows}):
        for w in WINDOWS:
            e = mean_std(rows, "eps_l2", case=case, condition="clean", window_h=w)[0]
            c = mean_std(rows, "core_ratio", case=case, condition="clean", window_h=w)[0]
            print(f"  {case:8s} W={w:g}h eps {e:.3e} core ratio {c:.3f}")
    return rows


if __name__ == "__main__":
    main()
