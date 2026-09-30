# piv_pressure/studies/run_method.py
"""Integration-method study (Figure 8; combined with the BC study).

Compares direct integration (path 1, path 2, their average), the Poisson
solver (div-G and velocity-gradient source forms) and omni-directional
integration on clean and noisy data, and records the wall-clock time and
the path-dependence indicators.

Usage:  python studies/run_method.py [--quick] [--cases tg kov20 cfd]
Output: results/method.csv, figures fig08a_method and (if results/bc.csv
        exists) the combined fig08_bc_method.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402

from analysis.plots import save_figure, set_style, style_of  # noqa: E402
from studies.common import (BASE_N, BASE_SIGMA, BASE_WINDOW, CASE_LABELS,  # noqa: E402
                            available_cases, base_parser, evaluate, make_case, mean_std,
                            read_csv, seeds_for, write_csv)
from studies.run_bc import plot_bc  # noqa: E402

VARIANTS: dict[str, dict] = {
    "direct_path1": dict(method="direct", direct_path="path1"),
    "direct_path2": dict(method="direct", direct_path="path2"),
    "direct": dict(method="direct", direct_path="average"),
    "poisson": dict(method="poisson"),
    "poisson_velocity": dict(method="poisson", source="velocity"),
    "omni": dict(method="omni"),
}


def run(args) -> list[dict]:
    """Run the comparison and return all result rows."""
    rows: list[dict] = []
    for case in available_cases(args.cases):
        for condition in ("clean", "noisy"):
            window, sigma = (0.0, 0.0) if condition == "clean" else (BASE_WINDOW, BASE_SIGMA)
            for seed in ([0] if condition == "clean" else seeds_for(args)):
                fld = make_case(case, BASE_N, window, sigma, seed)
                for variant, kw in VARIANTS.items():
                    row = evaluate(fld, case, **kw)
                    row.update(case=case, condition=condition, seed=seed, variant=variant)
                    rows.append(row)
        print(f"  {case:6s} done")
    return rows


def plot_methods(ax: plt.Axes, rows: list[dict]) -> None:
    """Dot plot of eps_L2 per method (x = case; clean hollow, noisy filled)."""
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    names = list(VARIANTS)
    width = 0.8 / len(names)
    for v_i, variant in enumerate(names):
        st = style_of(variant)
        for c_i, case in enumerate(cases):
            xpos = c_i + (v_i - (len(names) - 1) / 2) * width
            m_c, _, _ = mean_std(rows, "eps_l2", case=case, condition="clean", variant=variant)
            m_n, s_n, _ = mean_std(rows, "eps_l2", case=case, condition="noisy", variant=variant)
            ax.plot([xpos], [m_c], marker=st["marker"], mfc="white", mec=st["color"], ls="", ms=6)
            ax.errorbar([xpos], [m_n], yerr=[s_n], marker=st["marker"], color=st["color"], ls="",
                        ms=6, capsize=2, label=st["label"] if c_i == 0 else None)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([CASE_LABELS[c] for c in cases], fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
    ax.set_title("Integration methods: hollow = clean, filled = noisy mean $\\pm$ std")
    ax.legend(loc="best", fontsize=7, ncol=2)


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg", "kov20", "cfd"))
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "method.csv")
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    plot_methods(ax, rows)
    fig.tight_layout()
    save_figure(fig, "fig08a_method")
    bc_csv = args.outdir / "bc.csv"
    if bc_csv.is_file():
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.4))
        plot_methods(a1, rows)
        plot_bc(a2, read_csv(bc_csv))
        a1.set_title("(a) " + a1.get_title())
        a2.set_title("(b) " + a2.get_title())
        fig.tight_layout()
        save_figure(fig, "fig08_bc_method")
    else:
        print("[info] results/bc.csv not found - run studies/run_bc.py for the combined figure 8")
    for case in sorted({r["case"] for r in rows}):
        for variant in VARIANTS:
            mn = mean_std(rows, "eps_l2", case=case, condition="noisy", variant=variant)[0]
            t = mean_std(rows, "time_s", case=case, variant=variant)[0]
            print(f"  {case:6s} {variant:17s} noisy {mn:.3e}  time {t * 1e3:7.1f} ms")
    return rows


if __name__ == "__main__":
    main()
