# piv_pressure/studies/run_bc.py
"""Boundary-condition study for the Poisson method (Figure 8, panel b).

Variants
--------
``neumann``      dp/dn = n.G on all edges, zero-mean constraint (baseline)
``neumann_node`` all-Neumann, one node pinned to p = 0
``bernoulli_1``  Dirichlet p = p_inf + rho (U_inf^2 - |u|^2)/2 on ONE edge,
                 Neumann elsewhere
``exact_1``      Dirichlet from the true pressure on the same edge (an oracle,
                 e.g. pressure taps) - what a *correct* Dirichlet edge can do
``bernoulli_2``, ``exact_2``  the same on TWO opposite edges

With a single Dirichlet edge the Neumann-data error is "anchored" on one
side and accumulates towards the opposite side, so even exact Dirichlet data
can raise the error; two opposite edges bracket the domain.

Bernoulli is only valid where the flow is steady and irrotational.  Neither
analytic case has such an edge (both are rotational everywhere), so the
Bernoulli edge is *expected* to add error there; for the cylinder wake the
top/bottom edges (|y| = 2 D) are nearly irrotational but unsteady.

Usage:  python studies/run_bc.py [--quick] [--cases tg kov20 cfd]
Output: results/bc.csv, figure fig08b_bc.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402

from analysis.plots import save_figure, set_style, style_of  # noqa: E402
from studies.common import (BASE_N, BASE_SIGMA, BASE_WINDOW, CASE_LABELS,  # noqa: E402
                            available_cases, base_parser, evaluate, make_case, mean_std,
                            seeds_for, write_csv)

VARIANTS = ("neumann", "neumann_node", "bernoulli_1", "exact_1", "bernoulli_2", "exact_2")
#: Dirichlet edges per case: (one edge, two opposite edges).  For the wake the
#: top/bottom edges (|y| = 2 D) are the most free-stream-like; for Kovasznay
#: the outflow edge (x = 1.5) is closest to uniform flow.
EDGES = {
    "tg": (("top",), ("top", "bottom")),
    "kov20": (("right",), ("left", "right")),
    "cfd": (("top",), ("top", "bottom")),
}


def variant_kwargs(variant: str, case: str) -> dict:
    """reconstruct() keyword arguments for a BC variant."""
    one, two = EDGES.get(case, (("top",), ("top", "bottom")))
    return {
        "neumann": dict(bc="neumann", pin="mean"),
        "neumann_node": dict(bc="neumann", pin="node"),
        "bernoulli_1": dict(bc="bernoulli", dirichlet_edges=one),
        "exact_1": dict(bc="exact", dirichlet_edges=one),
        "bernoulli_2": dict(bc="bernoulli", dirichlet_edges=two),
        "exact_2": dict(bc="exact", dirichlet_edges=two),
    }[variant]


def run(args) -> list[dict]:
    """Run the comparison and return all result rows."""
    rows: list[dict] = []
    for case in available_cases(args.cases):
        for condition in ("clean", "noisy"):
            window, sigma = (0.0, 0.0) if condition == "clean" else (BASE_WINDOW, BASE_SIGMA)
            for seed in ([0] if condition == "clean" else seeds_for(args)):
                fld = make_case(case, BASE_N, window, sigma, seed)
                for variant in VARIANTS:
                    row = evaluate(fld, case, method="poisson", **variant_kwargs(variant, case))
                    kw = variant_kwargs(variant, case)
                    row.update(case=case, condition=condition, seed=seed, variant=variant,
                               edges="+".join(kw.get("dirichlet_edges", ())) or "none")
                    rows.append(row)
        print(f"  {case:6s} done")
    return rows


def plot_bc(ax: plt.Axes, rows: list[dict]) -> None:
    """Dot plot of eps_L2 per BC variant (x = case; clean hollow, noisy filled)."""
    cases = sorted({r["case"] for r in rows}, key=list(CASE_LABELS).index)
    width = 0.8 / len(VARIANTS)
    for v_i, variant in enumerate(VARIANTS):
        st = style_of(variant)
        for c_i, case in enumerate(cases):
            xpos = c_i + (v_i - (len(VARIANTS) - 1) / 2) * width
            m_c, _, _ = mean_std(rows, "eps_l2", case=case, condition="clean", variant=variant)
            m_n, s_n, _ = mean_std(rows, "eps_l2", case=case, condition="noisy", variant=variant)
            ax.plot([xpos], [m_c], marker=st["marker"], mfc="white", mec=st["color"], ls="", ms=6)
            ax.errorbar([xpos], [m_n], yerr=[s_n], marker=st["marker"], color=st["color"], ls="",
                        ms=6, capsize=2, label=st["label"] if c_i == 0 else None)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([CASE_LABELS[c] for c in cases], fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
    ax.set_title("Boundary conditions (Poisson): hollow = clean, filled = noisy mean $\\pm$ std")
    ax.legend(loc="best", fontsize=7)


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[1], ("tg", "kov20", "cfd"))
    args = ap.parse_args(argv)
    set_style()
    rows = run(args)
    write_csv(rows, args.outdir / "bc.csv")
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    plot_bc(ax, rows)
    fig.tight_layout()
    save_figure(fig, "fig08b_bc")
    for case in sorted({r["case"] for r in rows}):
        for variant in VARIANTS:
            mc = mean_std(rows, "eps_l2", case=case, condition="clean", variant=variant)[0]
            mn = mean_std(rows, "eps_l2", case=case, condition="noisy", variant=variant)[0]
            print(f"  {case:6s} {variant:13s} clean {mc:.3e} noisy {mn:.3e}")
    return rows


if __name__ == "__main__":
    main()
