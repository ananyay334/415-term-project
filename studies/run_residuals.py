# piv_pressure/studies/run_residuals.py
"""Iteration errors of the Poisson solve (ME-415 Lecture 11; Figures 11-15).

All runs solve the *same* discrete system (Taylor-Green, clean data,
Eulerian G, all-Neumann BCs with the zero-mean constraint, second order) so
only the solution method differs.

Figure 11  scaled residual L2_n(r)/S vs iteration for Jacobi, classic Jacobi,
           Gauss-Seidel, SOR and CG, with the tolerance line.
Figure 12  (a) one solver (SOR): L1_n, L2_n, Linf and scaled residual;
           (b) unnormalised vs normalised norms of the initial residual on
           grids 32^2-256^2 - why normalisation is needed.
Figure 13  Jacobi: scaled residual, change norm ||p^{k+1}-p^k|| and the true
           error ||p^k - p*|| - a small change per iteration is not a small
           residual (or error) when convergence is slow.
Figure 14  iterations to reach tol vs N = n^2 with fitted growth exponents.
Figure 15  SOR: iterations vs omega (omega_opt marked) and divergence for
           omega >= 2 (capped run).

Usage:  python studies/run_residuals.py [--quick] [--n 64] [--tol 1e-8]
                                        [--slow-max-n 128]
Output: results/residual_history.csv, residual_norms.csv, norm_scaling.csv,
        residual_vs_change.csv, iterations_vs_grid.csv, sor_omega.csv,
        sor_divergence.csv; figures fig11-fig15.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import LogLocator, ScalarFormatter  # noqa: E402

from analysis.metrics import convergence_slope  # noqa: E402
from analysis.norms import l1, l1_normalised, l2, l2_normalised  # noqa: E402
from analysis.plots import INK_2, plot_history, plot_series, save_figure, set_style, style_of  # noqa: E402
from core import iterative  # noqa: E402
from core.material_derivative import pressure_gradient  # noqa: E402
from core.poisson import linear_system, poisson_system, solve_poisson, source_div_G  # noqa: E402
from fields.analytic import taylor_green  # noqa: E402
from studies.common import base_parser, write_csv  # noqa: E402

SOLVERS = ("jacobi", "gs", "sor", "cg")


def baseline_system(n: int, order: int = 2):
    """(LinearSystem, dx, dy, direct solution p*) for clean Taylor-Green on n x n."""
    fld = taylor_green(n)
    g = pressure_gradient(fld, order=order)
    f = source_div_G(g.Gx, g.Gy, fld.dx, fld.dy, order)
    A, b, w, fixed, p_fixed, _ = poisson_system(f, g.Gx, g.Gy, fld.dx, fld.dy, order=order)
    system, _ = linear_system(A, b, w, fixed, p_fixed, (n, n))
    p_direct, _ = solve_poisson(f, g.Gx, g.Gy, fld.dx, fld.dy, order=order)
    return system, fld.dx, fld.dy, p_direct.ravel()


def run_solver(system, method: str, dx: float, dy: float, tol: float, omega=None,
               max_iter: int = iterative.DEFAULT_MAX_ITER, callback=None):
    """Solve with one method; returns (p, history)."""
    return iterative.solve(system, method, tol, max_iter, omega, spacing=(dx, dy), callback=callback)


def legend_text(key: str, h) -> str:
    """Legend entry with iterations and time (or 'not converged')."""
    label = style_of(key)["label"]
    if key == "sor":
        label = rf"SOR ($\omega$ = {h.omega:.3f})"
    if h.converged:
        return f"{label}: {h.iterations} it., {h.time_s:.2f} s"
    state = "diverged" if h.diverged else f"not converged (cap {h.iterations} it.)"
    return f"{label}: {state}"


# --------------------------------------------------------------------------
def part_history(args, system, dx, dy):
    """Figure 11 data: histories of all solvers on the baseline system."""
    hist = {}
    for m in SOLVERS:
        _, hist[m] = run_solver(system, m, dx, dy, args.tol)
        print(f"  {m:6s} {hist[m].iterations:7d} it.  {hist[m].time_s:.2f} s  converged={hist[m].converged}")
    cap = 2 * hist["jacobi"].iterations
    _, hist["jacobi_classic"] = run_solver(system, "jacobi", dx, dy, args.tol, omega=1.0, max_iter=cap)
    h = hist["jacobi_classic"]
    print(f"  jacobi(omega=1) {h.iterations} it. converged={h.converged} "
          f"final scaled residual {h.residual_scaled[-1]:.2e}")
    rows = []
    for key, h in hist.items():
        every = max(1, h.iterations // 4000)
        rows += [dict(r, solver=key, n=system.shape[1], time_s=h.time_s, converged=h.converged)
                 for r in h.to_rows(every)]
    write_csv(rows, args.outdir / "residual_history.csv")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.0))
    for key in ("jacobi", "jacobi_classic", "gs", "sor", "cg"):
        h = hist[key]
        k = np.arange(h.iterations + 1)
        ls = "--" if key == "jacobi_classic" else "-"
        plot_history(a1, h.residual_scaled, key, x=k, label=legend_text(key, h), ls=ls)
        plot_history(a2, h.residual_scaled[1:], key, x=k[1:], label=legend_text(key, h), ls=ls)
    for ax in (a1, a2):
        ax.axhline(args.tol, color=INK_2, ls="--", lw=1.0, label=f"tolerance {args.tol:g}")
        ax.set_ylabel(r"scaled residual L2$_n$(r) / S (-)")
    a1.set_xlim(0, max(hist["sor"].iterations, hist["cg"].iterations) * 1.6)
    a1.set_xlabel("iteration k")
    a1.set_title(f"(a) first iterations (linear axis), Taylor-Green {system.shape[1]}$^2$")
    a2.set_xscale("log")
    a2.set_xlabel("iteration k (log axis)")
    a2.set_title("(b) full history")
    a2.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    save_figure(fig, "fig11_residual_history")
    return hist


def part_norms(args, hist):
    """Figure 12: norms of one history + normalised vs unnormalised across grids."""
    h = hist["sor"]
    rows = h.to_rows(1)
    write_csv([dict(r, n=args.n) for r in rows], args.outdir / "residual_norms.csv")
    scaling = []
    for n in args.ns:
        system, *_ = baseline_system(n)
        r0 = system.residual(system.p0)             # initial residual = b (zero initial guess)
        interior = r0.reshape(n, n)[1:-1, 1:-1]     # rows of the 5-point Laplacian only
        scaling.append({"n": n, "N": interior.size, "l1": l1(interior), "l2": l2(interior),
                        "l1n": l1_normalised(interior), "l2n": l2_normalised(interior),
                        "l2n_with_boundary_rows": l2_normalised(r0)})
    write_csv(scaling, args.outdir / "norm_scaling.csv")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.0))
    k = np.arange(h.iterations + 1)
    for key, y in (("l1n", h.residual_l1n), ("l2n", h.residual_l2n), ("linf", h.residual_linf),
                   ("scaled", h.residual_scaled)):
        plot_history(a1, y, key, x=k)
    a1.set_xlabel("iteration k")
    a1.set_ylabel(r"residual norm (units of $\nabla^2 p$; scaled: -)")
    a1.set_title(rf"(a) SOR ($\omega$ = {h.omega:.3f}), Taylor-Green {args.n}$^2$")
    a1.legend(loc="upper right", fontsize=7)
    N = [s["N"] for s in scaling]
    for key in ("l1", "l2", "l1n", "l2n"):
        vals = [s[key] for s in scaling]
        plot_series(a2, N, vals, key, label=f"{style_of(key)['label']} (slope {convergence_slope(N, vals):.2f})")
    vals = [s["l2n_with_boundary_rows"] for s in scaling]
    plot_series(a2, N, vals, "l2n", ls=":", mfc="white",
                label=f"L2$_n$ incl. ghost-point boundary rows (slope {convergence_slope(N, vals):.2f})")
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlabel("number of unknowns N = n$^2$ (-)")
    a2.set_ylabel(r"norm of the initial residual $r_0 = b$ (interior rows)")
    a2.set_title("(b) same physical problem on finer grids")
    a2.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    save_figure(fig, "fig12_norm_comparison")


def part_change(args, system, dx, dy, p_direct):
    """Figure 13: Jacobi residual vs change norm vs true error."""
    errors = [l2_normalised(system.p0 - p_direct)]

    def track(p, k):
        errors.append(l2_normalised(p - p_direct))

    _, h = run_solver(system, "jacobi", dx, dy, args.tol, callback=track)
    rows = []
    every = max(1, h.iterations // 4000)
    for k in list(range(0, h.iterations + 1, every)) + ([h.iterations] if h.iterations % every else []):
        rows.append({"iteration": k, "residual_scaled": float(h.residual_scaled[k]),
                     "change_l2n": float(h.change_l2n[k - 1]) if k > 0 else float("nan"),
                     "error_l2n": errors[k]})
    write_csv(rows, args.outdir / "residual_vs_change.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    k = np.arange(h.iterations + 1)
    plot_history(ax, h.residual_scaled, "scaled", x=k)
    plot_history(ax, h.change_l2n, "change", x=k[1:])
    plot_history(ax, errors, "error", x=k)
    eps = args.change_eps
    below = np.nonzero(h.change_l2n < eps)[0]
    if below.size:
        k_eps = int(below[0]) + 1
        ax.axvline(k_eps, color=INK_2, ls=":", lw=1.0)
        ax.axhline(eps, color=INK_2, ls="--", lw=1.0, label=rf"naive $\varepsilon$ = {eps:g} on the change")
        ax.annotate(f"change < {eps:g} at k = {k_eps}\nbut scaled residual = {h.residual_scaled[k_eps]:.1e}"
                    f"\nand true error = {errors[k_eps]:.1e}", xy=(k_eps, eps), xytext=(0.42, 0.60),
                    textcoords="axes fraction", fontsize=7, color=INK_2,
                    arrowprops={"arrowstyle": "->", "color": INK_2},
                    bbox={"boxstyle": "round", "fc": "white", "ec": "#d9d8d4"})
    ax.set_xlabel("iteration k")
    ax.set_ylabel("norm (pressure in units of $\\rho U^2$)")
    ax.set_title(rf"Jacobi ($\omega$ = {h.omega:g}): residual vs change per iteration")
    ax.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    save_figure(fig, "fig13_residual_vs_change")


def part_grid(args):
    """Figure 14: iterations to reach tol vs N."""
    rows = []
    for n in args.ns:
        system, dx, dy, _ = baseline_system(n)
        for m in SOLVERS:
            if m in ("jacobi", "gs") and n > args.slow_max_n:
                continue
            _, h = run_solver(system, m, dx, dy, args.tol)
            rows.append({"n": n, "N": n * n, **h.summary()})
            print(f"  n={n:4d} {m:6s} {h.iterations:7d} it.  {h.time_s:7.2f} s")
    write_csv(rows, args.outdir / "iterations_vs_grid.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for m in SOLVERS:
        sub = [r for r in rows if r["method"] == m and r["converged"]]
        if not sub:
            continue
        N = [r["N"] for r in sub]
        it = [r["iterations"] for r in sub]
        a = convergence_slope(N, it)
        plot_series(ax, N, it, m, label=f"{style_of(m)['label']}: iterations $\\propto N^{{{a:.2f}}}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of unknowns N = n$^2$ (-)")
    ax.set_ylabel(f"iterations to reach tol = {args.tol:g}")
    ax.set_title("Cost of the iterative solvers vs grid size")
    ax.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    save_figure(fig, "fig14_iterations_vs_grid")


def part_sor(args, system, dx, dy):
    """Figure 15: SOR iterations vs omega and divergence for omega >= 2."""
    n = system.shape[1]
    w_neu = iterative.omega_opt(n, n, dx, dy, neumann=True)
    w_dir = iterative.omega_opt_dirichlet_formula(dx, (n - 1) * dx)
    omegas = sorted(set(args.omegas) | {round(w_neu, 4), round(w_dir, 4)})
    rows = []
    for w in omegas:
        _, h = iterative.sor(system, w, args.tol, args.omega_max_iter)
        rows.append({"omega": w, **{k: v for k, v in h.summary().items() if k != "omega"}})
    write_csv(rows, args.outdir / "sor_omega.csv")
    div_rows, div_hist = [], {}
    for w in args.diverge_omegas:
        _, h = iterative.sor(system, w, args.tol, args.diverge_cap)
        div_hist[w] = h
        div_rows += h.to_rows(1)
    write_csv(div_rows, args.outdir / "sor_divergence.csv")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.0))
    conv = [r for r in rows if r["converged"]]
    plot_series(a1, [r["omega"] for r in conv], [r["iterations"] for r in conv], "sor",
                label="SOR (red-black)")
    best = min(conv, key=lambda r: r["iterations"])
    a1.axvline(w_neu, color=style_of("sor")["color"], ls="--", lw=1.0,
               label=rf"$\omega_{{opt}}$ Neumann estimate = {w_neu:.3f}")
    a1.axvline(w_dir, color=INK_2, ls=":", lw=1.2,
               label=rf"$\omega_{{opt}}$ = 2/(1+sin($\pi h/L$)) = {w_dir:.3f}")
    a1.plot([best["omega"]], [best["iterations"]], marker="*", ms=12, color="#0b0b0b", ls="",
            label=rf"measured best $\omega$ = {best['omega']:g} ({best['iterations']} it.)")
    a1.set_yscale("log")
    a1.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    a1.yaxis.set_major_formatter(ScalarFormatter())
    a1.set_xlabel(r"relaxation factor $\omega$ (-)")
    a1.set_ylabel(f"iterations to reach tol = {args.tol:g}")
    a1.set_title(f"(a) SOR convergence vs $\\omega$, Taylor-Green {n}$^2$")
    a1.legend(loc="upper left", fontsize=7)
    for i, (w, h) in enumerate(div_hist.items()):
        state = "diverged" if h.diverged else "no convergence"
        plot_history(a2, h.residual_scaled, "sor", color=["#e34948", "#4a3aa7", "#eda100"][i % 3],
                     marker="", label=rf"$\omega$ = {w:g}: {state} ({h.iterations} it.)")
    a2.axhline(1.0, color=INK_2, ls="--", lw=1.0, label="initial residual")
    a2.set_xlabel("iteration k")
    a2.set_ylabel(r"scaled residual L2$_n$(r) / S (-)")
    a2.set_title(f"(b) $\\omega \\geq 2$: residual grows (run capped at {args.diverge_cap} it.)")
    a2.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    save_figure(fig, "fig15_sor_stability")
    print(f"  omega_opt Neumann {w_neu:.4f}, formula {w_dir:.4f}, measured best {best['omega']} "
          f"({best['iterations']} it.)")


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[0], ("tg",))
    ap.add_argument("--n", type=int, default=64, help="baseline grid (n x n)")
    ap.add_argument("--ns", type=int, nargs="+", default=[32, 64, 128, 256],
                    help="grids for the norm-scaling and iterations-vs-N studies")
    ap.add_argument("--slow-max-n", type=int, default=128,
                    help="largest grid for Jacobi / Gauss-Seidel (256 takes ~25 min)")
    ap.add_argument("--tol", type=float, default=iterative.DEFAULT_TOL)
    ap.add_argument("--change-eps", type=float, default=1e-6,
                    help="naive steady-state threshold on the change norm (Figure 13)")
    ap.add_argument("--omegas", type=float, nargs="+",
                    default=[0.5, 0.75, 1.0, 1.25, 1.5, 1.7, 1.8, 1.85, 1.9, 1.92, 1.94, 1.95,
                             1.96, 1.97, 1.98, 1.99])
    ap.add_argument("--omega-max-iter", type=int, default=200_000)
    ap.add_argument("--diverge-omegas", type=float, nargs="+", default=[2.0, 2.1])
    ap.add_argument("--diverge-cap", type=int, default=300)
    args = ap.parse_args(argv)
    if args.quick:
        args.n, args.ns, args.slow_max_n = 32, [16, 32, 64], 64
        args.omegas = [0.75, 1.0, 1.5, 1.8, 1.9, 1.95, 1.99]
    set_style()
    system, dx, dy, p_direct = baseline_system(args.n)
    hist = part_history(args, system, dx, dy)
    part_norms(args, hist)
    part_change(args, system, dx, dy, p_direct)
    part_sor(args, system, dx, dy)
    part_grid(args)


if __name__ == "__main__":
    main()
