# piv_pressure/studies/run_cfd_convergence.py
"""CFD solver convergence and grid independence (Figures 17 and 18).

Figure 17  (top) initial residuals of Ux, Uy and p vs simulation time for the
           coarse, medium and fine meshes (from the OpenFOAM logs);
           (bottom) Cl(t) on the finest mesh with the cycle-to-cycle
           periodicity check annotated, the relative cycle-to-cycle changes
           of amplitude and period for every mesh, and - if the optional
           steady Re = 20 case was run - its residuals and ||u^{n+1} - u^n||
           with the tolerance eps and the time at which the run stopped.
Figure 18  St and mean Cd vs number of cells with the Richardson-extrapolated
           value (dashed) and GCI error bars on the medium and fine meshes.

Needs the OpenFOAM runs (``sh cfd/run_cases.sh``); every part whose data
are missing is skipped with a message.  Run ``python cfd/grid_independence.py``
first for Figure 18 (run_all.py does this).

Usage:  python studies/run_cfd_convergence.py [--tol 1e-3] [--eps 1e-6]
Output: results/cfd/residuals_<case>.csv, periodicity_<mesh>.csv,
        change_norm_window.csv, steady_re20_change.csv, time_loop_summary.json;
        figures fig17_cfd_residuals, fig18_cfd_grid_convergence.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from analysis.plots import INK_2, PALETTE, RESULTS_DIR, save_figure, set_style, style_of  # noqa: E402
from cfd import change_norm, grid_independence, parse_log  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MESHES = ("coarse", "medium", "fine")


def _read(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _col(rows: list[dict], key: str) -> np.ndarray:
    return np.array([float(r[key]) if r[key] not in ("", "nan") else np.nan for r in rows])


def figure_residuals(residuals: dict, cases: Path, cfd_dir: Path, summary: dict, eps: float) -> bool:
    """Figure 17; returns False if there is nothing to plot."""
    meshes = [m for m in MESHES if m in residuals]
    if not meshes:
        print("[skip] no OpenFOAM logs for the Re = 100 meshes - Figure 17 not drawn.")
        return False
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.2))
    for ax, mesh in zip(axes[0], MESHES):
        if mesh not in residuals:
            ax.axis("off")
            ax.text(0.5, 0.5, f"{mesh}: no log", ha="center", transform=ax.transAxes)
            continue
        res = residuals[mesh]
        for f, key in (("Ux", "Ux"), ("Uy", "Uy"), ("p", "p_res")):
            ax.plot(res["t"], res[f], color=style_of(key)["color"], lw=0.8, label=f)
        ax.set_yscale("log")
        ax.set_xlabel(r"time $t U/D$ (-)")
        ax.set_ylabel("initial residual (-)")
        ax.set_title(f"({'abc'[MESHES.index(mesh)]}) {mesh} mesh, Re = 100")
        ax.legend(loc="upper right", fontsize=7)
    # Cl(t) with periodicity on the finest available mesh
    ax = axes[1, 0]
    last = meshes[-1]
    files = grid_independence.coefficient_files(cases / last)
    if files:
        t, _, cl = grid_independence.read_coefficients(files)
        per = change_norm.periodicity_check(t, cl, summary["tol"])
        ax.plot(t, cl, color=style_of(last)["color"], lw=0.8, label=f"Cl, {last} mesh")
        if per["satisfied"]:
            ax.axvline(per["t_periodic"], color=INK_2, ls="--", lw=1.0,
                       label=f"periodic from t = {per['t_periodic']:.1f}")
        ax.set_title(f"(d) lift coefficient; periodic when |dA/A|, |dT/T| < {summary['tol']:g}\n"
                     f"for {per['n_consecutive']} consecutive cycles", fontsize=8.5)
        ax.legend(loc="lower right", fontsize=7, frameon=True, framealpha=0.9)
    ax.set_xlabel(r"time $t U/D$ (-)")
    ax.set_ylabel(r"$C_L$ (-)")
    # cycle-to-cycle changes for every mesh
    ax = axes[1, 1]
    for mesh in MESHES:
        path = cfd_dir / f"periodicity_{mesh}.csv"
        if not path.is_file():
            continue
        rows = _read(path)
        st = style_of(mesh)
        ax.plot(_col(rows, "t_start"), _col(rows, "rel_change_amplitude"), color=st["color"],
                marker=st["marker"], ms=3, lw=0.9, label=f"{mesh}: |dA/A|")
        ax.plot(_col(rows, "t_start"), _col(rows, "rel_change_period"), color=st["color"],
                marker=st["marker"], ms=3, lw=0.9, ls="--", mfc="white", label=f"{mesh}: |dT/T|")
    ax.axhline(summary["tol"], color=INK_2, ls=":", lw=1.2, label=f"tolerance {summary['tol']:g}")
    ax.set_yscale("log")
    ax.set_xlabel(r"cycle start time $t U/D$ (-)")
    ax.set_ylabel("cycle-to-cycle relative change (-)")
    ax.set_title("(e) periodicity criterion", fontsize=9)
    ax.legend(loc="upper right", fontsize=6, ncol=2)
    # optional steady case
    ax = axes[1, 2]
    steady_csv = cfd_dir / "steady_re20_change.csv"
    if steady_csv.is_file() and "steady" in summary:
        rows = _read(steady_csv)
        ax.plot(_col(rows, "t"), _col(rows, "change_l2n"), color=PALETTE[1], lw=1.2, marker="s", ms=3,
                markevery=max(1, len(rows) // 12), label=r"$\|u^{n+1}-u^n\|$ / U (L2$_n$, window)")
        if "steady_re20" in residuals:
            r = residuals["steady_re20"]
            for f, key in (("Ux", "Ux"), ("Uy", "Uy"), ("p", "p_res")):
                ax.plot(r["t"], r[f], color=style_of(key)["color"], lw=0.7, alpha=0.8,
                        label=f"initial residual {f}")
        ax.axhline(eps, color=INK_2, ls="--", lw=1.0, label=rf"$\varepsilon$ = {eps:g}")
        s = summary["steady"]
        ax.axvline(s["t_stop"], color="#0b0b0b", ls=":", lw=1.2, label=f"run stopped, t = {s['t_stop']:.1f}")
        if np.isfinite(s["t_below_eps"]):
            ax.axvline(s["t_below_eps"], color=PALETTE[1], ls=":", lw=1.2,
                       label=rf"change < $\varepsilon$ at t = {s['t_below_eps']:.1f}")
        ax.set_yscale("log")
        ax.set_xlabel(r"time $t U/D$ (-)")
        ax.set_ylabel("residual / change (-)")
        ax.set_title("(f) steady Re = 20 case: steady-state criterion", fontsize=9)
        ax.legend(loc="upper right", fontsize=6)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "optional steady Re = 20 case not run\n(sh cfd/run_cases.sh N --steady)",
                ha="center", va="center", transform=ax.transAxes, color=INK_2)
    fig.tight_layout()
    save_figure(fig, "fig17_cfd_residuals")
    return True


def figure_grid(cfd_dir: Path) -> bool:
    """Figure 18 from grid_independence.csv and richardson.csv."""
    gi, rc = cfd_dir / "grid_independence.csv", cfd_dir / "richardson.csv"
    if not (gi.is_file() and rc.is_file()):
        print("[skip] grid_independence.csv / richardson.csv missing - Figure 18 not drawn.")
        return False
    rows = {r["level"]: r for r in _read(gi)}
    rich = {r["quantity"]: r for r in _read(rc)}
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))
    for ax, q, label in ((axes[0], "St", "Strouhal number St (-)"), (axes[1], "Cd_mean", r"mean $C_D$ (-)")):
        levels = [m for m in MESHES if m in rows]
        N = np.array([float(rows[m]["cells"]) for m in levels])
        val = np.array([float(rows[m][q]) for m in levels])
        R = rich.get(q, {})
        gci = {"medium": float(R.get("GCI_medium", "nan")), "fine": float(R.get("GCI_fine", "nan"))}
        ax.plot(N, val, color=PALETTE[0], lw=1.2, marker="", zorder=1)
        for n_, v_, m in zip(N, val, levels):
            st = style_of(m)
            err = gci.get(m, np.nan) * v_
            ax.errorbar([n_], [v_], yerr=[err] if np.isfinite(err) else None, marker=st["marker"],
                        color=st["color"], ls="", ms=7, capsize=4,
                        label=f"{m}" + (f" (GCI {100 * gci[m]:.2f} %)" if np.isfinite(gci.get(m, np.nan)) else ""))
        ext = float(R.get("extrapolated", "nan"))
        if np.isfinite(ext):
            ax.axhline(ext, color=INK_2, ls="--", lw=1.0,
                       label=f"Richardson extrapolated {ext:.4f} (p = {float(R['order']):.2f})")
        ax.set_xscale("log")
        ax.set_xlabel("number of cells (-)")
        ax.set_ylabel(label)
        ax.legend(loc="best", fontsize=7)
    axes[0].set_title("Strouhal number vs mesh size")
    axes[1].set_title("Mean drag coefficient vs mesh size")
    fig.tight_layout()
    save_figure(fig, "fig18_cfd_grid_convergence")
    return True


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", type=Path, default=ROOT / "cfd" / "cases")
    ap.add_argument("--outdir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--tol", type=float, default=1e-3, help="periodicity tolerance (cycle to cycle)")
    ap.add_argument("--eps", type=float, default=1e-6, help="steady-state tolerance on the change norm")
    ap.add_argument("--quick", action="store_true", help="accepted for run_all.py; no effect")
    ap.add_argument("--seeds", type=int, default=0, help="accepted for run_all.py; no effect")
    args = ap.parse_args(argv)
    set_style()
    cfd_dir = args.outdir / "cfd"
    residuals = parse_log.main(["--cases", str(args.cases), "--outdir", str(cfd_dir)])
    summary = change_norm.main(["--cases", str(args.cases), "--outdir", str(cfd_dir),
                                "--tol", str(args.tol), "--eps", str(args.eps)])
    figure_residuals(residuals, args.cases, cfd_dir, summary, args.eps)
    figure_grid(cfd_dir)


if __name__ == "__main__":
    main()
