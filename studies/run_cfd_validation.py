# piv_pressure/studies/run_cfd_validation.py
"""CFD validation: reconstruct the cylinder-wake pressure from the exported
CFD velocity window and compare with the CFD pressure (Figures 5, 6, 10).

* Figure 5 - p_CFD, p_reconstructed (Poisson, Eulerian) and the error map at
  one instant, first two on one shared colour scale.
* Figure 6 - pressure along the horizontal line through a vortex core
  (the pressure minimum inside x in [2, 4.8], |y| < 1.5): CFD vs Poisson
  (Eulerian and Lagrangian), direct and omni integration.
* Error history eps_L2(t) over one shedding period (Eulerian, Lagrangian).
* Figure 10 (``--animate``) - GIF of CFD vs reconstructed pressure over one
  shedding cycle.

The data are the "clean" CFD samples at 0.05 D (no synthetic PIV errors);
the parametric studies add those.  The script exits quietly if the CFD
export is missing.

Usage:  python studies/run_cfd_validation.py [--time-index N] [--animate]
Output: results/cfd_validation.csv, fig05_cfd_wake, fig06_vortex_core_line,
        fig05b_cfd_error_history, fig10_cfd_animation.gif.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

from analysis.metrics import core_pressure, demean, error_map, line_profile, standard_metrics  # noqa: E402
from analysis.plots import (CMAP_DIVERGING, FIG_DIR, RESULTS_DIR, error_panel, field_row,  # noqa: E402
                            plot_series, save_figure, set_style)
from core.pipeline import reconstruct  # noqa: E402
from fields.cfd_loader import DEFAULT_CFD_NPZ, cfd_available  # noqa: E402
from studies.common import CFD_CORE_REGION, base_parser, get_cfd, write_csv  # noqa: E402

XLAB, YLAB = "x / D (-)", "y / D (-)"


def shedding_period(fld) -> float:
    """Period from cfd_summary.json (grid study) or from the FFT of v at a probe."""
    summary = RESULTS_DIR / "cfd" / "cfd_summary.json"
    if summary.is_file():
        return float(json.loads(summary.read_text())["period"])
    j, i = fld.ny // 2, fld.nx // 2
    v = fld.v[:, j, i] - fld.v[:, j, i].mean()
    spec = np.abs(np.fft.rfft(v * np.hanning(v.size), 8 * v.size))
    freq = np.fft.rfftfreq(8 * v.size, fld.dt)
    return 1.0 / freq[1 + np.argmax(spec[1:])]


def figure_wake(fld, n: int) -> dict:
    """Figures 5 and 6 at snapshot n; returns the metrics of each method."""
    recs = {
        "poisson": reconstruct(fld, method="poisson", time_index=n),
        "lagrangian": reconstruct(fld, method="poisson", material_derivative="lagrangian", time_index=n),
        "direct": reconstruct(fld, method="direct", time_index=n),
        "omni": reconstruct(fld, method="omni", time_index=n),
    }
    p_ref = demean(fld.p_ref[n])
    p_rec = demean(recs["poisson"].p)
    vmax = max(np.abs(p_ref).max(), np.abs(p_rec).max())
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6))
    field_row(axes[:2], fld.x, fld.y, [p_ref, p_rec],
              [f"CFD pressure, t = {fld.t[n]:.2f}", "reconstructed (Poisson, Eulerian)"],
              clim=(-vmax, vmax), cbar_label=r"$(p - \bar p)/\rho U^2$", xlabel=XLAB, ylabel=YLAB)
    error_panel(axes[2], fld.x, fld.y, error_map(p_rec, p_ref), "error map",
                label=r"$|p - p_{CFD}|/\rho U^2$", xlabel=XLAB)
    save_figure(fig, "fig05_cfd_wake")

    core = core_pressure(recs["poisson"].p, fld.p_ref[n], fld.x, fld.y, region=CFD_CORE_REGION)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    x, pr, y_used = line_profile(p_ref, fld.x, fld.y, core["y_core"])
    ax.plot(x, pr, color="#0b0b0b", lw=2.0, label="CFD")
    for key, label in (("poisson", "Poisson, Eulerian"), ("lagrangian", "Poisson, Lagrangian"),
                       ("direct", "Direct (path avg.)"), ("omni", "Omni-directional")):
        _, pk, _ = line_profile(demean(recs[key].p), fld.x, fld.y, core["y_core"])
        plot_series(ax, x, pk, key, label=label, markevery=6, lw=1.2)
    ax.axvline(core["x_core"], color="#52514e", ls=":", lw=1.0, label="vortex core")
    ax.set_xlabel(XLAB)
    ax.set_ylabel(r"$(p - \bar p)/\rho U^2$")
    ax.set_title(f"Pressure along y/D = {y_used:.2f} through a vortex core (t = {fld.t[n]:.2f})")
    ax.legend(loc="best", fontsize=7)
    save_figure(fig, "fig06_vortex_core_line")

    out = {}
    for key, rec in recs.items():
        m = standard_metrics(rec.p, fld.p_ref[n], fld.rho, fld.U_ref)
        m.update(core_pressure(rec.p, fld.p_ref[n], fld.x, fld.y, region=CFD_CORE_REGION))
        m["time_s"] = rec.diagnostics["time_total_s"]
        out[key] = m
    return out


def error_history(fld, indices: list[int]) -> list[dict]:
    """eps_L2(t) for the Eulerian and Lagrangian Poisson reconstructions."""
    rows = []
    for n in indices:
        for scheme in ("eulerian", "lagrangian"):
            rec = reconstruct(fld, method="poisson", material_derivative=scheme, time_index=n)
            row = standard_metrics(rec.p, fld.p_ref[n], fld.rho, fld.U_ref)
            row.update(core_pressure(rec.p, fld.p_ref[n], fld.x, fld.y, region=CFD_CORE_REGION))
            row.update(kind="history", t=float(fld.t[n]), time_index=n, scheme=scheme,
                       compatibility_relative=rec.diagnostics["compatibility_relative"])
            rows.append(row)
    return rows


def animate(fld, indices: list[int], path: Path) -> None:
    """GIF of CFD vs reconstructed pressure over the given snapshots."""
    frames = [(demean(fld.p_ref[n]), demean(reconstruct(fld, time_index=n).p), fld.t[n]) for n in indices]
    vmax = max(max(np.abs(a).max(), np.abs(b).max()) for a, b, _ in frames)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6))
    ims = []
    for ax, title in zip(axes, ("CFD", "reconstructed (Poisson)")):
        ims.append(ax.pcolormesh(fld.x, fld.y, frames[0][0], cmap=CMAP_DIVERGING, vmin=-vmax, vmax=vmax,
                                 shading="nearest"))
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlabel(XLAB)
        ax.grid(False)
    axes[0].set_ylabel(YLAB)
    fig.colorbar(ims[1], ax=axes, shrink=0.85, label=r"$(p - \bar p)/\rho U^2$")
    stamp = fig.suptitle("")

    def update(k: int):
        a, b, t = frames[k]
        ims[0].set_array(a.ravel())
        ims[1].set_array(b.ravel())
        stamp.set_text(f"t U/D = {t:.2f}")
        return ims

    path.parent.mkdir(parents=True, exist_ok=True)
    FuncAnimation(fig, update, frames=len(frames), blit=False).save(path, writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"[gif] {path}")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = base_parser(__doc__.splitlines()[0], ("cfd",))
    ap.add_argument("--time-index", type=int, default=None, help="snapshot for Figs 5/6 (default: middle)")
    ap.add_argument("--animate", action="store_true", help="write the GIF (Figure 10)")
    args = ap.parse_args(argv)
    if not cfd_available():
        print(f"[skip] CFD data not found at {DEFAULT_CFD_NPZ} - run the CFD case first.")
        return []
    set_style()
    fld = get_cfd()
    n = fld.nt // 2 if args.time_index is None else args.time_index
    rows = []
    for key, m in figure_wake(fld, n).items():
        rows.append(dict(m, kind="snapshot", method=key, t=float(fld.t[n]), time_index=n))
    period = shedding_period(fld)
    n_period = int(round(period / fld.dt))
    stride = 4 if args.quick else 1
    idx = [k for k in range(n, min(n + n_period, fld.nt - 1), stride) if k >= 1]
    hist = error_history(fld, idx)
    rows += hist
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    for scheme in ("eulerian", "lagrangian"):
        tt = [r["t"] for r in hist if r["scheme"] == scheme]
        ee = [r["eps_l2"] for r in hist if r["scheme"] == scheme]
        plot_series(ax, tt, ee, scheme, markevery=max(1, len(tt) // 12))
    ax.set_xlabel(r"time $t U/D$ (-)")
    ax.set_ylabel(r"$\epsilon_{L2}$ (-)")
    ax.set_title(f"Reconstruction error over one shedding period (T = {period:.2f} D/U)")
    ax.legend(loc="best")
    save_figure(fig, "fig05b_cfd_error_history")
    write_csv(rows, args.outdir / "cfd_validation.csv")
    if args.animate:
        animate(fld, idx if not args.quick else idx[::2], FIG_DIR / "fig10_cfd_animation.gif")
    for key in ("poisson", "lagrangian", "direct", "omni"):
        r = next(r for r in rows if r.get("method") == key)
        print(f"  {key:10s} eps_L2 {r['eps_l2']:.3e} eps_max {r['eps_max']:.3e} core ratio {r['core_ratio']:.3f}")
    return rows


if __name__ == "__main__":
    main()
