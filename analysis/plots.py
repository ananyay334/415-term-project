# piv_pressure/analysis/plots.py
"""Shared Matplotlib style and plotting helpers (every figure uses these).

Design rules (applied everywhere):
* One fixed colour per *entity* (method, scheme, BC, ...), never re-assigned
  by rank; every series also has its own marker so identity never relies on
  colour alone (the colour-blind-safe palette was checked for CVD separation).
* Pressure fields: diverging blue-grey-red map; error magnitudes: single-hue
  sequential map.  Fields that are compared share one colour scale.
* Figures are saved as PNG (300 dpi) and PDF in ``results/figures``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

if "ipykernel" not in sys.modules:  # scripts: headless, figures only go to files
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
#: Output root; the environment variable PIV_RESULTS_DIR redirects all
#: CSVs and figures (used e.g. for smoke tests that must not touch results/).
RESULTS_DIR = Path(os.environ.get("PIV_RESULTS_DIR", PROJECT_ROOT / "results"))
FIG_DIR = RESULTS_DIR / "figures"

#: Categorical palette in fixed slot order (CVD-validated, light surface).
PALETTE: tuple[str, ...] = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                            "#4a3aa7", "#e34948")
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d9d8d4"

#: Entity -> (colour, marker, label).  Colour follows the entity in every figure.
STYLE: dict[str, tuple[str, str, str]] = {
    "poisson": (PALETTE[0], "o", "Poisson"),
    "direct": (PALETTE[1], "s", "Direct (path avg.)"),
    "omni": (PALETTE[2], "^", "Omni-directional"),
    "direct_path1": (PALETTE[3], "v", "Direct path 1 (x then y)"),
    "direct_path2": (PALETTE[4], "D", "Direct path 2 (y then x)"),
    "poisson_velocity": (PALETTE[5], "P", "Poisson (velocity source)"),
    "eulerian": (PALETTE[0], "o", "Eulerian"),
    "lagrangian": (PALETTE[1], "s", "Lagrangian"),
    "neumann": (PALETTE[0], "o", "All-Neumann (pin mean)"),
    "neumann_node": (PALETTE[3], "v", "All-Neumann (pin node)"),
    "bernoulli_1": (PALETTE[1], "s", "Bernoulli on 1 edge"),
    "exact_1": (PALETTE[2], "^", "Exact p on 1 edge"),
    "bernoulli_2": (PALETTE[4], "D", "Bernoulli on 2 opposite edges"),
    "exact_2": (PALETTE[5], "P", "Exact p on 2 opposite edges"),
    "none": (PALETTE[0], "o", "No denoising"),
    "gaussian": (PALETTE[1], "s", "Gaussian"),
    "savgol": (PALETTE[2], "^", "Savitzky-Golay"),
    "clean": (PALETTE[0], "o", "Clean"),
    "noisy": (PALETTE[1], "s", "Noisy"),
    "reference": (INK, "", "Reference"),
    # test cases (used when the series of a panel are cases)
    "case:tg": (PALETTE[0], "o", "Taylor-Green"),
    "case:tg_core": (PALETTE[0], "o", "Taylor-Green (core window)"),
    "case:tgc": (PALETTE[4], "D", "Taylor-Green, convected (A=0.4)"),
    "case:kov20": (PALETTE[1], "s", "Kovasznay Re=20"),
    "case:kov40": (PALETTE[3], "v", "Kovasznay Re=40"),
    "case:cfd": (PALETTE[2], "^", "CFD cylinder wake"),
    # iterative Poisson solvers
    "jacobi": (PALETTE[0], "o", r"Jacobi ($\omega$ = 0.95)"),
    "gs": (PALETTE[1], "s", "Gauss-Seidel"),
    "sor": (PALETTE[2], "^", "SOR"),
    "cg": (PALETTE[3], "D", "Conjugate gradient"),
    "jacobi_classic": (PALETTE[6], "v", r"Jacobi, classic ($\omega$ = 1)"),
    # residual / error norms
    "l1n": (PALETTE[0], "o", r"L1$_n$ (mean $|r|$)"),
    "l2n": (PALETTE[1], "s", r"L2$_n$ (rms $r$)"),
    "linf": (PALETTE[2], "^", r"L$_\infty$ (max $|r|$)"),
    "scaled": (PALETTE[3], "D", r"scaled L2$_n$ = L2$_n$ / S"),
    "l1": (PALETTE[4], "v", "L1 (unnormalised)"),
    "l2": (PALETTE[5], "P", "L2 (unnormalised)"),
    "change": (PALETTE[1], "s", r"change $\|p^{k+1}-p^{k}\|$ (L2$_n$)"),
    "error": (PALETTE[2], "^", r"true error $\|p^{k}-p^{*}\|$ (L2$_n$)"),
    # scheme order
    "order1": (PALETTE[1], "s", "first order (order=1)"),
    "order2": (PALETTE[0], "o", "second order (order=2)"),
    # CFD residuals and meshes
    "Ux": (PALETTE[0], "", "Ux"),
    "Uy": (PALETTE[1], "", "Uy"),
    "p_res": (PALETTE[2], "", "p"),
    "coarse": (PALETTE[3], "v", "coarse"),
    "medium": (PALETTE[1], "s", "medium"),
    "fine": (PALETTE[0], "o", "fine"),
}

CMAP_DIVERGING = LinearSegmentedColormap.from_list(
    "pressure_div", ["#104281", "#3987e5", "#f0efec", "#e34948", "#8a1c1c"])
CMAP_SEQUENTIAL = LinearSegmentedColormap.from_list(
    "error_seq", ["#f4f8fd", "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"])


def set_style() -> None:
    """Apply the shared rcParams (call once at the start of every script)."""
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
        "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "lines.linewidth": 1.5, "lines.markersize": 5,
        "legend.frameon": False, "image.cmap": "viridis",
        "axes.prop_cycle": matplotlib.cycler(color=list(PALETTE)),
    })


def style_of(key: str) -> dict:
    """Keyword arguments (color, marker, label) for a plotted entity."""
    color, marker, label = STYLE.get(key, (INK_2, "x", key))
    return {"color": color, "marker": marker, "label": label}


def save_figure(fig: plt.Figure, name: str, outdir: Path | None = None) -> list[Path]:
    """Save ``fig`` as ``<name>.png`` (300 dpi) and ``<name>.pdf``; close it."""
    outdir = FIG_DIR if outdir is None else Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = [outdir / f"{name}.png", outdir / f"{name}.pdf"]
    for p in paths:
        fig.savefig(p)
    plt.close(fig)
    return paths


def add_slope_triangle(ax: plt.Axes, h: Sequence[float], err: Sequence[float], order: float = 2.0,
                       label: str | None = None) -> None:
    """Draw a reference line err ~ h^order through the finest clean point."""
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    k = int(np.argmin(h))
    href = np.array([h.min(), h.max()])
    ax.plot(href, err[k] * (href / h[k]) ** order * 0.6, color=INK_2, ls="--", lw=1.0,
            marker="", label=label or f"slope {order:g}")


def draw_slope_triangle(ax: plt.Axes, x0: float, y0: float, order: float, factor: float = 2.0,
                        color: str = INK_2) -> None:
    """Right-angled slope triangle on log-log axes with the order written on it.

    Legs: (x0, y0) -> (x0*factor, y0) -> (x0*factor, y0*factor**order).  In
    log-log coordinates the hypotenuse has slope ``order``.  Works on reversed
    axes too (the label is offset in screen points).
    """
    x1, y1 = x0 * factor, y0 * factor**order
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color=color, lw=1.0, marker="", zorder=3)
    ax.annotate(f"{order:g}", xy=(x1, np.sqrt(y0 * y1)), xytext=(4, 0), textcoords="offset points",
                color=color, fontsize=8, va="center", ha="left")
    ax.annotate("1", xy=(np.sqrt(x0 * x1), y0), xytext=(0, -9), textcoords="offset points",
                color=color, fontsize=8, ha="center")


def threshold_guides(ax: plt.Axes, x_cross: float, level: float, color: str,
                     label: str | None = None, extrapolated: bool = False) -> None:
    """Dotted guide from the curve's crossing of ``level`` down to the x axis."""
    y_bottom = ax.get_ylim()[0]
    ls = (0, (1, 3)) if extrapolated else ":"
    ax.plot([x_cross, x_cross], [y_bottom, level], color=color, ls=ls, lw=1.2, marker="", label=label)
    ax.plot([x_cross], [level], marker="o", mfc="white", mec=color, ms=7, ls="")
    ax.set_ylim(bottom=y_bottom)


def plot_history(ax: plt.Axes, y: Sequence[float], key: str, x: Sequence[float] | None = None,
                 label: str | None = None, n_markers: int = 8, **kw) -> None:
    """Residual/iteration history on a log y axis (markers thinned out)."""
    st = style_of(key)
    if label is not None:
        st["label"] = label
    st.update(kw)
    y = np.asarray(y, dtype=float)
    x = np.arange(y.size) if x is None else np.asarray(x, dtype=float)
    st.setdefault("markevery", max(1, y.size // n_markers))
    st.setdefault("lw", 1.3)
    ax.plot(x, y, **st)
    ax.set_yscale("log")


def plot_series(ax: plt.Axes, x: Sequence[float], y: Sequence[float], key: str,
                yerr: Sequence[float] | None = None, label: str | None = None, **kw) -> None:
    """Line + markers (+ optional +/- std band) in the entity's fixed style."""
    st = style_of(key)
    if label is not None:
        st["label"] = label
    st.update(kw)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ax.plot(x, y, **st)
    if yerr is not None:
        yerr = np.asarray(yerr, dtype=float)
        ax.fill_between(x, np.maximum(y - yerr, 1e-300), y + yerr, color=st["color"], alpha=0.18, lw=0)


def field_row(
    axes: Iterable[plt.Axes],
    x: np.ndarray,
    y: np.ndarray,
    fields: Sequence[np.ndarray],
    titles: Sequence[str],
    clim: tuple[float, float] | None = None,
    cmap: LinearSegmentedColormap = CMAP_DIVERGING,
    cbar_label: str = "",
    xlabel: str = "x (non-dimensional)",
    ylabel: str = "y (non-dimensional)",
) -> None:
    """Plot several fields with ONE shared colour scale and one colour bar."""
    axes = list(axes)
    if clim is None:
        vmax = max(np.max(np.abs(f)) for f in fields)
        clim = (-vmax, vmax)
    im = None
    for ax, f, title in zip(axes, fields, titles):
        im = ax.pcolormesh(x, y, f, cmap=cmap, vmin=clim[0], vmax=clim[1], shading="nearest",
                           rasterized=True)
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlabel(xlabel)
        ax.grid(False)
    axes[0].set_ylabel(ylabel)
    cb = axes[0].figure.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label(cbar_label)


def error_panel(ax: plt.Axes, x: np.ndarray, y: np.ndarray, err: np.ndarray, title: str,
                label: str = "|p - p_ref| (non-dimensional)",
                xlabel: str = "x (non-dimensional)") -> None:
    """Error-magnitude map with its own sequential colour bar."""
    im = ax.pcolormesh(x, y, err, cmap=CMAP_SEQUENTIAL, vmin=0.0, shading="nearest", rasterized=True)
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.grid(False)
    cb = ax.figure.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(label)
