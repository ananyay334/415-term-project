# piv_pressure/studies/common.py
"""Helpers shared by the study scripts: test-case factory, evaluation, CSV I/O.

Case keys used throughout the studies
-------------------------------------
``tg``      stationary Taylor-Green vortex on the spec window [0.3, 2.5]^2
``tg_core`` Taylor-Green on [pi-1.1, pi+1.1]^2 (contains the vortex core at
            (pi, pi); used for vortex-core metrics)
``tgc``     Taylor-Green of amplitude A = 0.4 convected with c = (1, 0)
            (exact; du/dt = O(1); like a wake vortex whose own velocity is
            smaller than the speed at which it is carried downstream -
            an illustrative choice, not a fitted wake property)
``kov20``, ``kov40``  Kovasznay flow at Re = 20 / 40
``cfd``     OpenFOAM cylinder wake (needs results/cfd/cylinder_re100_window.npz)

Baseline synthetic-PIV parameters: n = 64 points per side, window W = 1 h,
noise sigma = 1 % of U_ref, dt = 0.05, truth sampled 4x finer than the PIV
grid.  Noise seeds are 0..n_seeds-1 in every study ("common random
numbers": the same noise pattern is reused across parameter values, which
makes comparisons between them less noisy).
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from analysis.metrics import core_pressure, standard_metrics
from analysis.plots import RESULTS_DIR
from core.pipeline import reconstruct
from fields.analytic import (KOV_X, KOV_Y, TG_CORE_WINDOW, TG_WINDOW, kovasznay_on_grid,
                             taylor_green_on_grid)
from fields.cfd_loader import DEFAULT_CFD_NPZ, cfd_available, load_cfd
from fields.field import FlowField
from piv.degrade import degrade, degrade_analytic

BASE_N = 64
BASE_WINDOW = 1.0
BASE_SIGMA = 0.01
BASE_DT = 0.05
REFINE = 4
T_CENTER = 1.0
#: Vortex amplitude of the convected Taylor-Green case (convection speed = 1).
TGC_AMPLITUDE = 0.4
#: Region searched for the CFD vortex core (skips the base region near x = 1).
CFD_CORE_REGION = (2.0, 4.8, -1.5, 1.5)

CASE_LABELS = {
    "tg": "Taylor-Green", "tg_core": "Taylor-Green (core window)",
    "tgc": "Taylor-Green, convected (A=0.4)", "kov20": "Kovasznay Re=20",
    "kov40": "Kovasznay Re=40", "cfd": "CFD cylinder wake",
}
ANALYTIC_CASES = ("tg", "tg_core", "tgc", "kov20", "kov40")

_CFD_CACHE: dict[str, FlowField] = {}


# --------------------------------------------------------------------------
# Case factory
# --------------------------------------------------------------------------
def analytic_truth(case: str, dt: float = BASE_DT):
    """Return (truth(x, y) -> FlowField, x_range, y_range) for an analytic case."""
    times = T_CENTER + dt * np.array([-1.0, 0.0, 1.0])
    if case in ("tg", "tg_core", "tgc"):
        rng = TG_CORE_WINDOW if case == "tg_core" else TG_WINDOW
        advect = (1.0, 0.0) if case == "tgc" else (0.0, 0.0)
        amp = TGC_AMPLITUDE if case == "tgc" else 1.0
        name = {"tg": "taylor_green", "tg_core": "taylor_green_core", "tgc": "taylor_green_convected"}[case]
        return (lambda x, y: taylor_green_on_grid(x, y, times, advect=advect, name=name,
                                                  amplitude=amp)), rng, rng
    if case in ("kov20", "kov40"):
        Re = 20.0 if case == "kov20" else 40.0
        return (lambda x, y: kovasznay_on_grid(x, y, Re)), KOV_X, KOV_Y
    raise ValueError(f"unknown analytic case {case!r}")


def get_cfd(path: Path = DEFAULT_CFD_NPZ) -> FlowField:
    """Load (once) and return the CFD window."""
    key = str(path)
    if key not in _CFD_CACHE:
        _CFD_CACHE[key] = load_cfd(path)
    return _CFD_CACHE[key]


def cfd_piv(
    stride: int = 1,
    window_h: float = 0.0,
    sigma: float = 0.0,
    seed: int = 0,
    dt_multiple: int = 1,
    time_index: int | None = None,
    crop: int = 0,
) -> FlowField:
    """Synthetic PIV from the CFD window: three snapshots n-k, n, n+k.

    ``stride`` sets the PIV spacing (h = stride * 0.05 D); the box width is
    ``window_h * stride`` export points.  Edge nodes use nearest-value
    padding in the filter unless ``crop`` removes them.
    """
    fld = get_cfd()
    k = dt_multiple
    n = fld.nt // 2 if time_index is None else time_index
    if n - k < 0 or n + k >= fld.nt:
        raise ValueError("not enough CFD snapshots for this time step")
    sub = fld.subset_time([n - k, n, n + k])
    return degrade(sub, window=int(round(window_h * stride)), stride=stride,
                   sigma=sigma, seed=seed, crop=crop)


def make_case(
    case: str,
    n: int = BASE_N,
    window_h: float = BASE_WINDOW,
    sigma: float = BASE_SIGMA,
    seed: int = 0,
    dt: float = BASE_DT,
    stride: int = 1,
    refine: int = REFINE,
) -> FlowField:
    """Synthetic PIV field for any case key (see module docstring).

    For analytic cases ``n`` sets the PIV grid; for ``cfd`` ``stride`` sets
    it and ``dt`` is rounded to a multiple of the export step.
    """
    if case == "cfd":
        fld = get_cfd()
        k = max(1, int(round(dt / fld.dt)))
        return cfd_piv(stride, window_h, sigma, seed, dt_multiple=k)
    truth, xr, yr = analytic_truth(case, dt)
    return degrade_analytic(truth, xr, yr, n, window_h, sigma, seed, refine)


def available_cases(requested: Iterable[str]) -> list[str]:
    """Drop ``cfd`` (with a message) when the CFD export is missing."""
    out = []
    for c in requested:
        if c == "cfd" and not cfd_available():
            print(f"[skip] CFD data not found at {DEFAULT_CFD_NPZ} - run the CFD case first.")
            continue
        out.append(c)
    return out


def core_region(case: str) -> tuple[float, float, float, float] | None:
    """Search region for the vortex-core minimum (None = whole window)."""
    return CFD_CORE_REGION if case == "cfd" else None


def has_core(case: str) -> bool:
    """True for cases whose window contains a vortex core."""
    return case in ("tg_core", "cfd")


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------
def evaluate(fld: FlowField, case: str, with_core: bool = False, **kwargs: Any) -> dict[str, Any]:
    """Reconstruct and return the metrics as a flat dict (one CSV row)."""
    rec = reconstruct(fld, **kwargs)
    row: dict[str, Any] = standard_metrics(rec.p, rec.p_ref, fld.rho, fld.U_ref)
    d = rec.diagnostics
    row.update(time_s=d["time_total_s"], time_integration_s=d["time_integration_s"],
               n_points=fld.nx, h=fld.dx)
    for key in ("compatibility_relative", "path_rms_difference", "direction_spread_rms",
                "fallback_fraction"):
        if key in d:
            row[key] = d[key]
    if with_core:
        row.update(core_pressure(rec.p, rec.p_ref, fld.x, fld.y, fld.rho, fld.U_ref,
                                 region=core_region(case)))
    return row


def mean_std(rows: list[dict], value: str, **match: Any) -> tuple[float, float, int]:
    """Mean and sample std of ``value`` over rows matching all ``match`` items."""
    vals = [float(r[value]) for r in rows
            if all(_eq(r.get(k), v) for k, v in match.items()) and _isnum(r.get(value))]
    if not vals:
        return float("nan"), float("nan"), 0
    a = np.array(vals)
    return float(a.mean()), float(a.std(ddof=1)) if a.size > 1 else 0.0, int(a.size)


def _isnum(v: Any) -> bool:
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _eq(a: Any, b: Any) -> bool:
    """Compare CSV values that may be strings or numbers."""
    if _isnum(a) and _isnum(b):
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)
    return str(a) == str(b)


# --------------------------------------------------------------------------
# CSV I/O and command line
# --------------------------------------------------------------------------
def write_csv(rows: list[dict], path: Path) -> Path:
    """Write rows (dicts) to CSV; the header is the union of all keys."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k, "")) for k in keys})
    print(f"[csv] wrote {len(rows)} rows -> {path}")
    return path


def _fmt(v: Any) -> Any:
    return f"{v:.6g}" if isinstance(v, float) else v


def read_csv(path: Path) -> list[dict]:
    """Read a CSV written by :func:`write_csv` (numbers stay strings)."""
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def base_parser(description: str, default_cases: tuple[str, ...]) -> argparse.ArgumentParser:
    """Common command-line options of every study script."""
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--cases", nargs="+", default=list(default_cases),
                    help=f"case keys (default: {' '.join(default_cases)})")
    ap.add_argument("--seeds", type=int, default=10, help="noise realisations per point (default 10)")
    ap.add_argument("--quick", action="store_true", help="fewer seeds / coarser grids for a fast check")
    ap.add_argument("--outdir", type=Path, default=RESULTS_DIR, help="output directory")
    return ap


def seeds_for(args: argparse.Namespace) -> list[int]:
    """Seed list: 0..seeds-1 (3 seeds in quick mode)."""
    return list(range(3 if args.quick else args.seeds))
