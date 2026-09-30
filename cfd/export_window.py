# piv_pressure/cfd/export_window.py
"""Convert the sampled PIV window of the fine OpenFOAM case into one .npz.

Input (written at run time by the ``sampleWindow`` function object)::

    cfd/cases/fine/postProcessing/sampleWindow/<time>/window_p.xy
    cfd/cases/fine/postProcessing/sampleWindow/<time>/window_U.xy

Raw set files have columns ``x y z <field components>``.  Depending on the
OpenFOAM version the scalar and vector fields are in separate files or in
one combined file such as ``window_p_U.xy`` [verify]; both layouts are
handled: the field names are read from the file name and the columns are
assigned in that order (p: 1 column, U: 3 columns).  Sample points are
mapped back to the 81 x 81 grid by their coordinates, so the point order in
the file does not matter.

Fallback for ANSYS Fluent (see FLUENT_SETTINGS.md): ``--fluent-dir`` reads
one ASCII export per time step (columns x-coordinate, y-coordinate,
x-velocity, y-velocity, pressure [verify header names]) and interpolates it
linearly onto the window grid.

Output: results/cfd/cylinder_re100_window.npz (see fields/cfd_loader.py).

Usage:  python cfd/export_window.py [--case cfd/cases/fine] [--fluent-dir DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fields.cfd_loader import DEFAULT_CFD_NPZ, save_cfd_npz  # noqa: E402

FIELD_WIDTH = {"p": 1, "U": 3}


def window_grid(info: dict) -> tuple[np.ndarray, np.ndarray]:
    """Grid coordinates from the ``window`` entry of mesh_info.json."""
    w = info["window"]
    return np.linspace(w["x0"], w["x1"], w["nx"]), np.linspace(w["y0"], w["y1"], w["ny"])


def fields_in_name(fname: str, set_name: str) -> list[str]:
    """Field names encoded in a set file name, e.g. window_p_U.xy -> [p, U]."""
    stem = Path(fname).name.split(".")[0]
    if not stem.startswith(set_name + "_"):
        return []
    return [f for f in stem[len(set_name) + 1:].split("_") if f]


def read_set_file(path: Path) -> np.ndarray:
    """Numeric table of a raw (whitespace) or csv set file; header lines skipped."""
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or re.match(r"^[A-Za-z]", line):
            continue
        rows.append([float(t) for t in re.split(r"[,\s]+", line)])
    return np.array(rows)


def to_grid_indices(xs: np.ndarray, ys: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest grid indices of sample coordinates (must lie on the grid)."""
    i = np.rint((xs - x[0]) / (x[1] - x[0])).astype(int)
    j = np.rint((ys - y[0]) / (y[1] - y[0])).astype(int)
    if (np.abs(x[np.clip(i, 0, x.size - 1)] - xs).max() > 1e-4
            or np.abs(y[np.clip(j, 0, y.size - 1)] - ys).max() > 1e-4):
        raise ValueError("sample points do not match the window grid")
    return j, i


def read_time_dir(tdir: Path, set_name: str, x: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray]:
    """u, v, p on the (ny, nx) grid from one sampling time directory."""
    out = {k: np.full((y.size, x.size), np.nan) for k in ("u", "v", "p")}
    for f in sorted(tdir.iterdir()):
        names = fields_in_name(f.name, set_name)
        if not names or any(n not in FIELD_WIDTH for n in names):
            continue
        data = read_set_file(f)
        expected = 3 + sum(FIELD_WIDTH[n] for n in names)
        if data.ndim != 2 or data.shape[1] != expected:
            raise ValueError(f"{f}: expected {expected} columns, found {data.shape}")
        j, i = to_grid_indices(data[:, 0], data[:, 1], x, y)
        col = 3
        for n in names:
            if n == "p":
                out["p"][j, i] = data[:, col]
            else:
                out["u"][j, i] = data[:, col]
                out["v"][j, i] = data[:, col + 1]
            col += FIELD_WIDTH[n]
    return out


def load_openfoam_samples(case: Path, set_name: str = "window", fo_name: str = "sampleWindow"):
    """Read every sampling time of the fine case; returns (x, y, t, u, v, p, info)."""
    info = json.loads((case / "mesh_info.json").read_text())
    x, y = window_grid(info)
    base = case / "postProcessing" / fo_name
    if not base.is_dir():
        raise FileNotFoundError(f"{base} not found - has the fine case been run?")
    tdirs = sorted((d for d in base.iterdir() if d.is_dir() and _is_number(d.name)),
                   key=lambda d: float(d.name))
    t, U, V, P = [], [], [], []
    for d in tdirs:
        g = read_time_dir(d, set_name, x, y)
        if any(np.isnan(a).any() for a in g.values()):
            print(f"[warn] incomplete samples in {d} - skipped")
            continue
        t.append(float(d.name))
        U.append(g["u"]); V.append(g["v"]); P.append(g["p"])
    if len(t) < 3:
        raise RuntimeError("fewer than 3 complete sample times found")
    return x, y, np.array(t), np.array(U), np.array(V), np.array(P), info


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def uniform_in_time(t: np.ndarray, *fields: np.ndarray, dt: float | None = None, tol: float = 1e-3):
    """Return data on a uniform time grid (linear interpolation if needed)."""
    dt = float(np.median(np.diff(t))) if dt is None else dt
    if np.all(np.abs(np.diff(t) - dt) <= tol * dt):
        return t[0] + dt * np.arange(t.size), fields  # already uniform: remove round-off
    print(f"[warn] non-uniform sample times; interpolating to dt = {dt:g}")
    tu = t[0] + dt * np.arange(int(np.floor((t[-1] - t[0]) / dt + 1e-9)) + 1)
    out = []
    for f in fields:
        flat = f.reshape(f.shape[0], -1)
        g = np.empty((tu.size, flat.shape[1]))
        for k in range(flat.shape[1]):
            g[:, k] = np.interp(tu, t, flat[:, k])
        out.append(g.reshape((tu.size,) + f.shape[1:]))
    return tu, tuple(out)


def load_fluent_exports(folder: Path, x: np.ndarray, y: np.ndarray, dt: float):
    """Read Fluent ASCII exports (one file per time, sorted by name) onto the grid."""
    from scipy.interpolate import griddata

    files = sorted(p for p in folder.iterdir() if p.is_file())
    X, Y = np.meshgrid(x, y)
    U, V, P = [], [], []
    for f in files:
        lines = f.read_text().splitlines()
        header = [h.strip().lower() for h in re.split(r"[,\s]+", lines[0].strip()) if h.strip()]
        data = np.array([[float(v) for v in re.split(r"[,\s]+", ln.strip()) if v]
                         for ln in lines[1:] if ln.strip()])
        col = {name: header.index(name) for name in
               ("x-coordinate", "y-coordinate", "x-velocity", "y-velocity", "pressure")}
        pts = data[:, [col["x-coordinate"], col["y-coordinate"]]]
        for store, key in ((U, "x-velocity"), (V, "y-velocity"), (P, "pressure")):
            store.append(griddata(pts, data[:, col[key]], (X, Y), method="linear"))
    t = dt * np.arange(len(files))
    return t, np.array(U), np.array(V), np.array(P)


def main(argv: list[str] | None = None) -> Path:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", type=Path, default=ROOT / "cfd" / "cases" / "fine")
    ap.add_argument("--out", type=Path, default=DEFAULT_CFD_NPZ)
    ap.add_argument("--fluent-dir", type=Path, default=None,
                    help="folder of Fluent ASCII exports (one per time step) instead of OpenFOAM")
    ap.add_argument("--fluent-dt", type=float, default=0.05, help="time between Fluent exports")
    args = ap.parse_args(argv)
    if args.fluent_dir is not None:
        info = json.loads((ROOT / "cfd" / "cases" / "fine" / "mesh_info.json").read_text())
        x, y = window_grid(info)
        t, U, V, P = load_fluent_exports(args.fluent_dir, x, y, args.fluent_dt)
        source = f"fluent:{args.fluent_dir}"
    else:
        x, y, t, U, V, P, info = load_openfoam_samples(args.case)
        source = f"openfoam:{args.case}"
    t, (U, V, P) = uniform_in_time(t, U, V, P)
    meta = {"source": source, "mesh": info.get("name"), "cells": info.get("cells"),
            "Re": info.get("Re"), "nu": info.get("nu"), "D": info.get("D")}
    path = save_cfd_npz(args.out, x, y, t, U, V, P, rho=1.0, nu=float(info.get("nu", 0.01)),
                        U_inf=float(info.get("U_inf", 1.0)), D=float(info.get("D", 1.0)), meta=meta)
    print(f"[export] {t.size} snapshots, grid {y.size} x {x.size}, dt = {t[1] - t[0]:g} -> {path}")
    return path


if __name__ == "__main__":
    main()
