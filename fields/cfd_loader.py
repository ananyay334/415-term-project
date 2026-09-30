# piv_pressure/fields/cfd_loader.py
"""Load the CFD cylinder-wake window (exported by ``cfd/export_window.py``).

File format (one ``.npz``)
--------------------------
x (nx,), y (ny,), t (nt,)          uniform grids (x/D, y/D, t U/D)
u, v, p (nt, ny, nx)               velocity and *kinematic* pressure p/rho
rho, nu, U_inf, D                  scalars (rho = 1 in the non-dimensional set-up)
meta_json                          JSON string with provenance (mesh, solver, ...)
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from fields.field import FlowField

PROJECT_ROOT = Path(__file__).resolve().parents[1]
#: Default location of the exported window; the environment variable
#: PIV_CFD_NPZ overrides it (e.g. to keep large data outside the project).
DEFAULT_CFD_NPZ = Path(os.environ.get(
    "PIV_CFD_NPZ", PROJECT_ROOT / "results" / "cfd" / "cylinder_re100_window.npz"))


def cfd_available(path: str | Path = DEFAULT_CFD_NPZ) -> bool:
    """True if the exported CFD window file exists."""
    return Path(path).is_file()


def save_cfd_npz(
    path: str | Path,
    x: np.ndarray,
    y: np.ndarray,
    t: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    p: np.ndarray,
    rho: float = 1.0,
    nu: float = 0.01,
    U_inf: float = 1.0,
    D: float = 1.0,
    meta: dict | None = None,
) -> Path:
    """Write the common CFD-window format (used by the export script)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, x=x, y=y, t=t, u=u, v=v, p=p, rho=rho, nu=nu, U_inf=U_inf, D=D,
        meta_json=json.dumps(meta or {}),
    )
    return path


def load_cfd(path: str | Path = DEFAULT_CFD_NPZ) -> FlowField:
    """Load the CFD window into the common :class:`FlowField` interface.

    Raises ``FileNotFoundError`` with instructions if the export is missing.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found. Run the OpenFOAM cases and then "
            "`python cfd/export_window.py` (see README, section 'CFD case')."
        )
    with np.load(path, allow_pickle=False) as d:
        rho = float(d["rho"])
        meta = json.loads(str(d["meta_json"])) if "meta_json" in d else {}
        meta.update(case="cfd_cylinder", source_file=str(path), D=float(d["D"]))
        return FlowField(
            x=d["x"], y=d["y"], t=d["t"], u=d["u"], v=d["v"], rho=rho,
            mu=rho * float(d["nu"]), p_ref=d["p"], U_ref=float(d["U_inf"]),
            name="cfd_cylinder_Re100", meta=meta,
        )
