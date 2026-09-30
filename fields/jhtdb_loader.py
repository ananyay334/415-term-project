# piv_pressure/fields/jhtdb_loader.py
"""Optional: 2-D slices of 3-D turbulence from the Johns Hopkins Turbulence
Database (JHTDB).  Nothing else in the project depends on this module.

Supported inputs
----------------
* ``.npz`` with keys ``x, y, t, u, v, p`` (shape ``[nt, ny, nx]``) and
  optionally ``w`` - the simplest route: convert whatever the JHTDB tools give
  you into this format.
* ``.h5`` cutout files (needs the optional package ``h5py``).  The dataset
  names assumed here - ``Velocity_NNNN`` with shape (nz, ny, nx, 3),
  ``Pressure_NNNN`` with shape (nz, ny, nx, 1), coordinates ``xcoor``,
  ``ycoor`` - follow the JHTDB cutout service [verify against your file with
  ``h5ls``].

Why a 2-D slice of 3-D turbulence is harder (the out-of-plane error)
-------------------------------------------------------------------
Planar PIV gives u, v on one plane.  The full x-momentum equation is::

    dp/dx = -rho (u_t + u u_x + v u_y + w u_z) + mu (u_xx + u_yy + u_zz)

The 2-D pipeline has no w, u_z or u_zz, so it drops
``M_x = -rho w u_z + mu u_zz`` (and similarly M_y).  In 3-D turbulence these
terms are as large as the in-plane ones, so the reconstructed pressure has an
O(1) model error that no grid refinement or denoising removes.  With a
three-plane cutout (z - dz, z, z + dz) :func:`out_of_plane_terms` computes
M_x, M_y so the size of this error can be quantified.

Isotropic-1024 dataset parameters (for reference only) [verify on the JHTDB
website]: nu = 0.000185, periodic domain 2 pi, stored time step 0.002.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from fields.field import FlowField

#: Reference parameters of JHTDB "isotropic1024coarse" [verify].
JHTDB_ISOTROPIC_NU = 0.000185
JHTDB_ISOTROPIC_DT = 0.002


def load_jhtdb_npz(path: str | Path, rho: float = 1.0, nu: float = JHTDB_ISOTROPIC_NU) -> FlowField:
    """Load a converted slice ``.npz`` (keys x, y, t, u, v, p[, w])."""
    with np.load(path, allow_pickle=False) as d:
        meta = {"case": "jhtdb", "source_file": str(path)}
        if "w" in d:
            meta["w"] = d["w"]
        return FlowField(x=d["x"], y=d["y"], t=d["t"], u=d["u"], v=d["v"], rho=rho,
                         mu=rho * nu, p_ref=rho * d["p"], U_ref=float(np.std(d["u"]) or 1.0),
                         name="jhtdb_slice", meta=meta)


def load_jhtdb_h5(
    path: str | Path,
    z_index: int = 0,
    dt: float = JHTDB_ISOTROPIC_DT,
    rho: float = 1.0,
    nu: float = JHTDB_ISOTROPIC_NU,
) -> FlowField:
    """Load plane ``z_index`` of a JHTDB cutout HDF5 file (needs ``h5py``)."""
    try:
        import h5py  # optional dependency
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("reading .h5 needs `pip install h5py`, or convert to .npz") from exc
    with h5py.File(path, "r") as h5:
        vel_names = sorted(k for k in h5 if re.match(r"Velocity_\d+", k))
        prs_names = sorted(k for k in h5 if re.match(r"Pressure_\d+", k))
        if not vel_names or len(vel_names) != len(prs_names):
            raise ValueError("expected matching Velocity_NNNN / Pressure_NNNN datasets")
        vel = np.stack([h5[k][z_index] for k in vel_names])     # (nt, ny, nx, 3)
        prs = np.stack([h5[k][z_index, ..., 0] for k in prs_names])
        x, y = np.asarray(h5["xcoor"]), np.asarray(h5["ycoor"])
    t = dt * np.arange(len(vel_names))
    return FlowField(x=x, y=y, t=t, u=vel[..., 0], v=vel[..., 1], rho=rho, mu=rho * nu,
                     p_ref=rho * prs, U_ref=float(np.std(vel[..., 0])),
                     name="jhtdb_slice", meta={"case": "jhtdb", "w": vel[..., 2]})


def load_jhtdb_slice(path: str | Path, **kwargs) -> FlowField:
    """Dispatch on the file suffix (``.npz`` or ``.h5``/``.hdf5``)."""
    suffix = Path(path).suffix.lower()
    if suffix == ".npz":
        return load_jhtdb_npz(path, **kwargs)
    if suffix in (".h5", ".hdf5"):
        return load_jhtdb_h5(path, **kwargs)
    raise ValueError(f"unsupported JHTDB file type {suffix!r}")


def out_of_plane_terms(
    u3: np.ndarray, v3: np.ndarray, w3: np.ndarray, dz: float, rho: float = 1.0, mu: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Terms of G missing from a planar reconstruction, on the middle plane.

    Parameters
    ----------
    u3, v3, w3 : arrays of shape (3, ny, nx) on planes z - dz, z, z + dz.
    dz : plane spacing.

    Returns
    -------
    M_x, M_y with M = -rho w du/dz + mu d2u/dz2 (central differences).
    """
    def dz1(f: np.ndarray) -> np.ndarray:
        return (f[2] - f[0]) / (2.0 * dz)

    def dz2(f: np.ndarray) -> np.ndarray:
        return (f[2] - 2.0 * f[1] + f[0]) / dz**2

    w = w3[1]
    return (-rho * w * dz1(u3) + mu * dz2(u3), -rho * w * dz1(v3) + mu * dz2(v3))
