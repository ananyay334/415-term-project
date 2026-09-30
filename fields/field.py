# piv_pressure/fields/field.py
"""Common data container used by every field source (analytic, CFD, JHTDB).

Every loader in :mod:`fields` returns a :class:`FlowField`.  The rest of the
pipeline (degradation, derivatives, pressure integration, metrics) only ever
talks to this interface, so a new data source only has to fill it in.

Array convention (used everywhere in the project)
-------------------------------------------------
* ``x`` has shape ``(nx,)`` and ``y`` has shape ``(ny,)``; both are uniformly
  spaced and increasing.
* Velocity and pressure arrays have shape ``(nt, ny, nx)``: index ``[n, j, i]``
  is the value at time ``t[n]``, height ``y[j]`` and streamwise position
  ``x[i]``.  Axis ``-1`` is therefore the x direction and axis ``-2`` the y
  direction.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np


@dataclass
class FlowField:
    """Planar velocity snapshots plus (optionally) the true pressure.

    Attributes
    ----------
    x, y : 1-D arrays of grid coordinates (uniform, increasing).
    t : 1-D array of snapshot times (uniform spacing if more than one).
    u, v : velocity components, shape ``(nt, ny, nx)``.
    rho : density (the project uses rho = 1, i.e. kinematic pressure).
    mu : dynamic viscosity (mu = rho * nu).
    p_ref : reference ("true") pressure, same shape as ``u``, or ``None``.
    U_ref : velocity scale used for noise levels and for epsilon_max.
    name : short label used in CSV files and plot legends.
    meta : free-form metadata (Reynolds number, source file, ...).
    """

    x: np.ndarray
    y: np.ndarray
    t: np.ndarray
    u: np.ndarray
    v: np.ndarray
    rho: float
    mu: float
    p_ref: np.ndarray | None = None
    U_ref: float = 1.0
    name: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        self.t = np.atleast_1d(np.asarray(self.t, dtype=float))
        self.u = np.asarray(self.u, dtype=float)
        self.v = np.asarray(self.v, dtype=float)
        if self.p_ref is not None:
            self.p_ref = np.asarray(self.p_ref, dtype=float)
        self.validate()

    # ------------------------------------------------------------------ checks
    def validate(self) -> None:
        """Raise ``ValueError`` if shapes or grid spacing are inconsistent."""
        shape = (self.t.size, self.y.size, self.x.size)
        for label, arr in (("u", self.u), ("v", self.v), ("p_ref", self.p_ref)):
            if arr is not None and arr.shape != shape:
                raise ValueError(f"{label} has shape {arr.shape}, expected {shape}")
        for label, c in (("x", self.x), ("y", self.y)):
            if c.size < 3:
                raise ValueError(f"need at least 3 points in {label}")
            d = np.diff(c)
            if np.any(d <= 0) or not np.allclose(d, d[0], rtol=1e-6, atol=0.0):
                raise ValueError(f"{label} must be uniform and increasing")
        if self.t.size > 1:
            d = np.diff(self.t)
            if np.any(d <= 0) or not np.allclose(d, d[0], rtol=1e-6, atol=0.0):
                raise ValueError("t must be uniform and increasing")

    # -------------------------------------------------------------- properties
    @property
    def nt(self) -> int:
        """Number of snapshots."""
        return int(self.t.size)

    @property
    def ny(self) -> int:
        """Number of grid points in y."""
        return int(self.y.size)

    @property
    def nx(self) -> int:
        """Number of grid points in x."""
        return int(self.x.size)

    @property
    def dx(self) -> float:
        """Grid spacing in x."""
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        """Grid spacing in y."""
        return float(self.y[1] - self.y[0])

    @property
    def dt(self) -> float:
        """Time between snapshots (``nan`` for a single steady snapshot)."""
        return float(self.t[1] - self.t[0]) if self.nt > 1 else float("nan")

    @property
    def steady(self) -> bool:
        """True when the field is a single steady snapshot (du/dt = 0)."""
        return self.nt == 1

    @property
    def nu(self) -> float:
        """Kinematic viscosity nu = mu / rho."""
        return self.mu / self.rho

    # ------------------------------------------------------------- utilities
    def copy_with(self, **changes: Any) -> "FlowField":
        """Return a copy with some attributes replaced (arrays are not copied)."""
        return replace(self, **changes)

    def subset_time(self, indices: list[int] | np.ndarray) -> "FlowField":
        """Return a field containing only the snapshots ``indices``."""
        idx = np.asarray(indices, dtype=int)
        return self.copy_with(
            t=self.t[idx],
            u=self.u[idx],
            v=self.v[idx],
            p_ref=None if self.p_ref is None else self.p_ref[idx],
        )

    def grids(self) -> tuple[np.ndarray, np.ndarray]:
        """Return 2-D coordinate arrays ``X, Y`` of shape ``(ny, nx)``."""
        return np.meshgrid(self.x, self.y, indexing="xy")
