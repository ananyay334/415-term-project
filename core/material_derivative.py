# piv_pressure/core/material_derivative.py
"""Material acceleration Du/Dt and the pressure gradient G from PIV data.

Navier-Stokes (incompressible, constant rho and mu) rearranged for the
pressure gradient (van Oudheusden 2013, Eq. 1)::

    grad p = G = -rho Du/Dt + mu lap(u),     Du/Dt = du/dt + (u . grad) u

Two estimates of Du/Dt are provided:

* **Eulerian** - du/dt by the central difference (u^{n+1} - u^{n-1})/(2 dt)
  and the convective term u du/dx + v du/dy by second-order central
  differences (one-sided second-order at the edges).
* **Lagrangian (pseudo-tracking)** - follow a straight path through the node
  using the local velocity: x_f = x + u dt, x_b = x - u dt, then
  Du/Dt = [u(x_f, t^{n+1}) - u(x_b, t^{n-1})] / (2 dt).  The velocity at the
  off-grid points is interpolated with ``RegularGridInterpolator`` (cubic if
  SciPy supports it, else linear).  Where x_f or x_b leaves the measurement
  domain the Eulerian value is used instead (de Kat & van Oudheusden 2012;
  van Gent et al. 2017 discuss both approaches).

For a steady field (one snapshot) du/dt = 0.  The Lagrangian option then
tracks through the same snapshot with a pseudo time step ``dt_track``.

``order`` (default 2) selects the *spatial* scheme of the convective and
viscous terms (see core/derivatives.py); ``order=1`` is only used to
demonstrate first- vs second-order grid convergence.  The time difference
stays central in both cases because the grid-convergence study varies dx at
fixed dt.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import spsolve

from core.derivatives import gradient, laplacian, time_derivative
from fields.field import FlowField


@dataclass
class PressureGradient:
    """Pressure gradient G = (Gx, Gy) at one time level plus diagnostics."""

    Gx: np.ndarray
    Gy: np.ndarray
    Dux: np.ndarray
    Duy: np.ndarray
    time_index: int
    info: dict = field(default_factory=dict)


def _check_time_index(fld: FlowField, n: int | None) -> int:
    """Return a valid time index (default: the central snapshot)."""
    if fld.steady:
        return 0
    if n is None:
        n = fld.nt // 2
    if not 0 <= n < fld.nt:
        raise IndexError(f"time index {n} out of range 0..{fld.nt - 1}")
    return n


def convective_acceleration(
    u: np.ndarray, v: np.ndarray, dx: float, dy: float, order: int = 2
) -> tuple[np.ndarray, np.ndarray]:
    """Return (u.grad)u = (u du/dx + v du/dy, u dv/dx + v dv/dy)."""
    ux, uy = gradient(u, dx, dy, order)
    vx, vy = gradient(v, dx, dy, order)
    return u * ux + v * uy, u * vx + v * vy


def eulerian_material_derivative(
    fld: FlowField, n: int | None = None, order: int = 2
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Eulerian Du/Dt at snapshot ``n`` (default: central snapshot).

    Returns ``(Dux, Duy, info)``.  The time derivative is central for interior
    snapshots and one-sided second-order at the first/last snapshot; the
    convective term uses spatial differences of the given ``order``.
    """
    n = _check_time_index(fld, n)
    ax, ay = convective_acceleration(fld.u[n], fld.v[n], fld.dx, fld.dy, order)
    if fld.steady:
        return ax, ay, {"time_scheme": "steady"}
    if fld.nt < 3:
        raise ValueError("unsteady Eulerian scheme needs at least 3 snapshots")
    ut = time_derivative(fld.u, fld.dt, n)
    vt = time_derivative(fld.v, fld.dt, n)
    scheme = "central" if 0 < n < fld.nt - 1 else "one-sided"
    return ut + ax, vt + ay, {"time_scheme": scheme}


def _make_interpolator(
    x: np.ndarray, y: np.ndarray, f: np.ndarray, method: str
) -> tuple[RegularGridInterpolator, str]:
    """Build a RegularGridInterpolator over (y, x); fall back to linear.

    ``f`` has shape (ny, nx) or (ny, nx, k) (k fields interpolated together).
    ``method="auto"`` tries cubic first (SciPy >= 1.10 and >= 4 points per
    axis) and falls back to linear if cubic is not available.

    Pitfall: from SciPy 1.13 the cubic spline coefficients are found with an
    *iterative* solver whose default tolerance (atol = 1e-6) limits the
    interpolation error to ~1e-5 however fine the grid is.  The Lagrangian
    formula divides that error by 2 dt, so the direct solver ``spsolve`` is
    requested whenever SciPy accepts the ``solver`` argument.
    """
    tries = ["cubic", "linear"] if method == "auto" else [method]
    last_error: Exception | None = None
    for m in tries:
        kwargs_list = [{"solver": spsolve}, {}] if m == "cubic" else [{}]
        for kwargs in kwargs_list:
            try:
                interp = RegularGridInterpolator(
                    (y, x), f, method=m, bounds_error=False, fill_value=np.nan, **kwargs
                )
                interp(np.array([[y[0], x[0]]]))  # force any lazy set-up errors now
                return interp, m
            except TypeError as exc:  # SciPy < 1.13: no 'solver' keyword
                last_error = exc
            except (ValueError, NotImplementedError) as exc:  # cubic unsupported
                last_error = exc
                break
    raise RuntimeError(f"could not build interpolator: {last_error}")


def lagrangian_material_derivative(
    fld: FlowField,
    n: int | None = None,
    interpolation: str = "auto",
    dt_track: float | None = None,
    order: int = 2,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Pseudo-tracking (Lagrangian) Du/Dt at snapshot ``n``.

    Parameters
    ----------
    fld : the flow field.
    n : snapshot index; needs ``0 < n < nt-1`` for unsteady data.
    interpolation : ``"auto"`` (cubic, else linear), ``"cubic"`` or ``"linear"``.
    dt_track : pseudo time step for *steady* data (default h / max|u|, so
        the tracked points move at most one grid spacing).
    order : spatial order of the Eulerian fall-back near the edges.

    Returns
    -------
    Dux, Duy, info : ``info`` holds the interpolation method used and the
        fraction of nodes that fell back to the Eulerian value.
    """
    n = _check_time_index(fld, n)
    if fld.steady:
        n_prev = n_next = 0
        speed_max = float(np.max(np.hypot(fld.u[0], fld.v[0])))
        dt = dt_track if dt_track is not None else min(fld.dx, fld.dy) / max(speed_max, 1e-12)
    else:
        if not 0 < n < fld.nt - 1:
            raise ValueError("Lagrangian scheme needs snapshots n-1 and n+1")
        n_prev, n_next = n - 1, n + 1
        dt = fld.dt

    X, Y = fld.grids()
    u0, v0 = fld.u[n], fld.v[n]
    # Straight-line particle positions one time step ahead and behind.
    xf, yf = X + u0 * dt, Y + v0 * dt
    xb, yb = X - u0 * dt, Y - v0 * dt

    # One interpolator per time level; u and v are stacked on a last axis.
    i_next, method = _make_interpolator(
        fld.x, fld.y, np.stack([fld.u[n_next], fld.v[n_next]], axis=-1), interpolation)
    i_prev, _ = _make_interpolator(
        fld.x, fld.y, np.stack([fld.u[n_prev], fld.v[n_prev]], axis=-1), method)

    pts_f = np.column_stack([yf.ravel(), xf.ravel()])  # interpolator order is (y, x)
    pts_b = np.column_stack([yb.ravel(), xb.ravel()])
    vel_f = i_next(pts_f)                               # shape (N, 2): (u, v)
    vel_b = i_prev(pts_b)
    uf, vf = vel_f[:, 0].reshape(X.shape), vel_f[:, 1].reshape(X.shape)
    ub, vb = vel_b[:, 0].reshape(X.shape), vel_b[:, 1].reshape(X.shape)

    Dux = (uf - ub) / (2.0 * dt)
    Duy = (vf - vb) / (2.0 * dt)

    # Points whose forward or backward position left the domain -> NaN from
    # the interpolator.  Replace them by the Eulerian estimate.
    outside = ~(np.isfinite(Dux) & np.isfinite(Duy))
    if np.any(outside):
        Ex, Ey, _ = eulerian_material_derivative(fld, n, order)
        Dux = np.where(outside, Ex, Dux)
        Duy = np.where(outside, Ey, Duy)
    info = {
        "interpolation": method,
        "fallback_fraction": float(outside.mean()),
        "dt_track": float(dt),
    }
    return Dux, Duy, info


def pressure_gradient(
    fld: FlowField,
    n: int | None = None,
    scheme: str = "eulerian",
    include_viscous: bool = True,
    interpolation: str = "auto",
    order: int = 2,
) -> PressureGradient:
    """Compute G = -rho Du/Dt + mu lap(u) at snapshot ``n``.

    Parameters
    ----------
    fld : the flow field.
    n : snapshot index (default: central snapshot).
    scheme : ``"eulerian"`` or ``"lagrangian"``.
    include_viscous : drop the mu lap(u) term if False (used only to show
        that the viscous term matters for Kovasznay flow).
    interpolation : interpolation used by the Lagrangian scheme.
    order : spatial order of the finite differences (2 = default, 1 = first
        order, for the grid-convergence demonstration only).
    """
    n = _check_time_index(fld, n)
    if scheme == "eulerian":
        Dux, Duy, info = eulerian_material_derivative(fld, n, order)
    elif scheme == "lagrangian":
        Dux, Duy, info = lagrangian_material_derivative(fld, n, interpolation, order=order)
    else:
        raise ValueError(f"unknown material-derivative scheme {scheme!r}")
    Gx = -fld.rho * Dux
    Gy = -fld.rho * Duy
    if include_viscous and fld.mu != 0.0:
        Gx = Gx + fld.mu * laplacian(fld.u[n], fld.dx, fld.dy, order)
        Gy = Gy + fld.mu * laplacian(fld.v[n], fld.dx, fld.dy, order)
    info = dict(info, scheme=scheme, order=order)
    return PressureGradient(Gx=Gx, Gy=Gy, Dux=Dux, Duy=Duy, time_index=n, info=info)
