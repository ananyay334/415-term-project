# piv_pressure/fields/analytic.py
"""Exact Navier-Stokes solutions used as validation cases.

1. Taylor-Green vortex (2-D, decaying, unsteady)::

       u =  cos x sin y F(t),   v = -sin x cos y F(t),   F = exp(-2 nu t)
       p = -(rho/4) (cos 2x + cos 2y) F(t)^2

   Each velocity component is an eigenfunction of the Laplacian
   (lap u = -2u), so du/dt = nu lap u and the pressure gradient balances
   the convective term only.

   Any amplitude A works: u -> A u, p -> A^2 p is still exact (the viscous
   decay does not depend on A).

   *Convected* variant (added for the time-step study, see EXPLANATIONS.md):
   the Navier-Stokes equations are Galilean invariant, so
   ``u'(x, t) = c + u(x - c t, t)`` and ``p'(x, t) = p(x - c t, t)`` is also
   an exact solution.  With c = (1, 0) the vortex pattern is advected
   through the window, so du/dt is O(1) (as in a real wake) instead of
   O(nu).

2. Kovasznay flow (2-D, steady; Kovasznay 1948)::

       lambda = Re/2 - sqrt(Re^2/4 + 4 pi^2)
       u = 1 - exp(lambda x) cos(2 pi y)
       v = (lambda / 2 pi) exp(lambda x) sin(2 pi y)
       p = (rho/2) (1 - exp(2 lambda x))

   with nu = 1/Re (U = 1, L = 1).
"""
from __future__ import annotations

import numpy as np

from fields.field import FlowField

#: Non-periodic Taylor-Green window required by the spec (both axes).
TG_WINDOW: tuple[float, float] = (0.3, 2.5)
#: Window centred on the vortex core at (pi, pi); used for vortex-core metrics
#: because TG_WINDOW contains no vortex core (see EXPLANATIONS.md).
TG_CORE_WINDOW: tuple[float, float] = (np.pi - 1.1, np.pi + 1.1)
#: Location of the Taylor-Green vortex core used for core metrics (t = 0).
TG_CORE: tuple[float, float] = (np.pi, np.pi)
#: Kovasznay domain required by the spec.
KOV_X: tuple[float, float] = (-0.5, 1.5)
KOV_Y: tuple[float, float] = (-0.5, 1.5)


# --------------------------------------------------------------------------
# Taylor-Green
# --------------------------------------------------------------------------
def taylor_green_uvp(
    X: np.ndarray,
    Y: np.ndarray,
    t: float,
    nu: float = 0.01,
    rho: float = 1.0,
    advect: tuple[float, float] = (0.0, 0.0),
    amplitude: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact Taylor-Green velocity and pressure at time ``t``.

    Parameters
    ----------
    X, Y : coordinate arrays (any matching shape).
    t : time.
    nu, rho : kinematic viscosity and density.
    advect : uniform advection velocity c = (cx, cy) of the Galilean-shifted
        solution (``(0, 0)`` gives the classical stationary vortex array).
    amplitude : velocity amplitude A of the vortex array (spec: A = 1).

    Returns
    -------
    u, v, p : arrays with the shape of ``X``.
    """
    cx, cy = advect
    xi = X - cx * t          # coordinates in the frame moving with the pattern
    eta = Y - cy * t
    F = amplitude * np.exp(-2.0 * nu * t)
    u = cx + np.cos(xi) * np.sin(eta) * F
    v = cy - np.sin(xi) * np.cos(eta) * F
    p = -0.25 * rho * (np.cos(2.0 * xi) + np.cos(2.0 * eta)) * F**2
    return u, v, p


def taylor_green_grad_p(
    X: np.ndarray,
    Y: np.ndarray,
    t: float,
    nu: float = 0.01,
    rho: float = 1.0,
    advect: tuple[float, float] = (0.0, 0.0),
    amplitude: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact pressure gradient (dp/dx, dp/dy) of the Taylor-Green vortex.

    Differentiating ``p = -(rho/4)(cos 2 xi + cos 2 eta) F^2`` gives
    ``dp/dx = (rho/2) sin(2 xi) F^2`` and ``dp/dy = (rho/2) sin(2 eta) F^2``.
    """
    cx, cy = advect
    xi = X - cx * t
    eta = Y - cy * t
    F2 = amplitude**2 * np.exp(-4.0 * nu * t)
    return 0.5 * rho * np.sin(2.0 * xi) * F2, 0.5 * rho * np.sin(2.0 * eta) * F2


def taylor_green_on_grid(
    x: np.ndarray,
    y: np.ndarray,
    t: np.ndarray,
    nu: float = 0.01,
    rho: float = 1.0,
    advect: tuple[float, float] = (0.0, 0.0),
    name: str = "taylor_green",
    amplitude: float = 1.0,
) -> FlowField:
    """Sample the Taylor-Green solution on given 1-D grids and times."""
    X, Y = np.meshgrid(x, y, indexing="xy")
    t = np.atleast_1d(np.asarray(t, dtype=float))
    u = np.empty((t.size,) + X.shape)
    v = np.empty_like(u)
    p = np.empty_like(u)
    for n, tn in enumerate(t):
        u[n], v[n], p[n] = taylor_green_uvp(X, Y, tn, nu, rho, advect, amplitude)
    return FlowField(
        x=x, y=y, t=t, u=u, v=v, rho=rho, mu=rho * nu, p_ref=p, U_ref=1.0,
        name=name, meta={"case": "taylor_green", "nu": nu, "advect": tuple(advect),
                         "amplitude": amplitude},
    )


def taylor_green(
    n: int = 64,
    x_range: tuple[float, float] = TG_WINDOW,
    y_range: tuple[float, float] = TG_WINDOW,
    t_center: float = 1.0,
    dt: float = 0.05,
    nt: int = 3,
    nu: float = 0.01,
    rho: float = 1.0,
    advect: tuple[float, float] = (0.0, 0.0),
    amplitude: float = 1.0,
) -> FlowField:
    """Taylor-Green snapshots on an ``n x n`` grid over a (non-periodic) window.

    The ``nt`` snapshots are centred on ``t_center`` and spaced by ``dt``
    (``nt = 3`` gives t - dt, t, t + dt as in the spec).
    """
    x = np.linspace(*x_range, n)
    y = np.linspace(*y_range, n)
    t = t_center + dt * (np.arange(nt) - (nt - 1) / 2.0)
    name = "taylor_green" if advect == (0.0, 0.0) else "taylor_green_convected"
    return taylor_green_on_grid(x, y, t, nu, rho, advect, name=name, amplitude=amplitude)


# --------------------------------------------------------------------------
# Kovasznay
# --------------------------------------------------------------------------
def kovasznay_lambda(Re: float) -> float:
    """lambda = Re/2 - sqrt(Re^2/4 + 4 pi^2)  (negative for Re > 0)."""
    return Re / 2.0 - np.sqrt(Re**2 / 4.0 + 4.0 * np.pi**2)


def kovasznay_uvp(
    X: np.ndarray, Y: np.ndarray, Re: float = 20.0, rho: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact Kovasznay velocity and pressure (steady)."""
    lam = kovasznay_lambda(Re)
    e = np.exp(lam * X)
    u = 1.0 - e * np.cos(2.0 * np.pi * Y)
    v = lam / (2.0 * np.pi) * e * np.sin(2.0 * np.pi * Y)
    p = 0.5 * rho * (1.0 - np.exp(2.0 * lam * X))
    return u, v, p


def kovasznay_grad_p(
    X: np.ndarray, Y: np.ndarray, Re: float = 20.0, rho: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Exact pressure gradient: dp/dx = -rho lambda exp(2 lambda x), dp/dy = 0."""
    lam = kovasznay_lambda(Re)
    return -rho * lam * np.exp(2.0 * lam * X), np.zeros_like(Y)


def kovasznay_on_grid(
    x: np.ndarray, y: np.ndarray, Re: float = 20.0, rho: float = 1.0
) -> FlowField:
    """Sample the Kovasznay solution on given 1-D grids (single steady snapshot)."""
    X, Y = np.meshgrid(x, y, indexing="xy")
    u, v, p = kovasznay_uvp(X, Y, Re, rho)
    return FlowField(
        x=x, y=y, t=np.array([0.0]), u=u[None], v=v[None], rho=rho,
        mu=rho / Re, p_ref=p[None], U_ref=1.0, name=f"kovasznay_Re{Re:g}",
        meta={"case": "kovasznay", "Re": Re},
    )


def kovasznay(
    n: int = 64,
    Re: float = 20.0,
    x_range: tuple[float, float] = KOV_X,
    y_range: tuple[float, float] = KOV_Y,
    rho: float = 1.0,
) -> FlowField:
    """Kovasznay flow on an ``n x n`` grid over the spec domain."""
    return kovasznay_on_grid(np.linspace(*x_range, n), np.linspace(*y_range, n), Re, rho)
