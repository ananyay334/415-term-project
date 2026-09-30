# piv_pressure/core/derivatives.py
"""Second-order finite differences on a uniform grid.

All operators are second-order accurate *everywhere*, including the boundary
nodes, because a first-order boundary stencil would pollute the Neumann data
of the Poisson problem and drop the global convergence rate to first order.

Stencils (h = grid spacing, f_i = f(x_i))
-----------------------------------------
First derivative
    interior      f'_i   = (f_{i+1} - f_{i-1}) / (2h)                 error -h^2 f'''/6
    left end      f'_0   = (-3 f_0 + 4 f_1 - f_2) / (2h)              error +h^2 f'''/3
    right end     f'_N   = ( 3 f_N - 4 f_{N-1} + f_{N-2}) / (2h)      error +h^2 f'''/3
Second derivative
    interior      f''_i  = (f_{i+1} - 2 f_i + f_{i-1}) / h^2          error -h^2 f''''/12
    left end      f''_0  = (2 f_0 - 5 f_1 + 4 f_2 - f_3) / h^2        error 11 h^2 f''''/12
    right end     mirror image of the left-end stencil.

First-order option (``order=1``; used only to show O(h) vs O(h^2) grid
convergence - the default is ``order=2`` everywhere)
    first derivative   f'_i  = (f_{i+1} - f_i) / h            error +h f''/2
                       (forward; backward at the last node)
    second derivative  f''_i = (f_{i+2} - 2 f_{i+1} + f_i) / h^2   error +h f'''
                       (forward; backward for the last two nodes)

Arrays are laid out as ``[..., ny, nx]`` so axis -1 is x and axis -2 is y.
"""
from __future__ import annotations

import numpy as np

ORDERS: tuple[int, ...] = (1, 2)


def _moveaxis_last(f: np.ndarray, axis: int) -> np.ndarray:
    """Return a view of ``f`` with ``axis`` moved to the last position."""
    return np.moveaxis(f, axis, -1)


def check_order(order: int) -> None:
    """Only first- and second-order schemes are implemented."""
    if order not in ORDERS:
        raise ValueError(f"order must be 1 or 2, got {order!r}")


def first_derivative(f: np.ndarray, h: float, axis: int = -1, order: int = 2) -> np.ndarray:
    """df/dx along ``axis``.

    Parameters
    ----------
    f : array with at least 3 points along ``axis``.
    h : uniform grid spacing along ``axis``.
    axis : axis to differentiate along (-1 = x, -2 = y).
    order : 2 = central interior + one-sided second-order ends (default);
        1 = forward differences (backward at the last node).
    """
    check_order(order)
    g = _moveaxis_last(np.asarray(f, dtype=float), axis)
    if g.shape[-1] < 3:
        raise ValueError("need at least 3 points for a second-order derivative")
    d = np.empty_like(g)
    if order == 1:
        d[..., :-1] = (g[..., 1:] - g[..., :-1]) / h
        d[..., -1] = (g[..., -1] - g[..., -2]) / h
        return np.moveaxis(d, -1, axis)
    d[..., 1:-1] = (g[..., 2:] - g[..., :-2]) / (2.0 * h)
    d[..., 0] = (-3.0 * g[..., 0] + 4.0 * g[..., 1] - g[..., 2]) / (2.0 * h)
    d[..., -1] = (3.0 * g[..., -1] - 4.0 * g[..., -2] + g[..., -3]) / (2.0 * h)
    return np.moveaxis(d, -1, axis)


def second_derivative(f: np.ndarray, h: float, axis: int = -1, order: int = 2) -> np.ndarray:
    """d2f/dx2 along ``axis``.

    ``order=2``: central interior + one-sided second-order ends (default);
    ``order=1``: forward one-sided stencil (backward for the last two nodes).
    """
    check_order(order)
    g = _moveaxis_last(np.asarray(f, dtype=float), axis)
    if g.shape[-1] < 4:
        raise ValueError("need at least 4 points for a second-order 2nd derivative")
    d = np.empty_like(g)
    if order == 1:
        d[..., :-2] = (g[..., 2:] - 2.0 * g[..., 1:-1] + g[..., :-2]) / h**2
        d[..., -2:] = (g[..., -2:] - 2.0 * g[..., -3:-1] + g[..., -4:-2]) / h**2
        return np.moveaxis(d, -1, axis)
    d[..., 1:-1] = (g[..., 2:] - 2.0 * g[..., 1:-1] + g[..., :-2]) / h**2
    d[..., 0] = (2.0 * g[..., 0] - 5.0 * g[..., 1] + 4.0 * g[..., 2] - g[..., 3]) / h**2
    d[..., -1] = (2.0 * g[..., -1] - 5.0 * g[..., -2] + 4.0 * g[..., -3] - g[..., -4]) / h**2
    return np.moveaxis(d, -1, axis)


def gradient(f: np.ndarray, dx: float, dy: float, order: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Return (df/dx, df/dy) of an array laid out as ``[..., ny, nx]``."""
    return first_derivative(f, dx, -1, order), first_derivative(f, dy, -2, order)


def laplacian(f: np.ndarray, dx: float, dy: float, order: int = 2) -> np.ndarray:
    """Return d2f/dx2 + d2f/dy2 of an array laid out as ``[..., ny, nx]``."""
    return second_derivative(f, dx, -1, order) + second_derivative(f, dy, -2, order)


def divergence(fx: np.ndarray, fy: np.ndarray, dx: float, dy: float, order: int = 2) -> np.ndarray:
    """Return d(fx)/dx + d(fy)/dy."""
    return first_derivative(fx, dx, -1, order) + first_derivative(fy, dy, -2, order)


def time_derivative(f: np.ndarray, dt: float, n: int) -> np.ndarray:
    """df/dt at snapshot ``n`` of a stack ``f[nt, ...]``.

    Uses the central difference (f^{n+1} - f^{n-1}) / (2 dt) for interior
    snapshots and the second-order one-sided formulas at the first/last
    snapshot, e.g. (-3 f^0 + 4 f^1 - f^2) / (2 dt).
    """
    nt = f.shape[0]
    if nt < 3:
        raise ValueError("need at least 3 snapshots for a second-order time derivative")
    if not 0 <= n < nt:
        raise IndexError(f"time index {n} out of range 0..{nt - 1}")
    if 0 < n < nt - 1:
        return (f[n + 1] - f[n - 1]) / (2.0 * dt)
    if n == 0:
        return (-3.0 * f[0] + 4.0 * f[1] - f[2]) / (2.0 * dt)
    return (3.0 * f[-1] - 4.0 * f[-2] + f[-3]) / (2.0 * dt)
