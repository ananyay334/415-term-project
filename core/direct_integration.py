# piv_pressure/core/direct_integration.py
"""Direct spatial integration of the pressure gradient (baseline method).

Pressure differences are line integrals of G = grad p::

    p(B) - p(A) = integral_A^B G . dl

Starting from p = 0 at a reference corner, the integral is marched with the
cumulative trapezoidal rule (second-order):

* path 1 - along the first row in x, then up every column in y
  p(x_i, y_0) = int_{x_0}^{x_i} Gx(x, y_0) dx
  p(x_i, y_j) = p(x_i, y_0) + int_{y_0}^{y_j} Gy(x_i, y) dy
* path 2 - along the first column in y, then along every row in x.

For an exact gradient field both paths agree; for measured (noisy,
non-conservative) G they differ, and the difference measures the path
dependence (van Oudheusden 2013, Sec. 3.1; Charonko et al. 2010).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

CORNERS: tuple[str, ...] = ("SW", "SE", "NW", "NE")


def _integrate_from_first_corner(
    Gx: np.ndarray, Gy: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Path-1 and path-2 pressures with p = 0 at node [0, 0] (x[0], y[0])."""
    # Path 1: x along row j = 0, then y along each column.
    p_row0 = cumulative_trapezoid(Gx[0, :], x, initial=0.0)                # (nx,)
    p1 = p_row0[None, :] + cumulative_trapezoid(Gy, y, axis=0, initial=0.0)
    # Path 2: y along column i = 0, then x along each row.
    p_col0 = cumulative_trapezoid(Gy[:, 0], y, initial=0.0)                # (ny,)
    p2 = p_col0[:, None] + cumulative_trapezoid(Gx, x, axis=1, initial=0.0)
    return p1, p2


def direct_integration(
    Gx: np.ndarray,
    Gy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    corner: str = "SW",
) -> dict[str, np.ndarray]:
    """Integrate G from a reference corner along two orthogonal paths.

    Parameters
    ----------
    Gx, Gy : pressure-gradient components, shape (ny, nx).
    x, y : 1-D grid coordinates.
    corner : reference corner where p = 0: ``"SW"`` (x[0], y[0]), ``"SE"``,
        ``"NW"`` or ``"NE"``.

    Returns
    -------
    dict with ``"path1"`` (x first, then y), ``"path2"`` (y first, then x)
    and ``"average"`` = (path1 + path2)/2.
    """
    if corner not in CORNERS:
        raise ValueError(f"corner must be one of {CORNERS}")
    # Reverse the axes so that the chosen corner becomes node [0, 0].  The
    # coordinates are reversed too, so dx < 0 and cumulative_trapezoid
    # integrates in the correct (negative) direction automatically.
    flip_x = corner in ("SE", "NE")
    flip_y = corner in ("NW", "NE")
    sx = slice(None, None, -1) if flip_x else slice(None)
    sy = slice(None, None, -1) if flip_y else slice(None)
    p1, p2 = _integrate_from_first_corner(Gx[sy, sx], Gy[sy, sx], x[sx], y[sy])
    p1, p2 = p1[sy, sx], p2[sy, sx]  # undo the flips (reversal is its own inverse)
    return {"path1": p1, "path2": p2, "average": 0.5 * (p1 + p2)}
