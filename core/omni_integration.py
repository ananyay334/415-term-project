# piv_pressure/core/omni_integration.py
"""Omni-directional integration of the pressure gradient.

Idea (Liu & Katz 2006): a single integration path picks up the error of the
measured G along that path.  Averaging the line integrals of G over *many*
straight paths that enter the domain from the boundary in all directions
averages these errors out and makes the result nearly path independent.

Implementation used here (a rotating-ray variant; see EXPLANATIONS.md for how
it relates to the original "virtual boundary" algorithm [verify details
against Liu & Katz 2006]):

For every direction theta_k = 2 pi k / N (k = 0..N-1) with unit vector
e_k = (cos theta_k, sin theta_k), and every node P, trace the ray backwards
from P along -e_k until it hits the boundary at E_k(P).  Then::

    p(P) ~ p(E_k(P)) + integral_{E_k}^{P} G . e_k ds

1. **Boundary stage.**  For the M boundary nodes the entry points E are
   themselves on the boundary; p(E) is linearly interpolated between the two
   neighbouring boundary nodes along the perimeter.  Averaging over the valid
   (inward) directions gives M linear equations

       p_i = mean_k [ sum_m A_im^(k) p_m + I_i^(k) ]   ->   (I - A) p_b = c.

   (I - A) is singular (constants), and the data are not exactly compatible
   because the measured G has curl, so the augmented system
   [[I - A, 1], [1^T, 0]] [p_b, mu] = [c, 0] is solved (zero mean, uniform
   correction mu) - the direct equivalent of iterating to convergence.
2. **Interior stage.**  Each interior node is the average over all N
   directions of p_b(E_k) + integral along the ray.

Line integrals use the trapezoidal rule with n = ceil(L / (step * h))
equal sub-steps and bilinear interpolation of G (both second order).
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


def bilinear(F: np.ndarray, x: np.ndarray, y: np.ndarray, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Bilinear interpolation of F[ny, nx] (on uniform x, y) at points (X, Y).

    Points slightly outside the grid (round-off) are clamped onto it.
    """
    nx, ny = x.size, y.size
    fx = np.clip((X - x[0]) / (x[1] - x[0]), 0.0, nx - 1.0)
    fy = np.clip((Y - y[0]) / (y[1] - y[0]), 0.0, ny - 1.0)
    i = np.minimum(fx.astype(int), nx - 2)
    j = np.minimum(fy.astype(int), ny - 2)
    tx = fx - i
    ty = fy - j
    return ((1 - tx) * (1 - ty) * F[j, i] + tx * (1 - ty) * F[j, i + 1]
            + (1 - tx) * ty * F[j + 1, i] + tx * ty * F[j + 1, i + 1])


def boundary_nodes(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Boundary node indices in counter-clockwise order from the SW corner.

    Returns ``(j, i, s)``: row/column indices and the arc-length coordinate
    s along the perimeter (s = 0 at the SW corner).
    """
    nx, ny = x.size, y.size
    i = np.concatenate([
        np.arange(nx),                 # bottom edge, west -> east
        np.full(ny - 1, nx - 1),       # right edge, south -> north
        np.arange(nx - 2, -1, -1),     # top edge, east -> west
        np.zeros(ny - 2, int),         # left edge, north -> south
    ])
    j = np.concatenate([
        np.zeros(nx, int),
        np.arange(1, ny),
        np.full(nx - 1, ny - 1),
        np.arange(ny - 2, 0, -1),
    ])
    s = arc_length(x[i], y[j], x, y)
    return j, i, s


def arc_length(xe: np.ndarray, ye: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Perimeter coordinate s of boundary points (xe, ye), counter-clockwise.

    bottom: s = x - x0;  right: Lx + (y - y0);  top: Lx + Ly + (x1 - x);
    left: 2 Lx + Ly + (y1 - y).  The edge is chosen by the nearest side.
    """
    x0, x1, y0, y1 = x[0], x[-1], y[0], y[-1]
    Lx, Ly = x1 - x0, y1 - y0
    d = np.stack([np.abs(ye - y0), np.abs(xe - x1), np.abs(ye - y1), np.abs(xe - x0)])
    side = np.argmin(d, axis=0)  # 0 bottom, 1 right, 2 top, 3 left
    s = np.select(
        [side == 0, side == 1, side == 2, side == 3],
        [xe - x0, Lx + (ye - y0), Lx + Ly + (x1 - xe), 2 * Lx + Ly + (y1 - ye)],
    )
    return np.mod(s, 2.0 * (Lx + Ly))


def _perimeter_weights(
    s_query: np.ndarray, s_nodes: np.ndarray, perimeter: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Linear interpolation weights along the (periodic) perimeter.

    Returns node indices k0, k1 and weights w0, w1 with
    value(s) = w0 * value[k0] + w1 * value[k1].
    """
    M = s_nodes.size
    s = np.mod(s_query, perimeter)
    k0 = np.clip(np.searchsorted(s_nodes, s, side="right") - 1, 0, M - 1)
    k1 = (k0 + 1) % M
    s1 = np.where(k1 == 0, perimeter, s_nodes[k1])
    t = (s - s_nodes[k0]) / (s1 - s_nodes[k0])
    return k0, k1, 1.0 - t, t


def _exit_distance(
    Px: np.ndarray, Py: np.ndarray, dx_: float, dy_: float, x: np.ndarray, y: np.ndarray
) -> np.ndarray:
    """Distance from P along direction (dx_, dy_) to the rectangle boundary."""
    with np.errstate(divide="ignore", invalid="ignore"):
        if dx_ > 0:
            Lx = (x[-1] - Px) / dx_
        elif dx_ < 0:
            Lx = (x[0] - Px) / dx_
        else:
            Lx = np.full_like(Px, np.inf)
        if dy_ > 0:
            Ly = (y[-1] - Py) / dy_
        elif dy_ < 0:
            Ly = (y[0] - Py) / dy_
        else:
            Ly = np.full_like(Py, np.inf)
    return np.maximum(np.minimum(Lx, Ly), 0.0)


def _ray_integrals(
    Ge: np.ndarray, x: np.ndarray, y: np.ndarray,
    Px: np.ndarray, Py: np.ndarray, ex: float, ey: float, L: np.ndarray,
    ds: float, max_samples: int = 2_000_000,
) -> np.ndarray:
    """Trapezoidal integral of Ge = G.e along the segments from P - L e to P.

    ``Ge`` is the gradient already projected on the ray direction e (one
    grid field per direction), so only one interpolation per sample is needed.
    """
    n = np.maximum(1, np.ceil(L / ds).astype(int))   # sub-steps per ray
    out = np.empty(L.size)
    order = np.argsort(n)                            # group rays of similar length
    rows = max(1, max_samples // (int(n.max()) + 1))  # rays per chunk (bounded memory)
    for start in range(0, L.size, rows):
        idx = order[start:start + rows]
        n_max = n[idx].max()
        k = np.arange(n_max + 1)[None, :]
        nk = n[idx][:, None]
        frac = np.minimum(k / nk, 1.0)               # 0 at P, 1 at the entry point E
        Qx = Px[idx, None] - L[idx, None] * frac * ex
        Qy = Py[idx, None] - L[idx, None] * frac * ey
        g = bilinear(Ge, x, y, Qx, Qy)
        wts = np.where(k < nk, 1.0, 0.0)             # trapezoid weights 1/2,1,...,1,1/2
        wts[:, 0] = 0.5
        wts[np.arange(idx.size), n[idx]] = 0.5
        out[idx] = np.sum(wts * g, axis=1) * (L[idx] / n[idx])
    return out


def omni_integration(
    Gx: np.ndarray,
    Gy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    n_directions: int = 32,
    step: float = 1.0,
) -> tuple[np.ndarray, dict]:
    """Omni-directional (rotating-ray) integration of G.

    Parameters
    ----------
    Gx, Gy : pressure gradient, shape (ny, nx).
    x, y : uniform 1-D grids.
    n_directions : number of ray directions N (uniform in [0, 2 pi)).
    step : integration sub-step as a fraction of the grid spacing.

    Returns
    -------
    p : zero-mean pressure, shape (ny, nx).
    info : dict with the direction count, the boundary-stage correction mu,
        the RMS spread of the single-direction estimates (a path-dependence
        indicator) and the wall-clock time.
    """
    t0 = time.perf_counter()
    ny, nx = Gx.shape
    h = min(x[1] - x[0], y[1] - y[0])
    ds = step * h
    thetas = 2.0 * np.pi * np.arange(n_directions) / n_directions
    jb, ib, s_nodes = boundary_nodes(x, y)
    M = s_nodes.size
    perimeter = 2.0 * ((x[-1] - x[0]) + (y[-1] - y[0]))
    Pbx, Pby = x[ib], y[jb]

    # ---------------- stage 1: boundary nodes -----------------------------
    rows, cols, vals = [], [], []
    c = np.zeros(M)
    count = np.zeros(M)
    for th in thetas:
        ex, ey = np.cos(th), np.sin(th)
        ex = 0.0 if abs(ex) < 1e-12 else ex    # exact axis directions
        ey = 0.0 if abs(ey) < 1e-12 else ey
        L = _exit_distance(Pbx, Pby, -ex, -ey, x, y)
        ok = L > 1e-9 * h                      # ray points into the domain
        if not np.any(ok):
            continue
        idx = np.nonzero(ok)[0]
        Ex, Ey = Pbx[idx] - L[idx] * ex, Pby[idx] - L[idx] * ey
        k0, k1, w0, w1 = _perimeter_weights(arc_length(Ex, Ey, x, y), s_nodes, perimeter)
        I = _ray_integrals(Gx * ex + Gy * ey, x, y, Pbx[idx], Pby[idx], ex, ey, L[idx], ds)
        rows += [idx, idx]
        cols += [k0, k1]
        vals += [w0, w1]
        c[idx] += I
        count[idx] += 1
    if np.any(count == 0):
        raise RuntimeError("a boundary node received no rays; increase n_directions")
    rows_a, cols_a, vals_a = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    A = sp.csr_matrix((vals_a / count[rows_a], (rows_a, cols_a)), shape=(M, M))
    c /= count
    ones = sp.csr_matrix(np.ones((M, 1)))
    K = sp.bmat([[sp.identity(M) - A, ones], [ones.T, None]], format="csc")
    sol = spsolve(K, np.append(c, 0.0))
    p_b, mu = sol[:M], sol[M]

    # ---------------- stage 2: interior nodes -----------------------------
    p = np.zeros((ny, nx))
    p[jb, ib] = p_b
    J, I_ = np.meshgrid(np.arange(1, ny - 1), np.arange(1, nx - 1), indexing="ij")
    Px, Py = x[I_.ravel()], y[J.ravel()]
    acc = np.zeros(Px.size)
    acc2 = np.zeros(Px.size)
    for th in thetas:
        ex, ey = np.cos(th), np.sin(th)
        ex = 0.0 if abs(ex) < 1e-12 else ex
        ey = 0.0 if abs(ey) < 1e-12 else ey
        L = _exit_distance(Px, Py, -ex, -ey, x, y)
        Ex, Ey = Px - L * ex, Py - L * ey
        k0, k1, w0, w1 = _perimeter_weights(arc_length(Ex, Ey, x, y), s_nodes, perimeter)
        I = _ray_integrals(Gx * ex + Gy * ey, x, y, Px, Py, ex, ey, L, ds)
        est = w0 * p_b[k0] + w1 * p_b[k1] + I
        acc += est
        acc2 += est**2
    mean_est = acc / n_directions
    p[J.ravel(), I_.ravel()] = mean_est
    spread = np.sqrt(np.maximum(acc2 / n_directions - mean_est**2, 0.0))
    p -= p.mean()
    info = {
        "n_directions": n_directions,
        "boundary_correction": float(mu),
        "direction_spread_rms": float(np.sqrt(np.mean(spread**2))),
        "time_s": time.perf_counter() - t0,
    }
    return p, info
