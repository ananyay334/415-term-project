# piv_pressure/piv/degrade.py
"""Synthetic PIV degradation of an exact (or CFD) velocity field.

Model, applied in this order (Charonko et al. 2010; van Gent et al. 2017):

1. **Interrogation-window averaging** - a box (moving-average) filter of
   width W fine-grid points, applied to u and v of every snapshot with
   ``scipy.ndimage.uniform_filter``.
2. **Subsampling** to the PIV vector spacing h (every ``stride``-th node).
3. **Measurement noise** - independent Gaussian noise with standard
   deviation sigma * U_ref added to u and v of every snapshot (seeded).

The reference pressure is only subsampled (never filtered or perturbed):
the reconstruction is always compared with the *true* pressure.

Note on even window widths: a moving average over an even number of points
is not centred on a node (``uniform_filter`` shifts it by half a cell).
For even W the average of the two half-cell-shifted windows is used, which
is the centred trapezoidal box [1/2, 1, ..., 1, 1/2] / W - exactly a box of
physical width W * h_fine.
"""
from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.ndimage import uniform_filter

from fields.field import FlowField


def box_filter(f: np.ndarray, width: int, mode: str = "nearest") -> np.ndarray:
    """Centred 2-D moving average of ``f[ny, nx]`` over ``width`` points.

    ``width <= 1`` returns an unchanged copy.  Odd widths use a plain
    ``uniform_filter``; even widths average the four half-shifted windows
    (origins 0 and -1 on each axis), giving a centred trapezoidal box.
    """
    f = np.asarray(f, dtype=float)
    if width <= 1:
        return f.copy()
    if width % 2 == 1:
        return uniform_filter(f, size=width, mode=mode)
    acc = np.zeros_like(f)
    for oy in (0, -1):
        for ox in (0, -1):
            acc += uniform_filter(f, size=width, origin=(oy, ox), mode=mode)
    return acc / 4.0


def subsample(fld: FlowField, stride: int, crop: int = 0) -> FlowField:
    """Keep every ``stride``-th node, then drop ``crop`` nodes on each side."""
    sl = slice(None, None, stride)
    x, y = fld.x[sl], fld.y[sl]
    u, v = fld.u[:, sl, sl], fld.v[:, sl, sl]
    p = None if fld.p_ref is None else fld.p_ref[:, sl, sl]
    if crop > 0:
        c = slice(crop, -crop)
        x, y, u, v = x[c], y[c], u[:, c, c], v[:, c, c]
        p = None if p is None else p[:, c, c]
    return fld.copy_with(x=x, y=y, u=u, v=v, p_ref=p)


def add_noise(fld: FlowField, sigma: float, seed: int | None = 0) -> FlowField:
    """Add N(0, (sigma * U_ref)^2) noise independently to u and v of every snapshot."""
    if sigma <= 0.0:
        return fld.copy_with(u=fld.u.copy(), v=fld.v.copy())
    rng = np.random.default_rng(seed)
    s = sigma * fld.U_ref
    return fld.copy_with(u=fld.u + rng.normal(0.0, s, fld.u.shape),
                         v=fld.v + rng.normal(0.0, s, fld.v.shape))


def degrade(
    fld: FlowField,
    window: int = 1,
    stride: int = 1,
    sigma: float = 0.0,
    seed: int | None = 0,
    crop: int = 0,
) -> FlowField:
    """Apply the three-step PIV model: box filter -> subsample -> noise.

    Parameters
    ----------
    fld : fine-grid "true" field.
    window : box-filter width W in *fine-grid* points (<= 1: no filtering).
    stride : subsampling factor (PIV spacing h = stride * fine spacing).
    sigma : noise level as a fraction of ``fld.U_ref`` (0.01 = 1 %).
    seed : random seed for the noise (fixed for reproducibility).
    crop : PIV nodes removed from every side after subsampling (use it to
        discard nodes whose filter window reached past the data).
    """
    u = np.stack([box_filter(s, window) for s in fld.u])
    v = np.stack([box_filter(s, window) for s in fld.v])
    out = subsample(fld.copy_with(u=u, v=v), stride, crop)
    out = add_noise(out, sigma, seed)
    out.meta = dict(fld.meta, window_fine=window, stride=stride, sigma=sigma, seed=seed)
    return out


def degrade_analytic(
    truth: Callable[[np.ndarray, np.ndarray], FlowField],
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    n: int,
    window_h: float = 1.0,
    sigma: float = 0.0,
    seed: int | None = 0,
    refine: int = 4,
) -> FlowField:
    """Synthetic PIV field on an ``n x n`` grid from an analytic solution.

    The truth is sampled on a grid ``refine`` times finer than the PIV grid
    and extended by a margin so that the box filter (width ``window_h`` PIV
    spacings) never reaches past the data; after filtering and subsampling
    the margin is cropped, so the result covers exactly the requested window
    with the PIV nodes at ``linspace(x0, x1, n)``.

    Parameters
    ----------
    truth : function ``(x, y) -> FlowField`` evaluating the exact solution on
        1-D grids (e.g. ``lambda x, y: taylor_green_on_grid(x, y, t)``).
    x_range, y_range : window limits.
    n : PIV grid points per direction.
    window_h : interrogation-window width in units of the PIV spacing h
        (0 = ideal point sampling, no filtering).
    sigma : noise level (fraction of U_ref).
    seed : noise seed.
    refine : fine-grid points per PIV spacing.
    """
    hx = (x_range[1] - x_range[0]) / (n - 1)
    hy = (y_range[1] - y_range[0]) / (n - 1)
    width = int(round(window_h * refine))
    if width <= 1:                        # point sampling: no fine grid needed
        fld = truth(np.linspace(*x_range, n), np.linspace(*y_range, n))
        return degrade(fld, window=1, stride=1, sigma=sigma, seed=seed)
    margin = math.ceil(width / (2 * refine)) + 1          # in PIV spacings
    nf = (n - 1 + 2 * margin) * refine + 1
    xf = x_range[0] - margin * hx + np.arange(nf) * (hx / refine)
    yf = y_range[0] - margin * hy + np.arange(nf) * (hy / refine)
    out = degrade(truth(xf, yf), window=width, stride=refine, sigma=sigma, seed=seed, crop=margin)
    if out.nx != n or out.ny != n:
        raise RuntimeError("internal error: wrong PIV grid size after cropping")
    out.meta.update(window_h=window_h, refine=refine)
    return out
