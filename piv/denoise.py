# piv_pressure/piv/denoise.py
"""Spatial denoising filters applied to PIV velocity before differentiation.

* **Gaussian filter** - convolution with a Gaussian of standard deviation
  ``size`` grid points (``scipy.ndimage.gaussian_filter``).  Attenuates a
  Fourier mode of wavenumber k by exp(-(k s h)^2 / 2).
* **2-D Savitzky-Golay filter** - at every node, a least-squares fit of a
  2-D polynomial of total degree ``order`` to the ``size x size``
  neighbourhood, evaluated at the node.  It preserves polynomials up to that
  degree exactly, so it removes noise with less bias on smooth peaks than a
  box or Gaussian filter.  Near the edges the window is shifted to stay
  inside the data and the fit is evaluated off-centre (the 2-D analogue of
  ``scipy.signal.savgol_filter(mode="interp")``).
"""
from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import gaussian_filter

from fields.field import FlowField

DENOISERS: tuple[str, ...] = ("none", "gaussian", "savgol")


def gaussian_denoise(f: np.ndarray, size: float = 1.0) -> np.ndarray:
    """Gaussian smoothing of a 2-D array; ``size`` is the std in grid points."""
    return gaussian_filter(np.asarray(f, dtype=float), sigma=size, mode="nearest")


def savgol_exponents(order: int) -> list[tuple[int, int]]:
    """Monomials x^a y^b with a + b <= order, lowest degree first."""
    return [(a, d - a) for d in range(order + 1) for a in range(d, -1, -1)]


def savgol_coefficients(size: int, order: int) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Least-squares matrix C mapping window values to polynomial coefficients.

    For window offsets (xo, yo) in [-m, m]^2 (m = size // 2) the Vandermonde
    matrix V has one row per window point and one column per monomial; the
    fitted coefficients are c = pinv(V) @ values, so C = pinv(V).
    """
    if size % 2 == 0 or size < 3:
        raise ValueError("Savitzky-Golay window size must be odd and >= 3")
    m = size // 2
    off = np.arange(-m, m + 1)
    Yo, Xo = np.meshgrid(off, off, indexing="ij")        # row = y offset, col = x offset
    exps = savgol_exponents(order)
    if len(exps) > size * size:
        raise ValueError("polynomial order too high for the window")
    V = np.stack([Xo.ravel() ** a * Yo.ravel() ** b for a, b in exps], axis=1)
    return np.linalg.pinv(V), exps


def savitzky_golay_2d(f: np.ndarray, size: int = 5, order: int = 2) -> np.ndarray:
    """2-D Savitzky-Golay smoothing of ``f[ny, nx]`` with shifted edge windows."""
    f = np.asarray(f, dtype=float)
    ny, nx = f.shape
    m = size // 2
    if ny < size or nx < size:
        raise ValueError("array smaller than the Savitzky-Golay window")
    C, exps = savgol_coefficients(size, order)
    # Polynomial coefficients for every window that fits inside the array.
    windows = sliding_window_view(f, (size, size))              # (ny-2m, nx-2m, size, size)
    coef = windows.reshape(ny - 2 * m, nx - 2 * m, -1) @ C.T    # (ny-2m, nx-2m, n_terms)
    # Each node uses the window centred on it if possible, otherwise the
    # nearest window inside the array, evaluated at its offset (a_x, a_y).
    jc = np.clip(np.arange(ny), m, ny - 1 - m)
    ic = np.clip(np.arange(nx), m, nx - 1 - m)
    ay = (np.arange(ny) - jc).astype(float)[:, None]            # y offset of node in window
    ax = (np.arange(nx) - ic).astype(float)[None, :]
    c = coef[(jc - m)[:, None], (ic - m)[None, :]]              # (ny, nx, n_terms)
    basis = np.stack([ax ** a * ay ** b * np.ones((ny, nx)) for a, b in exps], axis=-1)
    return np.sum(c * basis, axis=-1)


def denoise_field(fld: FlowField, method: str = "none", size: float = 1.0, order: int = 2) -> FlowField:
    """Apply a denoiser to u and v of every snapshot.

    ``method``: ``"none"``, ``"gaussian"`` (``size`` = std in grid points) or
    ``"savgol"`` (``size`` = odd window width, ``order`` = polynomial degree).
    """
    if method == "none":
        return fld
    if method == "gaussian":
        filt = lambda a: gaussian_denoise(a, size)                   # noqa: E731
    elif method == "savgol":
        filt = lambda a: savitzky_golay_2d(a, int(size), order)      # noqa: E731
    else:
        raise ValueError(f"unknown denoiser {method!r}; use one of {DENOISERS}")
    u = np.stack([filt(s) for s in fld.u])
    v = np.stack([filt(s) for s in fld.v])
    return fld.copy_with(u=u, v=v, meta=dict(fld.meta, denoise=method, denoise_size=size))
