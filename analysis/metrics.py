# piv_pressure/analysis/metrics.py
"""Error metrics for reconstructed pressure fields.

Pressure from velocity is only defined up to an additive constant, so every
comparison first subtracts the spatial mean of *both* fields::

    p'     = p - mean(p),          p_ref' = p_ref - mean(p_ref)
    eps_L2 = ||p' - p_ref'||_2 / ||p_ref'||_2          (error in % = 100 eps_L2)
    eps_max = max |p' - p_ref'| / (rho U^2 / 2)

(The discrete 2-norm over all grid nodes is used; with a uniform grid this
is proportional to the continuous L2 norm and the ratio is identical.)

The error e = p' - p_ref' is also reported in the lecture's norms
(:mod:`analysis.norms`): L1_n, L2_n and Linf, in units of rho U^2.
"""
from __future__ import annotations

import numpy as np

from analysis.norms import l1_normalised, l2, l2_normalised, linf


def demean(p: np.ndarray) -> np.ndarray:
    """Subtract the spatial mean."""
    return p - np.mean(p)


def eps_l2(p: np.ndarray, p_ref: np.ndarray) -> float:
    """Relative L2 error after removing both means: L2(e) / L2(p_ref')."""
    e = demean(p) - demean(p_ref)
    return l2(e) / l2(demean(p_ref))


def eps_max(p: np.ndarray, p_ref: np.ndarray, rho: float = 1.0, U: float = 1.0) -> float:
    """Maximum error normalised by the dynamic pressure: Linf(e) / (rho U^2 / 2)."""
    return linf(demean(p) - demean(p_ref)) / (0.5 * rho * U**2)


def error_norms(p: np.ndarray, p_ref: np.ndarray) -> dict[str, float]:
    """Error e = p' - p_ref' in the lecture's norms plus the error in %.

    Returns ``err_l1n``, ``err_l2n``, ``err_linf`` (pressure units) and
    ``err_pct`` = 100 * relative L2 error (= 100 eps_L2).
    """
    e = demean(p) - demean(p_ref)
    return {"err_l1n": l1_normalised(e), "err_l2n": l2_normalised(e), "err_linf": linf(e),
            "err_pct": 100.0 * eps_l2(p, p_ref)}


def error_map(p: np.ndarray, p_ref: np.ndarray) -> np.ndarray:
    """Point-wise absolute error |p' - p_ref'|."""
    return np.abs(demean(p) - demean(p_ref))


def locate_minimum(
    p: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    region: tuple[float, float, float, float] | None = None,
    edge_band: int = 2,
) -> tuple[int, int]:
    """Indices (j, i) of the minimum of ``p`` inside ``region``.

    ``region = (xmin, xmax, ymin, ymax)`` restricts the search (e.g. to skip
    the low-pressure base region right behind a cylinder); ``edge_band``
    nodes next to the window edges are always excluded.
    """
    X, Y = np.meshgrid(x, y, indexing="xy")
    mask = np.ones(p.shape, dtype=bool)
    if edge_band > 0:
        mask[:edge_band, :] = mask[-edge_band:, :] = False
        mask[:, :edge_band] = mask[:, -edge_band:] = False
    if region is not None:
        x0, x1, y0, y1 = region
        mask &= (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
    if not mask.any():
        raise ValueError("empty search region")
    q = np.where(mask, p, np.inf)
    j, i = np.unravel_index(np.argmin(q), p.shape)
    return int(j), int(i)


def line_profile(
    p: np.ndarray, x: np.ndarray, y: np.ndarray, y_line: float
) -> tuple[np.ndarray, np.ndarray, float]:
    """Pressure along the grid row closest to ``y = y_line``.

    Returns ``(x, p(x), y_used)``.
    """
    j = int(np.argmin(np.abs(y - y_line)))
    return x, p[j, :], float(y[j])


def core_pressure(
    p: np.ndarray,
    p_ref: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    rho: float = 1.0,
    U: float = 1.0,
    region: tuple[float, float, float, float] | None = None,
    radius: float | None = None,
) -> dict[str, float]:
    """Vortex-core pressure minimum: reconstruction vs reference.

    The core is the minimum of the (demeaned) reference pressure inside
    ``region``.  The reconstructed minimum is searched within ``radius``
    (default 3 grid spacings) of that point, so a slightly displaced core is
    still found.

    Returns a dict with the core position, both minima (demeaned), their
    ratio (1 = perfect, < 1 = under-predicted suction) and the error
    normalised by rho U^2 / 2.
    """
    pr, pp = demean(p_ref), demean(p)
    j, i = locate_minimum(pr, x, y, region)
    h = float(x[1] - x[0])
    r = 3.0 * h if radius is None else radius
    X, Y = np.meshgrid(x, y, indexing="xy")
    near = np.hypot(X - x[i], Y - y[j]) <= r + 1e-12
    p_rec_min = float(pp[near].min())
    p_ref_min = float(pr[j, i])
    return {
        "x_core": float(x[i]), "y_core": float(y[j]),
        "p_core_ref": p_ref_min, "p_core_rec": p_rec_min,
        "core_ratio": p_rec_min / p_ref_min if p_ref_min != 0 else float("nan"),
        "core_error": (p_rec_min - p_ref_min) / (0.5 * rho * U**2),
    }


def convergence_slope(h: np.ndarray, err: np.ndarray) -> float:
    """Least-squares slope of log(err) vs log(h) (the observed order)."""
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    ok = (h > 0) & (err > 0) & np.isfinite(err)
    if ok.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(h[ok]), np.log(err[ok]), 1)[0])


def standard_metrics(
    p: np.ndarray, p_ref: np.ndarray, rho: float = 1.0, U: float = 1.0
) -> dict[str, float]:
    """eps_L2, eps_max and the error norms (L1_n, L2_n, Linf, %) in one dict."""
    return {"eps_l2": eps_l2(p, p_ref), "eps_max": eps_max(p, p_ref, rho, U),
            **error_norms(p, p_ref)}
