# piv_pressure/analysis/norms.py
"""Vector norms and scaled residuals (ME-415, Lecture 11).

For a residual (or error) vector r with n components::

    L1   = sum |r_i|                 L1_n = (1/n) sum |r_i|
    L2   = sqrt( sum |r_i|^2 )       L2_n = sqrt( (1/n) sum |r_i|^2 )
    Linf = max |r_i|

* L1 weights every component equally; L2 gives more weight to the large
  components; Linf sees only the worst cell.
* L1 and L2 grow with the number of cells n (a constant residual c gives
  L1 = n c and L2 = sqrt(n) c), so a fixed tolerance on them means something
  different on every grid.  The normalised forms L1_n, L2_n are per-cell
  averages and do not depend on n.
* A *scaled* residual r_s = r / S divides by a characteristic scale S of the
  equation, so that r_s is a relative ("percentage") imbalance per cell.

Scale options for the pressure-Poisson residual (units of lap p, i.e.
pressure / length^2):

``"rhs"``       S = L2_n(b): the typical magnitude of the source term, the
                Poisson analogue of the "flux through a face" used to scale
                finite-volume residuals.  With a zero initial guess the first
                scaled residual is exactly 1.
``"physical"``  S = rho U^2 / L^2: the fundamental pressure-Laplacian scale
                of the flow (dynamic pressure over length squared).
"""
from __future__ import annotations

import numpy as np

SCALES: tuple[str, ...] = ("rhs", "physical")


def _flat(r: np.ndarray) -> np.ndarray:
    """Contiguous 1-D view of an array of any shape."""
    return np.asarray(r).ravel()


def l1(r: np.ndarray) -> float:
    """L1 norm: sum of absolute values."""
    return float(np.sum(np.abs(_flat(r))))


def l2(r: np.ndarray) -> float:
    """L2 (Euclidean) norm: sqrt of the sum of squares."""
    return float(np.linalg.norm(_flat(r)))


def linf(r: np.ndarray) -> float:
    """L-infinity norm: largest absolute value."""
    return float(np.max(np.abs(_flat(r))))


def l1_normalised(r: np.ndarray) -> float:
    """Normalised L1 norm (1/n) sum |r_i| (mean absolute value)."""
    v = _flat(r)
    return float(np.sum(np.abs(v)) / v.size)


def l2_normalised(r: np.ndarray) -> float:
    """Normalised L2 norm sqrt((1/n) sum r_i^2) (root mean square)."""
    v = _flat(r)
    return float(np.linalg.norm(v) / np.sqrt(v.size))


def scaled(r: np.ndarray, S: float) -> np.ndarray:
    """Scaled residual r_s = r / S (S > 0 is a characteristic scale)."""
    if not S > 0:
        raise ValueError(f"the scale S must be positive, got {S}")
    return np.asarray(r) / S


def residual_scale(
    b: np.ndarray | None = None,
    kind: str = "rhs",
    rho: float = 1.0,
    U: float = 1.0,
    L: float = 1.0,
) -> float:
    """Characteristic scale S of the Poisson residual.

    ``kind="rhs"``: S = L2_n(b) (needs ``b``); ``kind="physical"``:
    S = rho U^2 / L^2.  Falls back to 1 if the right-hand side is zero.
    """
    if kind == "rhs":
        if b is None:
            raise ValueError("kind='rhs' needs the right-hand side b")
        S = l2_normalised(b)
        return S if S > 0 else 1.0
    if kind == "physical":
        return rho * U**2 / L**2
    raise ValueError(f"unknown scale {kind!r}; use one of {SCALES}")


def norm_summary(r: np.ndarray) -> dict[str, float]:
    """All five norms of ``r`` in one dict (l1, l2, linf, l1n, l2n)."""
    return {"l1": l1(r), "l2": l2(r), "linf": linf(r),
            "l1n": l1_normalised(r), "l2n": l2_normalised(r)}
