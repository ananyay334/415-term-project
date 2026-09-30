# piv_pressure/core/iterative.py
"""Iterative solvers for the pressure-Poisson system with residual tracking.

All solvers work on the *same* assembled system as the direct solver in
core/poisson.py (same matrix A, right-hand side b, boundary conditions and
pinning), so the only difference is the solution method:

``jacobi``  p^{k+1} = p^k + w_J D^-1 (b - A p^k)              (D = diag A)
``gs``      Gauss-Seidel = SOR with omega = 1
``sor``     successive over-relaxation, red-black ordering:
            p_i <- p_i + omega (b_i - (A p)_i) / A_ii, first all "red"
            nodes (i + j even), then all "black" nodes, each half-sweep using
            the newest values.  On the 5-point stencil red nodes couple only
            to black nodes, so each half-sweep is one vectorised update and the
            ordering is "consistently ordered" (Young's SOR theory applies).
``cg``      Jacobi-preconditioned conjugate gradients on the symmetric form
            S = -W A (W = diag of the row weights that make W A symmetric).

Pseudo-time view of Jacobi (see EXPLANATIONS.md): with dx = dy = h,
D = -4/h^2 and Jacobi is exactly explicit Euler marching of
dp/dtau = lap_h p - f with dtau = w_J h^2 / 4.  Classic Jacobi (w_J = 1) sits
*at* the explicit stability limit dtau = h^2/4, where the checkerboard mode
(-1)^(i+j) has amplification factor -1.  For the all-Neumann problem that
mode is an exact eigenvector, so classic Jacobi stalls (the residual levels
off).  The default w_J = 0.95 (dtau just inside the limit) damps it.

Every iteration records the residual r = b - A p (over the unknowns that are
actually solved for) in L1_n, L2_n, Linf and scaled form r/S, and the change
norm L2_n(p^{k+1} - p^k).  Iteration stops when the scaled L2_n residual is
below ``tol`` (converged), when ``max_iter`` is reached, or when the residual
grows by a factor 1e10 or becomes non-finite (diverged).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import scipy.sparse as sp

from analysis.norms import l1_normalised, l2_normalised, linf

METHODS: tuple[str, ...] = ("jacobi", "gs", "sor", "cg")
#: Default Jacobi relaxation: dtau = 0.95 h^2/4, just inside the stability limit.
JACOBI_OMEGA = 0.95
DEFAULT_TOL = 1e-8
DEFAULT_MAX_ITER = 500_000
DIVERGENCE_FACTOR = 1e10
#: Optional per-iteration hook: callback(p, k) with the current iterate p (after
#: mean removal for the zero-mean problem) and the iteration number k >= 1.
Callback = Callable[[np.ndarray, int], None]


@dataclass
class LinearSystem:
    """The discrete Poisson system in the form the iterative solvers need.

    Attributes
    ----------
    A : system matrix (CSR); rows of prescribed unknowns are identity rows.
    b : right-hand side; for the zero-mean problem it has been made
        compatible (uniform source correction), so A p = b is solvable.
    fixed : boolean mask of prescribed unknowns (Dirichlet nodes, pinned node).
    p0 : initial guess: prescribed values on ``fixed``, zero elsewhere.
    zero_mean : True for the all-Neumann problem with the zero-mean
        constraint - the mean is removed after every iteration.
    w : row weights with W A symmetric (trapezoidal weights for the
        second-order ghost-point matrix, ones for the first-order one).
    shape : grid shape (ny, nx).
    """

    A: sp.csr_matrix
    b: np.ndarray
    fixed: np.ndarray
    p0: np.ndarray
    zero_mean: bool
    w: np.ndarray
    shape: tuple[int, int]

    @property
    def free(self) -> np.ndarray:
        """Mask of the unknowns that are iterated on."""
        return ~self.fixed

    def red_mask(self) -> np.ndarray:
        """Red nodes of the red-black ordering: (i + j) even."""
        ny, nx = self.shape
        j, i = np.divmod(np.arange(nx * ny), nx)
        return (i + j) % 2 == 0

    def residual(self, p: np.ndarray) -> np.ndarray:
        """r = b - A p on the free unknowns (prescribed rows excluded)."""
        return (self.b - self.A @ p)[self.free]


@dataclass
class IterationHistory:
    """Convergence record of one iterative solve.

    Arrays are indexed by iteration k = 0 .. iterations (k = 0 is the initial
    guess); ``change_l2n[k-1]`` is L2_n(p^k - p^{k-1}).
    """

    method: str
    omega: float
    tol: float
    scale: float
    scale_kind: str
    residual_l1n: np.ndarray = field(default_factory=lambda: np.empty(0))
    residual_l2n: np.ndarray = field(default_factory=lambda: np.empty(0))
    residual_linf: np.ndarray = field(default_factory=lambda: np.empty(0))
    change_l2n: np.ndarray = field(default_factory=lambda: np.empty(0))
    iterations: int = 0
    converged: bool = False
    diverged: bool = False
    time_s: float = 0.0

    @property
    def residual_scaled(self) -> np.ndarray:
        """Scaled residual history L2_n(r) / S."""
        return self.residual_l2n / self.scale

    def summary(self) -> dict:
        """Scalar summary (one CSV row)."""
        return {"method": self.method, "omega": self.omega, "tol": self.tol,
                "iterations": self.iterations, "converged": self.converged,
                "diverged": self.diverged, "time_s": self.time_s,
                "final_scaled_residual": float(self.residual_scaled[-1]),
                "scale": self.scale, "scale_kind": self.scale_kind}

    def to_rows(self, every: int = 1) -> list[dict]:
        """One dict per recorded iteration (every ``every``-th, plus the last)."""
        ks = list(range(0, self.iterations + 1, max(1, every)))
        if ks[-1] != self.iterations:
            ks.append(self.iterations)
        rows = []
        for k in ks:
            rows.append({
                "method": self.method, "omega": self.omega, "iteration": k,
                "residual_l1n": float(self.residual_l1n[k]),
                "residual_l2n": float(self.residual_l2n[k]),
                "residual_linf": float(self.residual_linf[k]),
                "residual_scaled": float(self.residual_l2n[k] / self.scale),
                "change_l2n": float(self.change_l2n[k - 1]) if k > 0 else float("nan"),
            })
        return rows


class _Recorder:
    """Collects residual/change norms and applies the stopping rule."""

    def __init__(self, method: str, omega: float, tol: float, scale: float, scale_kind: str,
                 max_iter: int) -> None:
        self.h = IterationHistory(method, omega, tol, scale, scale_kind)
        self.l1n, self.l2n, self.linf, self.change = [], [], [], []
        self.max_iter = max_iter
        self.t0 = time.perf_counter()

    def record(self, r: np.ndarray, change: float | None) -> bool:
        """Store the norms of residual ``r``; return True if iteration must stop."""
        l2n = l2_normalised(r) if r.size else 0.0
        self.l1n.append(l1_normalised(r) if r.size else 0.0)
        self.l2n.append(l2n)
        self.linf.append(linf(r) if r.size else 0.0)
        if change is not None:
            self.change.append(change)
        k = len(self.l2n) - 1
        if not np.isfinite(l2n) or l2n > DIVERGENCE_FACTOR * max(self.l2n[0], 1e-300):
            self.h.diverged = True
            return True
        if l2n / self.h.scale < self.h.tol:
            self.h.converged = True
            return True
        return k >= self.max_iter

    def finish(self) -> IterationHistory:
        """Freeze the arrays and the timing."""
        self.h.residual_l1n = np.array(self.l1n)
        self.h.residual_l2n = np.array(self.l2n)
        self.h.residual_linf = np.array(self.linf)
        self.h.change_l2n = np.array(self.change)
        self.h.iterations = len(self.l2n) - 1
        self.h.time_s = time.perf_counter() - self.t0
        return self.h


def _project(p: np.ndarray, system: LinearSystem) -> np.ndarray:
    """Remove the mean for the zero-mean (all-Neumann) problem."""
    return p - p.mean() if system.zero_mean else p


def _scale_of(system: LinearSystem, scale: float | None) -> tuple[float, str]:
    """Residual scale S: L2_n of the free right-hand side unless given."""
    if scale is not None:
        return float(scale), "given"
    S = l2_normalised(system.b[system.free])
    return (S if S > 0 else 1.0), "rhs"


# --------------------------------------------------------------------------
# SOR optimum
# --------------------------------------------------------------------------
def jacobi_spectral_radius(nx: int, ny: int, dx: float, dy: float, neumann: bool = True) -> float:
    """Estimate of the Jacobi spectral radius rho_J of the 5-point Laplacian.

    Eigenvalues of D^-1 A are mu = [4/dx^2 sin^2(a) + 4/dy^2 sin^2(b)] /
    (2/dx^2 + 2/dy^2).  rho_J = 1 - (smallest relevant mu):
    * Dirichlet walls (the model problem): lowest mode sin(pi x/Lx) sin(pi y/Ly),
      a = pi/(2(nx-1)), b = pi/(2(ny-1)); for a square with dx = dy this gives
      rho_J = cos(pi h / L).
    * All-Neumann walls: the constant mode is the null space; the lowest
      non-constant mode is cos(pi x/L) (one direction only), i.e. half the
      Dirichlet eigenvalue, so rho_J is closer to 1.
    """
    ax, ay = 1.0 / dx**2, 1.0 / dy**2
    diag = 2.0 * ax + 2.0 * ay
    sx = 4.0 * ax * np.sin(np.pi / (2.0 * (nx - 1))) ** 2
    sy = 4.0 * ay * np.sin(np.pi / (2.0 * (ny - 1))) ** 2
    mu = min(sx, sy) / diag if neumann else (sx + sy) / diag
    return float(1.0 - mu)


def omega_opt(nx: int, ny: int, dx: float, dy: float, neumann: bool = True) -> float:
    """Young's optimum SOR factor omega = 2 / (1 + sqrt(1 - rho_J^2)).

    Assumptions: 5-point Laplacian, consistently ordered unknowns (true for
    red-black ordering), Jacobi eigenvalues real with |mu_J| <= rho_J < 1
    apart from the null-space pair.  For the Dirichlet model problem on a
    square this is the textbook omega = 2 / (1 + sin(pi h / L)).
    """
    rho = jacobi_spectral_radius(nx, ny, dx, dy, neumann)
    return float(2.0 / (1.0 + np.sqrt(1.0 - rho**2)))


def omega_opt_dirichlet_formula(h: float, L: float) -> float:
    """The lecture/spec formula omega_opt = 2 / (1 + sin(pi h / L))."""
    return float(2.0 / (1.0 + np.sin(np.pi * h / L)))


# --------------------------------------------------------------------------
# Solvers
# --------------------------------------------------------------------------
def jacobi(system: LinearSystem, omega: float = JACOBI_OMEGA, tol: float = DEFAULT_TOL,
           max_iter: int = DEFAULT_MAX_ITER, scale: float | None = None,
           callback: Callback | None = None) -> tuple[np.ndarray, IterationHistory]:
    """(Weighted) Jacobi: p <- p + omega D^-1 (b - A p) for all free nodes at once."""
    S, kind = _scale_of(system, scale)
    rec = _Recorder("jacobi", omega, tol, S, kind, max_iter)
    D = system.A.diagonal()
    free = system.free
    p = _project(system.p0.copy(), system)
    r_full = system.b - system.A @ p
    stop = rec.record(r_full[free], None)
    while not stop:
        p_new = p.copy()
        p_new[free] += omega * r_full[free] / D[free]
        p_new = _project(p_new, system)
        change = l2_normalised(p_new - p)
        p = p_new
        r_full = system.b - system.A @ p
        stop = rec.record(r_full[free], change)
        if callback is not None:
            callback(p, len(rec.l2n) - 1)
    return p, rec.finish()


def sor(system: LinearSystem, omega: float, tol: float = DEFAULT_TOL,
        max_iter: int = DEFAULT_MAX_ITER, scale: float | None = None,
        name: str = "sor", callback: Callback | None = None) -> tuple[np.ndarray, IterationHistory]:
    """Red-black SOR (``omega = 1`` gives red-black Gauss-Seidel)."""
    S, kind = _scale_of(system, scale)
    rec = _Recorder(name, omega, tol, S, kind, max_iter)
    A, b = system.A, system.b
    D = A.diagonal()
    red = system.red_mask()
    half = [np.nonzero(red & system.free)[0], np.nonzero(~red & system.free)[0]]
    rows = [A[idx] for idx in half]                    # A restricted to red / black rows
    p = _project(system.p0.copy(), system)
    stop = rec.record(system.residual(p), None)
    while not stop:
        p_old = p.copy()
        for idx, A_half in zip(half, rows):            # red half-sweep, then black
            p[idx] += omega * (b[idx] - A_half @ p) / D[idx]
        p = _project(p, system)
        stop = rec.record(system.residual(p), l2_normalised(p - p_old))
        if callback is not None:
            callback(p, len(rec.l2n) - 1)
    return p, rec.finish()


def gauss_seidel(system: LinearSystem, tol: float = DEFAULT_TOL, max_iter: int = DEFAULT_MAX_ITER,
                 scale: float | None = None,
                 callback: Callback | None = None) -> tuple[np.ndarray, IterationHistory]:
    """Red-black Gauss-Seidel (SOR with omega = 1)."""
    return sor(system, 1.0, tol, max_iter, scale, name="gs", callback=callback)


def conjugate_gradient(system: LinearSystem, tol: float = DEFAULT_TOL, max_iter: int = DEFAULT_MAX_ITER,
                       scale: float | None = None,
                       callback: Callback | None = None) -> tuple[np.ndarray, IterationHistory]:
    """Jacobi-preconditioned CG on the symmetric form S = -W A (free unknowns).

    Prescribed unknowns are eliminated (their columns move to the right-hand
    side).  The recorded residual is that of the *original* system A p = b,
    so the history is directly comparable with the other solvers.
    """
    S_scale, kind = _scale_of(system, scale)
    rec = _Recorder("cg", float("nan"), tol, S_scale, kind, max_iter)
    free, fixed = system.free, system.fixed
    Sym = (-(sp.diags(system.w) @ system.A)).tocsr()
    S_ff = Sym[free][:, free]
    rhs = -(system.w * system.b)[free] - Sym[free][:, fixed] @ system.p0[fixed]
    M_inv = 1.0 / S_ff.diagonal()                      # Jacobi preconditioner

    p = system.p0.copy()
    x = p[free].copy()
    r = rhs - S_ff @ x
    z = M_inv * r
    d = z.copy()
    rz = r @ z
    p_report = _project(p, system)
    stop = rec.record(system.residual(p), None)
    while not stop:
        Sd = S_ff @ d
        dSd = d @ Sd
        if dSd <= 0.0:                                 # exact solution reached (d = 0)
            break
        alpha = rz / dSd
        x += alpha * d
        r -= alpha * Sd
        z = M_inv * r
        rz_new = r @ z
        d = z + (rz_new / rz) * d
        rz = rz_new
        p[free] = x
        p_new = _project(p, system)
        change = l2_normalised(p_new - p_report)
        p_report = p_new
        stop = rec.record(system.residual(p), change)
        if callback is not None:
            callback(p_report, len(rec.l2n) - 1)
    return p_report, rec.finish()


def solve(system: LinearSystem, method: str, tol: float = DEFAULT_TOL,
          max_iter: int = DEFAULT_MAX_ITER, omega: float | None = None,
          scale: float | None = None, neumann: bool = True,
          spacing: tuple[float, float] = (1.0, 1.0),
          callback: Callback | None = None) -> tuple[np.ndarray, IterationHistory]:
    """Dispatch to one of :data:`METHODS`.

    ``omega=None`` selects the default: 0.95 for Jacobi, Young's optimum
    (:func:`omega_opt` for the grid, walls ``neumann`` or Dirichlet) for SOR.
    """
    if method == "jacobi":
        return jacobi(system, JACOBI_OMEGA if omega is None else omega, tol, max_iter, scale, callback)
    if method == "gs":
        return gauss_seidel(system, tol, max_iter, scale, callback)
    if method == "sor":
        ny, nx = system.shape
        w_sor = omega_opt(nx, ny, spacing[0], spacing[1], neumann) if omega is None else omega
        return sor(system, w_sor, tol, max_iter, scale, callback=callback)
    if method == "cg":
        return conjugate_gradient(system, tol, max_iter, scale, callback)
    raise ValueError(f"unknown iterative method {method!r}; use one of {METHODS}")
