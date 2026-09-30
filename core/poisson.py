# piv_pressure/core/poisson.py
"""Pressure from the pressure-Poisson equation (the main method).

Taking the divergence of grad p = G gives (van Oudheusden 2013, Sec. 3.2)::

    lap p = div G                    in the domain
    dp/dn = n . G                    on every edge (Neumann)

Discretisation (uniform grid, spacing dx, dy; unknowns at the grid nodes)
-----------------------------------------------------------------------
* 5-point Laplacian
  (p[i+1,j] - 2p[i,j] + p[i-1,j])/dx^2 + (p[i,j+1] - 2p[i,j] + p[i,j-1])/dy^2 = f[i,j]
* Neumann edges by *ghost-point elimination*: on the left edge the centred
  difference (p[1] - p[-1])/(2 dx) = Gx gives p[-1] = p[1] - 2 dx Gx, which is
  substituted into the Laplacian, so the row becomes
  (2p[1] - 2p[0])/dx^2 + ... = f + 2 Gx/dx.
  At a corner there are two ghost nodes, one per direction; each is
  eliminated with *its own* edge condition (dp/dx = Gx and dp/dy = Gy), so
  both edge conditions are used - no averaging is needed.
* First-order option (``order=1``, only for the O(h) vs O(h^2) comparison):
  the ghost node is eliminated with the one-sided condition
  (p[0] - p[-1])/dx = Gx, i.e. p[-1] = p[0] - dx Gx, so the row becomes
  (p[1] - p[0])/dx^2 + ... = f + Gx/dx.  The matrix is then symmetric and its
  left null vector is 1 instead of the trapezoidal weights.
* The all-Neumann matrix A is singular (A @ 1 = 0).  Two fixes:
  - ``pin="node"`` replaces one row by p = 0;
  - ``pin="mean"`` solves the augmented system [[A, 1], [1^T, 0]] which
    imposes sum(p) = 0 and adds one unknown lambda that absorbs any
    incompatibility of the data as a uniform source correction.
* Compatibility (Gauss): integral(f dA) - contour(n.G ds) must be 0 for a
  solution to exist; it is computed with the trapezoidal rule and reported.
  (Multiplying edge rows by 1/2 and corner rows by 1/4 - the trapezoidal
  weights w - makes W A symmetric; w^T A = 0, so the discrete solvability
  condition is exactly sum(w b) = 0, i.e. the trapezoidal compatibility
  residual.)
* Optional Dirichlet edges, e.g. Bernoulli p = p_inf + rho (U_inf^2 - |u|^2)/2
  (valid only where the flow is steady and irrotational).
* Solvers: sparse LU (``spsolve``, default) or the iterative solvers of
  core/iterative.py (``"jacobi"``, ``"gs"``, ``"sor"``, ``"cg"``), which work on
  the same matrix, right-hand side, boundary conditions and pinning and
  record the residual history.

How errors in G propagate into p depends on the boundary-condition type and
the domain shape; Pan et al. (2016, Meas. Sci. Technol. 27, 084012) analyse
this error propagation for the pressure-Poisson solver [verify details].
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.integrate import trapezoid
from scipy.sparse.linalg import spsolve

from core import iterative
from core.derivatives import check_order, divergence, gradient

EDGES: tuple[str, ...] = ("left", "right", "bottom", "top")
SOLVERS: tuple[str, ...] = ("direct",) + iterative.METHODS


# --------------------------------------------------------------------------
# Matrix assembly
# --------------------------------------------------------------------------
def neumann_second_difference_1d(n: int, h: float, order: int = 2) -> sp.csr_matrix:
    """1-D second-difference matrix with ghost-point Neumann end rows.

    Interior rows are [1, -2, 1]/h^2.  ``order=2``: first row [-2, 2]/h^2 and
    last row [2, -2]/h^2 (ghost node = mirror, centred condition);
    ``order=1``: first row [-1, 1]/h^2 and last row [1, -1]/h^2 (ghost node
    from the one-sided condition).
    """
    check_order(order)
    main = -2.0 * np.ones(n)
    upper = np.ones(n - 1)
    lower = np.ones(n - 1)
    if order == 2:
        upper[0] = 2.0      # row 0: p[-1] = p[1] - 2 h g  ->  coefficient of p[1] doubles
        lower[-1] = 2.0     # row n-1: p[n] = p[n-2] + 2 h g -> coefficient of p[n-2] doubles
    else:
        main[0] = main[-1] = -1.0   # p[-1] = p[0] - h g  ->  -2 p[0] + p[-1] = -p[0] - h g
    return sp.diags([lower, main, upper], [-1, 0, 1], format="csr") / h**2


def trapezoid_weights_1d(n: int) -> np.ndarray:
    """Trapezoidal-rule weights (without h): 1/2 at the two ends, 1 inside."""
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    return w


def assemble_neumann_laplacian(
    nx: int, ny: int, dx: float, dy: float, order: int = 2
) -> tuple[sp.csr_matrix, np.ndarray]:
    """Return the ghost-point Neumann Laplacian A (N x N) and weights w (N,).

    Unknowns are ordered row by row: index k = j * nx + i, i.e. the C-order
    ravel of an array of shape (ny, nx).  A = I_y (x) D_xx + D_yy (x) I_x.
    ``w`` is the left null vector of A (w^T A = 0) and makes W A symmetric:
    trapezoidal weights for ``order=2``, ones for ``order=1``.
    """
    Dxx = neumann_second_difference_1d(nx, dx, order)
    Dyy = neumann_second_difference_1d(ny, dy, order)
    A = sp.kron(sp.identity(ny), Dxx) + sp.kron(Dyy, sp.identity(nx))
    if order == 2:
        w = np.outer(trapezoid_weights_1d(ny), trapezoid_weights_1d(nx)).ravel()
    else:
        w = np.ones(nx * ny)
    return A.tocsr(), w


def edge_mask(edge: str, ny: int, nx: int) -> np.ndarray:
    """Boolean (ny, nx) mask of the nodes on ``edge`` (corners included)."""
    m = np.zeros((ny, nx), dtype=bool)
    if edge == "left":
        m[:, 0] = True
    elif edge == "right":
        m[:, -1] = True
    elif edge == "bottom":
        m[0, :] = True
    elif edge == "top":
        m[-1, :] = True
    else:
        raise ValueError(f"unknown edge {edge!r}; use one of {EDGES}")
    return m


# --------------------------------------------------------------------------
# Right-hand side, compatibility, Bernoulli
# --------------------------------------------------------------------------
def source_div_G(Gx: np.ndarray, Gy: np.ndarray, dx: float, dy: float, order: int = 2) -> np.ndarray:
    """Poisson source f = div G (finite differences of the given order)."""
    return divergence(Gx, Gy, dx, dy, order)


def source_velocity(
    u: np.ndarray, v: np.ndarray, dx: float, dy: float, rho: float, order: int = 2
) -> np.ndarray:
    """2-D incompressible source f = -rho [u_x^2 + 2 u_y v_x + v_y^2].

    Valid only when div u = 0: then div(du/dt) = 0 and div(lap u) = 0, and
    div((u.grad)u) reduces to the velocity-gradient invariant above.
    """
    ux, uy = gradient(u, dx, dy, order)
    vx, vy = gradient(v, dx, dy, order)
    return -rho * (ux**2 + 2.0 * uy * vx + vy**2)


def neumann_rhs(
    f: np.ndarray, Gx: np.ndarray, Gy: np.ndarray, dx: float, dy: float, order: int = 2
) -> np.ndarray:
    """Right-hand side b of A p = b including the ghost-point Neumann terms.

    order 2 - left   : p[-1] = p[1] - 2 dx Gx   ->  b += 2 Gx / dx
              right  : p[n]  = p[n-2] + 2 dx Gx ->  b -= 2 Gx / dx
              bottom : b += 2 Gy / dy ;   top : b -= 2 Gy / dy
    order 1 - the same with factor 1 instead of 2 (p[-1] = p[0] - dx Gx).
    Corners receive both an x and a y contribution.
    """
    check_order(order)
    c = 2.0 if order == 2 else 1.0
    b = np.array(f, dtype=float, copy=True)
    b[:, 0] += c * Gx[:, 0] / dx
    b[:, -1] -= c * Gx[:, -1] / dx
    b[0, :] += c * Gy[0, :] / dy
    b[-1, :] -= c * Gy[-1, :] / dy
    return b


def compatibility_residual(
    f: np.ndarray, Gx: np.ndarray, Gy: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[float, float]:
    """Trapezoidal  R = integral(f dA) - contour(dp/dn ds)  and R / scale.

    With outward normals: contour(n.G ds) = int_right Gx dy - int_left Gx dy
    + int_top Gy dx - int_bottom Gy dx.  ``scale`` is integral(|f| dA) +
    contour(|n.G| ds), so the relative value is 0 for compatible data and at
    most 1.
    """
    area = trapezoid(trapezoid(f, x, axis=1), y)
    flux_terms = (
        trapezoid(Gx[:, -1], y), -trapezoid(Gx[:, 0], y),
        trapezoid(Gy[-1, :], x), -trapezoid(Gy[0, :], x),
    )
    flux = sum(flux_terms)
    scale = trapezoid(trapezoid(np.abs(f), x, axis=1), y) + (
        trapezoid(np.abs(Gx[:, -1]), y) + trapezoid(np.abs(Gx[:, 0]), y)
        + trapezoid(np.abs(Gy[-1, :]), x) + trapezoid(np.abs(Gy[0, :]), x)
    )
    R = float(area - flux)
    return R, float(abs(R) / scale) if scale > 0 else 0.0


def bernoulli_pressure(
    u: np.ndarray, v: np.ndarray, rho: float, U_inf: float = 1.0, p_inf: float = 0.0
) -> np.ndarray:
    """Steady irrotational Bernoulli: p = p_inf + rho (U_inf^2 - |u|^2) / 2."""
    return p_inf + 0.5 * rho * (U_inf**2 - u**2 - v**2)


# --------------------------------------------------------------------------
# System assembly shared by all solvers
# --------------------------------------------------------------------------
def poisson_system(
    f: np.ndarray,
    Gx: np.ndarray,
    Gy: np.ndarray,
    dx: float,
    dy: float,
    dirichlet_edges: tuple[str, ...] = (),
    dirichlet_values: np.ndarray | None = None,
    pin: str = "mean",
    pin_index: int = 0,
    order: int = 2,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    """Assemble A, b and the prescribed unknowns (the same for every solver).

    Returns ``(A, b, w, fixed, p_fixed, info)``: the unmodified ghost-point
    matrix and right-hand side, the weights w (w^T A = 0), the mask and
    values of prescribed unknowns (Dirichlet nodes or the pinned node) and a
    dict with the compatibility residual and the BC description.
    """
    ny, nx = f.shape
    N = nx * ny
    A, w = assemble_neumann_laplacian(nx, ny, dx, dy, order)
    b = neumann_rhs(f, Gx, Gy, dx, dy, order).ravel()
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    R, R_rel = compatibility_residual(f, Gx, Gy, x, y)
    info: dict = {
        "n_unknowns": N, "order": order,
        "compatibility_residual": R, "compatibility_relative": R_rel,
        "source_correction": 0.0,
    }
    fixed = np.zeros(N, dtype=bool)
    p_fixed = np.zeros(N)
    if dirichlet_edges:
        if dirichlet_values is None:
            raise ValueError("dirichlet_values required with dirichlet_edges")
        for edge in dirichlet_edges:
            m = edge_mask(edge, ny, nx).ravel()
            fixed |= m
            p_fixed[m] = np.asarray(dirichlet_values).ravel()[m]
        info["bc"] = "dirichlet:" + "+".join(dirichlet_edges)
        info["pin"] = "none"
    elif pin == "node":
        fixed[pin_index] = True
        info["bc"], info["pin"] = "neumann", "node"
    elif pin == "mean":
        info["bc"], info["pin"] = "neumann", "mean"
    else:
        raise ValueError(f"unknown pin option {pin!r}")
    return A, b, w, fixed, p_fixed, info


def linear_system(
    A: sp.csr_matrix, b: np.ndarray, w: np.ndarray, fixed: np.ndarray, p_fixed: np.ndarray,
    shape: tuple[int, int],
) -> tuple[iterative.LinearSystem, float]:
    """Wrap the assembled problem for the iterative solvers.

    Prescribed unknowns get identity rows (exactly as in the direct solver).
    Without them (zero-mean problem) the data are made compatible by the same
    uniform correction lambda = w^T b / w^T 1 that the augmented direct system
    applies.  Returns the system and the source correction (-lambda or 0).
    """
    if fixed.any():
        keep = sp.diags((~fixed).astype(float))
        A_mod = (keep @ A + sp.diags(fixed.astype(float))).tocsr()
        b_mod = np.where(fixed, p_fixed, b)
        return iterative.LinearSystem(A_mod, b_mod, fixed, p_fixed.copy(), False, w, shape), 0.0
    lam = float(w @ b / w.sum())
    system = iterative.LinearSystem(A.tocsr(), b - lam, fixed, np.zeros_like(b), True, w, shape)
    return system, -lam


# --------------------------------------------------------------------------
# Solvers
# --------------------------------------------------------------------------
def solve_poisson(
    f: np.ndarray,
    Gx: np.ndarray,
    Gy: np.ndarray,
    dx: float,
    dy: float,
    dirichlet_edges: tuple[str, ...] = (),
    dirichlet_values: np.ndarray | None = None,
    pin: str = "mean",
    pin_index: int = 0,
    solver: str = "direct",
    tol: float = iterative.DEFAULT_TOL,
    max_iter: int = iterative.DEFAULT_MAX_ITER,
    order: int = 2,
    omega: float | None = None,
    residual_scale: float | None = None,
) -> tuple[np.ndarray, dict]:
    """Solve lap p = f with dp/dn = n.G (and optional Dirichlet edges).

    Parameters
    ----------
    f : source term on the grid, shape (ny, nx).
    Gx, Gy : pressure gradient; only the edge values are used (Neumann data).
    dx, dy : grid spacings.
    dirichlet_edges : edges on which p is prescribed instead of dp/dn.
    dirichlet_values : (ny, nx) array holding the prescribed p on those edges.
    pin : ``"mean"`` (zero-mean constraint) or ``"node"`` (p = 0 at
        ``pin_index``); ignored when a Dirichlet edge is present.
    solver : ``"direct"`` (spsolve, default) or ``"jacobi"``, ``"gs"``,
        ``"sor"``, ``"cg"`` (core/iterative.py).
    tol, max_iter : iterative stopping rule (scaled L2_n residual < tol).
    order : 2 (default) or 1 - Neumann-BC treatment (see module docstring).
    omega : relaxation factor for Jacobi / SOR (None = default).
    residual_scale : scale S of the residual (None = L2_n of the right-hand side).

    Returns
    -------
    p : pressure, shape (ny, nx).
    info : dict with the compatibility residual, the uniform source
        correction actually applied and solver details; for iterative
        solvers ``info["history"]`` is the :class:`IterationHistory`.
    """
    ny, nx = f.shape
    N = nx * ny
    A, b, w, fixed, p_fixed, info = poisson_system(
        f, Gx, Gy, dx, dy, dirichlet_edges, dirichlet_values, pin, pin_index, order)
    info["solver"] = solver

    if solver == "direct":
        if fixed.any():
            # Replace each fixed row by the identity row "p_k = value".
            keep = sp.diags((~fixed).astype(float))
            A_mod = (keep @ A + sp.diags(fixed.astype(float))).tocsr()
            b_mod = np.where(fixed, p_fixed, b)
            p = spsolve(A_mod, b_mod)
        else:
            # Augmented system [[A, 1], [1^T, 0]] [p, lam] = [b, 0]:
            # A p + lam = b and sum(p) = 0.  lam = sum(w b)/sum(w) makes the
            # corrected data b - lam compatible (w^T A = 0).
            ones = np.ones((N, 1))
            K = sp.bmat([[A, sp.csr_matrix(ones)], [sp.csr_matrix(ones.T), None]], format="csc")
            sol = spsolve(K, np.append(b, 0.0))
            p, lam = sol[:N], sol[N]
            info["source_correction"] = float(-lam)
    elif solver in iterative.METHODS:
        system, correction = linear_system(A, b, w, fixed, p_fixed, (ny, nx))
        info["source_correction"] = correction
        all_dirichlet = set(dirichlet_edges) == set(EDGES)
        p, history = iterative.solve(system, solver, tol, max_iter, omega, residual_scale,
                                     neumann=not all_dirichlet, spacing=(dx, dy))
        info.update(history=history, iterations=history.iterations, converged=history.converged,
                    omega=history.omega)
        if solver == "cg":
            info["cg_iterations"] = history.iterations
    else:
        raise ValueError(f"unknown solver {solver!r}; use one of {SOLVERS}")
    return np.asarray(p).reshape(ny, nx), info


def poisson_pressure(
    Gx: np.ndarray,
    Gy: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    source: str = "div_G",
    u: np.ndarray | None = None,
    v: np.ndarray | None = None,
    rho: float = 1.0,
    dirichlet_edges: tuple[str, ...] = (),
    dirichlet_values: np.ndarray | None = None,
    pin: str = "mean",
    pin_index: int = 0,
    solver: str = "direct",
    tol: float = iterative.DEFAULT_TOL,
    max_iter: int = iterative.DEFAULT_MAX_ITER,
    order: int = 2,
    omega: float | None = None,
    residual_scale: float | None = None,
) -> tuple[np.ndarray, dict]:
    """Integrate G to pressure with the Poisson method.

    ``source="div_G"`` uses f = div G (general);
    ``source="velocity"`` uses the 2-D incompressible identity
    f = -rho [u_x^2 + 2 u_y v_x + v_y^2] (needs ``u`` and ``v``); the
    Neumann data always come from G.
    """
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    if source == "div_G":
        f = source_div_G(Gx, Gy, dx, dy, order)
    elif source == "velocity":
        if u is None or v is None:
            raise ValueError("source='velocity' needs u and v")
        f = source_velocity(u, v, dx, dy, rho, order)
    else:
        raise ValueError(f"unknown source {source!r}")
    p, info = solve_poisson(
        f, Gx, Gy, dx, dy, dirichlet_edges, dirichlet_values, pin, pin_index, solver, tol,
        max_iter, order, omega, residual_scale,
    )
    info["source"] = source
    return p, info
