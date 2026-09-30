# piv_pressure/tests/test_poisson_exact.py
"""Poisson solver: manufactured solution p = sin x cos y with exact Neumann data."""
import numpy as np
import pytest
import scipy.sparse as sp

from analysis.metrics import convergence_slope, eps_l2
from core.poisson import (assemble_neumann_laplacian, compatibility_residual, neumann_rhs,
                          solve_poisson)

X_RANGE, Y_RANGE = (0.3, 2.5), (0.1, 2.0)   # non-symmetric window, dx != dy


def manufactured(nx, ny):
    """p = sin x cos y, G = grad p, f = lap p = -2 p."""
    x, y = np.linspace(*X_RANGE, nx), np.linspace(*Y_RANGE, ny)
    X, Y = np.meshgrid(x, y)
    p = np.sin(X) * np.cos(Y)
    Gx, Gy = np.cos(X) * np.cos(Y), -np.sin(X) * np.sin(Y)
    return x, y, p, -2.0 * p, Gx, Gy


@pytest.mark.parametrize("pin", ["mean", "node"])
@pytest.mark.parametrize("solver", ["direct", "cg"])
def test_manufactured_second_order(pin, solver):
    hs, errs = [], []
    for n in [17, 33, 65, 129]:
        x, y, p, f, Gx, Gy = manufactured(n, n + 4)
        ph, info = solve_poisson(f, Gx, Gy, x[1] - x[0], y[1] - y[0], pin=pin, solver=solver)
        hs.append(x[1] - x[0])
        errs.append(eps_l2(ph, p))
        if pin == "node" and solver == "direct":
            assert ph.ravel()[0] == pytest.approx(0.0, abs=1e-12)   # the replaced row
    slope = convergence_slope(hs, errs)
    assert 1.8 < slope < 2.4, f"slope {slope:.2f}"
    assert errs[-1] < 1e-3


def test_direct_and_cg_agree():
    x, y, p, f, Gx, Gy = manufactured(40, 33)
    a, _ = solve_poisson(f, Gx, Gy, x[1] - x[0], y[1] - y[0], solver="direct")
    b, _ = solve_poisson(f, Gx, Gy, x[1] - x[0], y[1] - y[0], solver="cg", tol=1e-12)
    assert np.max(np.abs((a - a.mean()) - (b - b.mean()))) < 1e-8


def test_matrix_properties_and_compatibility():
    nx, ny, dx, dy = 9, 7, 0.2, 0.3
    A, w = assemble_neumann_laplacian(nx, ny, dx, dy)
    assert np.allclose(A @ np.ones(nx * ny), 0.0)            # constants in the null space
    assert np.allclose(w @ A, 0.0, atol=1e-12)                # left null vector = trapezoid weights
    WA = sp.diags(w) @ A
    assert abs(WA - WA.T).max() < 1e-12                       # scaled matrix is symmetric
    # trapezoidal compatibility residual == dx*dy*sum(w b)
    x, y, p, f, Gx, Gy = manufactured(nx, ny)
    dx, dy = x[1] - x[0], y[1] - y[0]
    R, _ = compatibility_residual(f, Gx, Gy, x, y)
    b = neumann_rhs(f, Gx, Gy, dx, dy).ravel()
    assert R == pytest.approx(dx * dy * np.sum(w * b), rel=1e-10, abs=1e-12)


def test_compatibility_residual_converges():
    res = []
    for n in [17, 33, 65]:
        x, y, p, f, Gx, Gy = manufactured(n, n)
        res.append(abs(compatibility_residual(f, Gx, Gy, x, y)[0]))
    assert res[2] < res[1] < res[0]
    assert np.log(res[0] / res[2]) / np.log(4.0) > 1.8       # O(h^2) (trapezoidal rule)


def test_dirichlet_edge_exact_values():
    x, y, p, f, Gx, Gy = manufactured(65, 65)
    ph, info = solve_poisson(f, Gx, Gy, x[1] - x[0], y[1] - y[0],
                             dirichlet_edges=("top",), dirichlet_values=p)
    assert np.allclose(ph[-1, :], p[-1, :])                   # prescribed values kept
    assert eps_l2(ph, p) < 1e-3
    ph2, _ = solve_poisson(f, Gx, Gy, x[1] - x[0], y[1] - y[0], dirichlet_edges=("top",),
                           dirichlet_values=p, solver="cg", tol=1e-12)
    assert np.max(np.abs(ph - ph2)) < 1e-8


def test_linear_pressure_is_exact():
    """5-point Laplacian + ghost-point Neumann reproduce a linear field exactly."""
    x, y = np.linspace(0, 1, 11), np.linspace(0, 2, 15)
    X, Y = np.meshgrid(x, y)
    p = 2.0 * X - 3.0 * Y
    ph, info = solve_poisson(np.zeros_like(p), 2.0 + 0 * X, -3.0 + 0 * X, 0.1, 2 / 14)
    assert np.max(np.abs((ph - ph.mean()) - (p - p.mean()))) < 1e-10
    assert abs(info["compatibility_residual"]) < 1e-12
