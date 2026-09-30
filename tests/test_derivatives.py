# piv_pressure/tests/test_derivatives.py
"""Second-order accuracy of the finite-difference operators (incl. boundaries)."""
import numpy as np
import pytest

from analysis.metrics import convergence_slope
from core.derivatives import (divergence, first_derivative, gradient, laplacian,
                              second_derivative, time_derivative)

NS = [16, 32, 64, 128]


def smooth(X, Y):
    """Smooth test function f = sin(1.3x + 0.4) cos(0.7y) exp(0.2x) and its derivatives."""
    a = 1.3 * X + 0.4
    f = np.sin(a) * np.cos(0.7 * Y) * np.exp(0.2 * X)
    fx = (1.3 * np.cos(a) + 0.2 * np.sin(a)) * np.cos(0.7 * Y) * np.exp(0.2 * X)
    fxx = ((0.04 - 1.69) * np.sin(a) + 2 * 0.2 * 1.3 * np.cos(a)) * np.cos(0.7 * Y) * np.exp(0.2 * X)
    fy = -0.7 * np.sin(a) * np.sin(0.7 * Y) * np.exp(0.2 * X)
    fyy = -0.49 * f
    return f, fx, fy, fxx, fyy


def grid(n):
    x = np.linspace(0.0, 2.0, n)
    y = np.linspace(-1.0, 1.5, n)
    X, Y = np.meshgrid(x, y)
    return x, y, X, Y


@pytest.mark.parametrize("which", ["fx", "fy", "fxx", "fyy", "lap"])
def test_second_order_everywhere(which):
    hs, errs, errs_edge = [], [], []
    for n in NS:
        x, y, X, Y = grid(n)
        dx, dy = x[1] - x[0], y[1] - y[0]
        f, fx, fy, fxx, fyy = smooth(X, Y)
        num, exact = {
            "fx": (first_derivative(f, dx, -1), fx),
            "fy": (first_derivative(f, dy, -2), fy),
            "fxx": (second_derivative(f, dx, -1), fxx),
            "fyy": (second_derivative(f, dy, -2), fyy),
            "lap": (laplacian(f, dx, dy), fxx + fyy),
        }[which]
        e = np.abs(num - exact)
        hs.append(dx)
        errs.append(e.max())
        errs_edge.append(max(e[:, 0].max(), e[:, -1].max(), e[0, :].max(), e[-1, :].max()))
    assert 1.85 < convergence_slope(hs, errs) < 2.3
    assert 1.85 < convergence_slope(hs, errs_edge) < 2.3   # boundary stencils are 2nd order too


def test_exact_for_polynomials():
    x = np.linspace(-1, 2, 9)
    q = 3 * x**2 - x + 2
    c = x**3 - 2 * x
    assert np.allclose(first_derivative(q, x[1] - x[0]), 6 * x - 1)
    assert np.allclose(second_derivative(c, x[1] - x[0]), 6 * x)


def test_gradient_divergence_consistency():
    x, y, X, Y = grid(40)
    dx, dy = x[1] - x[0], y[1] - y[0]
    f = np.sin(X) * np.cos(Y)
    fx, fy = gradient(f, dx, dy)
    assert np.allclose(divergence(fx, fy, dx, dy), divergence(*gradient(f, dx, dy), dx, dy))
    assert np.allclose(fx, first_derivative(f, dx, -1)) and np.allclose(fy, first_derivative(f, dy, -2))


def test_time_derivative_central_and_one_sided():
    t = np.linspace(0.0, 0.4, 5)
    stack = (2 * t**2 + 3 * t)[:, None, None] * np.ones((5, 3, 3))   # exact for quadratics
    dt = t[1] - t[0]
    for n in range(5):
        assert np.allclose(time_derivative(stack, dt, n), 4 * t[n] + 3)
