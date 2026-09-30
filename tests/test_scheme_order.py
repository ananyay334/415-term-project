# piv_pressure/tests/test_scheme_order.py
"""First- vs second-order schemes: fitted grid-convergence slopes ~1 and ~2."""
import numpy as np
import pytest

from analysis.metrics import convergence_slope, eps_l2
from core.derivatives import first_derivative, second_derivative
from core.pipeline import reconstruct
from core.poisson import assemble_neumann_laplacian
from fields.analytic import taylor_green

NS = [32, 45, 64, 91, 128]


@pytest.mark.parametrize("order, lo, hi", [(1, 0.8, 1.2), (2, 1.7, 2.3)])
def test_taylor_green_reconstruction_order(order, lo, hi):
    hs, errs = [], []
    for n in NS:
        fld = taylor_green(n)
        rec = reconstruct(fld, order=order)
        hs.append(fld.dx)
        errs.append(eps_l2(rec.p, rec.p_ref))
    slope = convergence_slope(hs, errs)
    assert lo < slope < hi, f"order {order}: slope {slope:.2f}"


@pytest.mark.parametrize("order", [1, 2])
def test_derivative_orders(order):
    hs, e1, e2 = [], [], []
    for n in (16, 32, 64, 128):
        x = np.linspace(0.0, 2.0, n)
        h = x[1] - x[0]
        hs.append(h)
        e1.append(np.abs(first_derivative(np.sin(x), h, order=order) - np.cos(x)).max())
        e2.append(np.abs(second_derivative(np.sin(x), h, order=order) + np.sin(x)).max())
    assert abs(convergence_slope(hs, e1) - order) < 0.2
    assert abs(convergence_slope(hs, e2) - order) < 0.2


def test_first_order_matrix_is_symmetric_with_unit_null_vector():
    A, w = assemble_neumann_laplacian(7, 6, 0.3, 0.25, order=1)
    assert abs(A - A.T).max() < 1e-12
    assert np.allclose(w, 1.0) and np.allclose(A @ np.ones(42), 0.0)


def test_default_order_unchanged():
    fld = taylor_green(32)
    a, b = reconstruct(fld), reconstruct(fld, order=2)
    assert np.array_equal(a.p, b.p)


def test_invalid_order_rejected():
    with pytest.raises(ValueError):
        first_derivative(np.arange(5.0), 1.0, order=3)
