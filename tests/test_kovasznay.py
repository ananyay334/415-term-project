# piv_pressure/tests/test_kovasznay.py
"""Kovasznay flow: second-order convergence including the viscous term."""
import numpy as np
import pytest

from analysis.metrics import convergence_slope, eps_l2
from core.material_derivative import pressure_gradient
from core.pipeline import reconstruct
from fields.analytic import kovasznay, kovasznay_grad_p


@pytest.mark.parametrize("Re", [20.0, 40.0])
def test_gradient_second_order(Re):
    hs, errs = [], []
    for n in [32, 64, 128]:
        fld = kovasznay(n, Re)
        g = pressure_gradient(fld)
        X, Y = fld.grids()
        px, py = kovasznay_grad_p(X, Y, Re)
        hs.append(fld.dx)
        errs.append(max(np.abs(g.Gx - px).max(), np.abs(g.Gy - py).max()))
    assert 1.8 < convergence_slope(hs, errs) < 2.3


@pytest.mark.parametrize("Re", [20.0, 40.0])
@pytest.mark.parametrize("source", ["div_G", "velocity"])
def test_reconstruction_second_order(Re, source):
    hs, errs = [], []
    for n in [32, 64, 128]:
        fld = kovasznay(n, Re)
        rec = reconstruct(fld, source=source)
        hs.append(fld.dx)
        errs.append(eps_l2(rec.p, rec.p_ref))
    assert 1.7 < convergence_slope(hs, errs) < 2.7
    assert errs[-1] < 5e-3


def test_viscous_term_is_needed():
    """Dropping mu lap(u) leaves an O(1) model error that does not converge."""
    errs_full, errs_inviscid = [], []
    for n in [32, 64, 128]:
        fld = kovasznay(n, 20.0)
        errs_full.append(eps_l2(reconstruct(fld).p, fld.p_ref[0]))
        errs_inviscid.append(eps_l2(reconstruct(fld, include_viscous=False).p, fld.p_ref[0]))
    assert errs_inviscid[-1] > 10 * errs_full[-1]
    assert errs_inviscid[-1] > 0.5 * errs_inviscid[0]          # no convergence without it


def test_lagrangian_steady_option_runs():
    fld = kovasznay(64, 40.0)
    rec = reconstruct(fld, material_derivative="lagrangian")
    assert eps_l2(rec.p, rec.p_ref) < 2e-2
