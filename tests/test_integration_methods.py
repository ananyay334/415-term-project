# piv_pressure/tests/test_integration_methods.py
"""Direct, Poisson and omni-directional integration agree on clean data."""
import numpy as np
import pytest

from analysis.metrics import demean, eps_l2
from core.direct_integration import CORNERS, direct_integration
from core.omni_integration import omni_integration
from core.pipeline import reconstruct
from core.poisson import poisson_pressure
from fields.analytic import kovasznay, taylor_green, taylor_green_grad_p, taylor_green_uvp


def linear_field(nx=21, ny=17):
    x, y = np.linspace(0.0, 1.0, nx), np.linspace(-0.5, 0.5, ny)
    X, Y = np.meshgrid(x, y)
    return x, y, 1.5 * X - 2.5 * Y, 1.5 + 0 * X, -2.5 + 0 * X


def test_all_methods_exact_for_linear_pressure():
    x, y, p, Gx, Gy = linear_field()
    for c in CORNERS:
        for q in direct_integration(Gx, Gy, x, y, corner=c).values():
            assert np.abs(demean(q) - demean(p)).max() < 1e-12
    po, _ = omni_integration(Gx, Gy, x, y, n_directions=12)
    assert np.abs(demean(po) - demean(p)).max() < 1e-10
    pp, _ = poisson_pressure(Gx, Gy, x, y)
    assert np.abs(demean(pp) - demean(p)).max() < 1e-10


def test_exact_gradient_integrates_to_exact_pressure():
    fld = taylor_green(64)
    X, Y = fld.grids()
    gx, gy = taylor_green_grad_p(X, Y, fld.t[1])
    p = fld.p_ref[1]
    d = direct_integration(gx, gy, fld.x, fld.y)
    assert eps_l2(d["path1"], p) < 2e-3 and eps_l2(d["path2"], p) < 2e-3
    po, info = omni_integration(gx, gy, fld.x, fld.y)
    assert eps_l2(po, p) < 2e-3
    assert info["direction_spread_rms"] < 1e-2      # nearly path independent for an exact gradient


@pytest.mark.parametrize(
    "make, tol",
    [(lambda: taylor_green(64), 5e-3),
     # Kovasznay has a steep pressure rise near x = -0.5; a single integration
     # path accumulates the O(h^2) error of G there, so the tolerance is wider.
     (lambda: kovasznay(64, 40.0), 3e-2)],
)
def test_methods_agree_on_clean_analytic_data(make, tol):
    fld = make()
    recs = {m: reconstruct(fld, method=m) for m in ("poisson", "direct", "omni")}
    for m, r in recs.items():
        assert eps_l2(r.p, r.p_ref) < tol, m
    ref = recs["poisson"].p
    for m in ("direct", "omni"):
        assert eps_l2(recs[m].p, ref) < tol, m


def test_path_dependence_reported_for_noisy_gradient():
    fld = taylor_green(48)
    X, Y = fld.grids()
    gx, gy = taylor_green_grad_p(X, Y, fld.t[1])
    rng = np.random.default_rng(0)
    gx, gy = gx + 0.05 * rng.standard_normal(gx.shape), gy + 0.05 * rng.standard_normal(gy.shape)
    d = direct_integration(gx, gy, fld.x, fld.y)
    assert np.abs(demean(d["path1"]) - demean(d["path2"])).max() > 1e-3
    # omni and Poisson average the noise: both beat a single path
    po, _ = omni_integration(gx, gy, fld.x, fld.y)
    pp, _ = poisson_pressure(gx, gy, fld.x, fld.y)
    p = fld.p_ref[1]
    assert eps_l2(po, p) < eps_l2(d["path1"], p)
    assert eps_l2(pp, p) < eps_l2(d["path1"], p)


def test_bernoulli_dirichlet_runs():
    fld = taylor_green(32)
    rec = reconstruct(fld, bc="bernoulli", dirichlet_edges=("top",))
    u, v, _ = taylor_green_uvp(*fld.grids(), fld.t[1])
    expected_top = 0.5 * (1.0 - u[-1] ** 2 - v[-1] ** 2)
    assert np.allclose(rec.p[-1], expected_top)
    assert rec.diagnostics["bc"] == "dirichlet:top"
