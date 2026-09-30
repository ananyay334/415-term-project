# piv_pressure/tests/test_taylor_green.py
"""Taylor-Green vortex: G = grad p to second order; full reconstruction converges.

Note: G is defined as G = -rho Du/Dt + mu lap u = +grad p, so the test compares
G with +grad p (the spec's wording "-grad p" is a sign slip).
"""
import numpy as np
import pytest

from analysis.metrics import convergence_slope, eps_l2
from core.material_derivative import pressure_gradient
from core.pipeline import reconstruct
from fields.analytic import taylor_green, taylor_green_grad_p, taylor_green_uvp

NS = [32, 64, 128]


@pytest.mark.parametrize("advect, amp", [((0.0, 0.0), 1.0), ((1.0, 0.0), 1.0), ((1.0, 0.0), 0.4)])
def test_gradient_matches_grad_p_second_order(advect, amp):
    hs, errs = [], []
    for n in NS:
        fld = taylor_green(n, advect=advect, dt=0.01, amplitude=amp)
        g = pressure_gradient(fld, scheme="eulerian")
        X, Y = fld.grids()
        px, py = taylor_green_grad_p(X, Y, fld.t[1], advect=advect, amplitude=amp)
        hs.append(fld.dx)
        errs.append(max(np.abs(g.Gx - px).max(), np.abs(g.Gy - py).max()))
    assert 1.8 < convergence_slope(hs, errs) < 2.3


@pytest.mark.parametrize("method", ["poisson", "direct", "omni"])
def test_full_reconstruction_second_order(method):
    ns = NS if method != "omni" else [16, 32, 64]
    hs, errs = [], []
    for n in ns:
        fld = taylor_green(n)
        rec = reconstruct(fld, method=method)
        hs.append(fld.dx)
        errs.append(eps_l2(rec.p, rec.p_ref))
    assert 1.7 < convergence_slope(hs, errs) < 2.4


def test_convected_taylor_green_is_exact_ns_solution():
    """The Galilean-shifted vortex satisfies du/dt + u.grad u = -grad p/rho + nu lap u."""
    nu, c, t, eps = 0.01, (1.0, 0.3), 0.7, 1e-5
    X, Y = np.meshgrid(np.linspace(0.2, 2.0, 7), np.linspace(0.1, 1.9, 7))
    u, v, _ = taylor_green_uvp(X, Y, t, nu, advect=c)
    up, vp, _ = taylor_green_uvp(X, Y, t + eps, nu, advect=c)
    um, vm, _ = taylor_green_uvp(X, Y, t - eps, nu, advect=c)
    d = lambda F, dx_, dy_: (F(X + dx_, Y + dy_) - F(X - dx_, Y - dy_)) / (2 * eps)   # noqa: E731
    U = lambda Xq, Yq: taylor_green_uvp(Xq, Yq, t, nu, advect=c)[0]                   # noqa: E731
    ut = (up - um) / (2 * eps)
    ux, uy = d(U, eps, 0), d(U, 0, eps)
    xi = X - c[0] * t
    eta = Y - c[1] * t
    lap_u = -2.0 * np.cos(xi) * np.sin(eta) * np.exp(-2 * nu * t)
    px, _ = taylor_green_grad_p(X, Y, t, nu, advect=c)
    residual = ut + u * ux + v * uy + px - nu * lap_u
    assert np.abs(residual).max() < 1e-6


def test_lagrangian_second_order_in_dt():
    """Pseudo-tracking error decreases as dt^2 (convected vortex, fine grid)."""
    dts, errs = [0.2, 0.1, 0.05], []
    for dt in dts:
        fld = taylor_green(128, advect=(1.0, 0.0), dt=dt)
        g = pressure_gradient(fld, scheme="lagrangian")
        X, Y = fld.grids()
        px, _ = taylor_green_grad_p(X, Y, fld.t[1], advect=(1.0, 0.0))
        # Central region only: at dt = 0.2 a particle moves up to 0.4 (~23 nodes),
        # and nodes closer to the edge use the Eulerian fall-back.
        errs.append(np.abs(g.Gx - px)[25:-25, 25:-25].max())
        assert g.info["interpolation"] in ("cubic", "linear")
    assert 1.8 < convergence_slope(dts, errs) < 2.3
