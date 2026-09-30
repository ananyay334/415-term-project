# piv_pressure/tests/test_iterative.py
"""Iterative Poisson solvers match spsolve and show the expected behaviour."""
import numpy as np
import pytest

from core import iterative
from core.pipeline import reconstruct
from core.poisson import linear_system, poisson_system, solve_poisson
from fields.analytic import taylor_green

X_RANGE, Y_RANGE = (0.3, 2.5), (0.1, 2.0)


def manufactured(n=33, m=29):
    """p = sin x cos y with exact Neumann data (as in test_poisson_exact.py)."""
    x, y = np.linspace(*X_RANGE, n), np.linspace(*Y_RANGE, m)
    X, Y = np.meshgrid(x, y)
    p = np.sin(X) * np.cos(Y)
    return x, y, p, -2.0 * p, np.cos(X) * np.cos(Y), -np.sin(X) * np.sin(Y)


@pytest.mark.parametrize("solver", iterative.METHODS)
@pytest.mark.parametrize("pin", ["mean", "node"])
def test_solvers_match_spsolve(solver, pin):
    x, y, p, f, Gx, Gy = manufactured()
    dx, dy = x[1] - x[0], y[1] - y[0]
    pd, _ = solve_poisson(f, Gx, Gy, dx, dy, pin=pin)
    pi, info = solve_poisson(f, Gx, Gy, dx, dy, pin=pin, solver=solver, tol=1e-10)
    assert info["converged"], solver
    assert np.max(np.abs((pi - pi.mean()) - (pd - pd.mean()))) < 1e-6
    h = info["history"]
    assert h.residual_scaled[-1] < 1e-10
    assert h.residual_scaled[-1] < 1e-6 * h.residual_scaled[0]      # decreases overall
    assert h.change_l2n.size == h.iterations


def test_dirichlet_edges_with_iterative_solvers():
    x, y, p, f, Gx, Gy = manufactured()
    dx, dy = x[1] - x[0], y[1] - y[0]
    kw = dict(dirichlet_edges=("top", "bottom"), dirichlet_values=p)
    pd, _ = solve_poisson(f, Gx, Gy, dx, dy, **kw)
    for s in iterative.METHODS:
        pi, info = solve_poisson(f, Gx, Gy, dx, dy, solver=s, tol=1e-10, **kw)
        assert np.allclose(pi[-1], p[-1]) and np.max(np.abs(pi - pd)) < 1e-6, s


def test_jacobi_is_explicit_pseudo_time_marching():
    """One Jacobi step equals p + dtau (lap_h p - f) with dtau = omega h^2/4 (dx = dy)."""
    fld = taylor_green(16)
    rec = reconstruct(fld)                     # just to get a consistent G
    from core.poisson import source_div_G
    f = source_div_G(rec.Gx, rec.Gy, fld.dx, fld.dy)
    A, b, w, fixed, pf, _ = poisson_system(f, rec.Gx, rec.Gy, fld.dx, fld.dy)
    system, _ = linear_system(A, b, w, fixed, pf, (16, 16))
    p0 = np.random.default_rng(0).standard_normal(256)
    p0 -= p0.mean()
    system.p0 = p0
    p1, _ = iterative.jacobi(system, omega=0.8, max_iter=1)
    dtau = 0.8 * fld.dx**2 / 4.0
    expected = p0 + dtau * (A @ p0 - system.b)
    assert np.allclose(p1, expected - expected.mean())


def test_classic_jacobi_stalls_on_all_neumann_problem():
    """omega = 1 leaves the checkerboard mode undamped (amplification -1)."""
    fld = taylor_green(24)
    r1 = reconstruct(fld, solver="jacobi", omega=1.0, max_iter=6000)
    r2 = reconstruct(fld, solver="jacobi", max_iter=6000)          # default omega = 0.95
    h1, h2 = r1.extra["history"], r2.extra["history"]
    assert not h1.converged and h2.converged
    tail = h1.residual_scaled[-500:]
    assert tail.max() / tail.min() < 1.01                            # levelled off


def test_sor_faster_and_diverges_for_omega_ge_2():
    fld = taylor_green(32)
    it = {s: reconstruct(fld, solver=s).extra["history"].iterations for s in ("jacobi", "gs", "sor", "cg")}
    assert it["sor"] < it["gs"] / 10 and it["cg"] < it["gs"] / 10 and it["gs"] < it["jacobi"]
    bad = reconstruct(fld, solver="sor", omega=2.1, max_iter=400).extra["history"]
    assert bad.diverged or bad.residual_scaled[-1] > bad.residual_scaled[0]


def test_omega_opt_formulas():
    n, h = 65, 2.2 / 64
    assert iterative.omega_opt(n, n, h, h, neumann=False) == pytest.approx(
        iterative.omega_opt_dirichlet_formula(h, 64 * h), rel=1e-12)
    assert 1.0 < iterative.omega_opt_dirichlet_formula(h, 64 * h) < iterative.omega_opt(n, n, h, h, True) < 2.0


def test_history_rows_and_physical_scale():
    fld = taylor_green(16)
    r = reconstruct(fld, solver="cg", residual_scale="physical", L_ref=2.0)
    h = r.extra["history"]
    assert h.scale == pytest.approx(1.0 * 1.0**2 / 2.0**2)
    rows = h.to_rows(every=5)
    assert rows[0]["iteration"] == 0 and rows[-1]["iteration"] == h.iterations
