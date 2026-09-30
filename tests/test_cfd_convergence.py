# piv_pressure/tests/test_cfd_convergence.py
"""Richardson extrapolation / GCI, OpenFOAM log parsing and time-loop criteria.

The log and signal inputs are small synthetic fixtures in the format the
parsers expect; they test the logic, not OpenFOAM itself.
"""
import numpy as np
import pytest

from cfd.change_norm import change_norms, cycle_statistics, first_below, periodicity_check
from cfd.grid_independence import representative_h, richardson
from cfd.parse_log import parse_log


@pytest.mark.parametrize("p, hs", [(2.0, (0.25, 0.5, 1.0)), (1.0, (0.1, 0.2, 0.4)),
                                   (1.5, (0.3, 0.45, 0.9))])
def test_richardson_exact_on_power_law(p, hs):
    """f(h) = f0 + C h^p  ->  exact order and limit, also for unequal ratios."""
    f0, C = 0.1649, -0.37
    h1, h2, h3 = hs
    f = [f0 + C * h**p for h in hs]
    R = richardson(f[0], f[1], f[2], h2 / h1, h3 / h2)
    assert R["order"] == pytest.approx(p, rel=1e-8)
    assert R["extrapolated"] == pytest.approx(f0, rel=1e-10)
    assert R["monotone"] and not R["oscillatory"]
    assert R["GCI_fine"] < R["GCI_medium"]
    assert R["pct_change_medium_fine"] < R["pct_change_coarse_medium"]


def test_richardson_flags_oscillatory_convergence():
    R = richardson(1.00, 1.02, 0.99)          # e21 > 0, e32 < 0
    assert R["oscillatory"] and np.isnan(R["order"])


def test_representative_h():
    assert representative_h(400, area=4.0) == pytest.approx(0.1)


def test_parse_log_esi_and_foundation(tmp_path):
    text = "\n".join([
        "Starting time loop", "Courant Number mean: 0.01 max: 0.3", "Time = 0.1", "",
        "smoothSolver:  Solving for Ux, Initial residual = 0.5, Final residual = 1e-09, No Iterations 3",
        "smoothSolver:  Solving for Uy, Initial residual = 0.25, Final residual = 1e-09, No Iterations 3",
        "GAMG:  Solving for p, Initial residual = 0.125, Final residual = 1e-05, No Iterations 4",
        "GAMG:  Solving for p, Initial residual = 0.001, Final residual = 1e-08, No Iterations 4",
        "Courant Number mean: 0.02 max: 0.4", "Time = 0.2s", "",
        "smoothSolver:  Solving for Ux, Initial residual = 0.05, Final residual = 1e-09, No Iterations 3",
        "GAMG:  Solving for p, Initial residual = 0.0125, Final residual = 1e-05, No Iterations 4",
    ])
    path = tmp_path / "log.pimpleFoam"
    path.write_text(text)
    res = parse_log(path)
    assert np.allclose(res["t"], [0.1, 0.2])
    assert np.allclose(res["Ux"], [0.5, 0.05]) and np.allclose(res["p"], [0.125, 0.0125])
    assert res["Uy"][0] == 0.25 and np.isnan(res["Uy"][1])
    assert np.allclose(res["Co_max"], [0.3, 0.4])


def test_periodicity_check_detects_settled_cycles():
    t = np.arange(0.0, 200.0, 0.01)
    amp = 0.33 * (1.0 - np.exp(-t / 15.0))                # growing, then settled
    cl = amp * np.sin(2 * np.pi * 0.165 * t)
    res = periodicity_check(t, cl, tol=1e-3, n_consecutive=3)
    assert res["satisfied"]
    assert 50.0 < res["t_periodic"] < 150.0
    c = cycle_statistics(t, cl)
    assert np.allclose(c["period"][-5:], 1 / 0.165, rtol=1e-3)
    assert np.allclose(c["amplitude"][-3:], 0.33, rtol=1e-3)


def test_change_norm_steady_and_periodic():
    t = np.linspace(0, 10, 101)
    X = np.linspace(0, 1, 20)[None, None, :] * np.ones((101, 20, 20))
    decay = np.exp(-t)[:, None, None]
    u_steady = 1.0 + decay * np.sin(3 * X)
    ch = change_norms(u_steady, 0 * u_steady)
    assert np.all(np.diff(ch) < 0)                        # falls monotonically to zero
    k, t_eps = first_below(t[1:], ch, 1e-4)
    assert k >= 0 and ch[k] < 1e-4 <= ch[k - 1]
    u_per = 1.0 + np.sin(2 * np.pi * 0.5 * t)[:, None, None] * np.sin(3 * X)
    ch_p = change_norms(u_per, 0 * u_per)
    assert ch_p[-20:].min() > 1e-3                        # never reaches a steady tolerance
