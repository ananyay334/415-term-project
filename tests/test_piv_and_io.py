# piv_pressure/tests/test_piv_and_io.py
"""PIV degradation, denoising, CFD file round-trip and OpenFOAM parsers.

(Beyond the five required test files: these guard the data path of the
studies.  The OpenFOAM files here are small synthetic fixtures written in
the format the parsers expect - they test the parsing logic, not OpenFOAM.)
"""
import json

import numpy as np
import pytest

from cfd import export_window, generate_cases, grid_independence
from fields.analytic import taylor_green_on_grid
from fields.cfd_loader import load_cfd, save_cfd_npz
from piv.degrade import add_noise, box_filter, degrade, degrade_analytic
from piv.denoise import denoise_field, savitzky_golay_2d


def test_box_filter_centred_and_normalised():
    d = np.zeros((15, 15))
    d[7, 7] = 1.0
    for w in (2, 3, 4, 5):
        k = box_filter(d, w)
        assert k.sum() == pytest.approx(1.0)
        assert np.allclose(k, k[::-1, ::-1])                  # centred (no half-cell shift)
    lin = np.add.outer(np.arange(15.0), 2 * np.arange(15.0))
    assert np.allclose(box_filter(lin, 4)[4:-4, 4:-4], lin[4:-4, 4:-4])   # preserves linear fields


def test_noise_statistics_and_seed():
    fld = taylor_green_on_grid(np.linspace(0, 1, 64), np.linspace(0, 1, 64), np.array([0.0, 0.1, 0.2]))
    a, b, c = add_noise(fld, 0.02, seed=1), add_noise(fld, 0.02, seed=1), add_noise(fld, 0.02, seed=2)
    assert np.array_equal(a.u, b.u) and not np.array_equal(a.u, c.u)
    assert np.std(a.u - fld.u) == pytest.approx(0.02, rel=0.05)


def test_degrade_analytic_grid_and_reference():
    truth = lambda x, y: taylor_green_on_grid(x, y, np.array([0.95, 1.0, 1.05]))   # noqa: E731
    f = degrade_analytic(truth, (0.3, 2.5), (0.3, 2.5), 32, window_h=2.0, sigma=0.0)
    assert f.nx == f.ny == 32 and np.allclose(f.x, np.linspace(0.3, 2.5, 32))
    exact = truth(f.x, f.y)
    assert np.allclose(f.p_ref, exact.p_ref)                  # reference pressure never filtered
    assert 0 < np.abs(f.u - exact.u).max() < 0.05             # velocity is filtered


def test_degrade_subsample_and_crop():
    fld = taylor_green_on_grid(np.linspace(0, 2, 41), np.linspace(0, 2, 41), np.array([0.0, 0.1, 0.2]))
    out = degrade(fld, window=1, stride=4, crop=1)
    assert out.nx == 9 and out.x[0] == pytest.approx(0.2)


def test_savgol_preserves_quadratics_everywhere():
    x = np.linspace(-1, 1, 20)
    X, Y = np.meshgrid(x, x)
    q = 1 + 2 * X - Y + 0.5 * X**2 + X * Y - 3 * Y**2
    assert np.allclose(savitzky_golay_2d(q, 5, 2), q, atol=1e-10)


def test_denoisers_reduce_noise():
    x = np.linspace(0.3, 2.5, 64)
    fld = taylor_green_on_grid(x, x, np.array([0.95, 1.0, 1.05]))
    noisy = add_noise(fld, 0.05, seed=0)
    base = np.std(noisy.u - fld.u)
    for method, size in (("gaussian", 1.0), ("savgol", 5)):
        assert np.std(denoise_field(noisy, method, size).u - fld.u) < 0.6 * base


def test_cfd_npz_round_trip(tmp_path):
    x, y, t = np.linspace(1, 5, 9), np.linspace(-2, 2, 9), np.arange(4) * 0.05
    u = np.random.default_rng(0).random((4, 9, 9))
    path = save_cfd_npz(tmp_path / "w.npz", x, y, t, u, 2 * u, 3 * u, nu=0.01, meta={"mesh": "test"})
    fld = load_cfd(path)
    assert np.allclose(fld.v, 2 * u) and np.allclose(fld.p_ref, 3 * u)
    assert fld.mu == pytest.approx(0.01) and fld.meta["mesh"] == "test"


def test_export_parser_handles_both_file_layouts(tmp_path):
    info = {"window": {"x0": 1.0, "x1": 2.0, "y0": 0.0, "y1": 1.0, "nx": 5, "ny": 4}}
    (tmp_path / "mesh_info.json").write_text(json.dumps(info))
    x, y = export_window.window_grid(info)
    X, Y = np.meshgrid(x, y)
    perm = np.random.default_rng(0).permutation(X.size)
    xyz = np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)][perm]
    for k, t in enumerate((10.0, 10.05, 10.1)):
        d = tmp_path / "postProcessing" / "sampleWindow" / f"{t:g}"
        d.mkdir(parents=True)
        p = (X + k).ravel()[perm]
        if k == 1:   # combined layout in one directory, separate files in the others
            np.savetxt(d / "window_p_U.xy", np.c_[xyz, p, 2 * p, 3 * p, 0 * p])
        else:
            np.savetxt(d / "window_p.xy", np.c_[xyz, p])
            np.savetxt(d / "window_U.xy", np.c_[xyz, 2 * p, 3 * p, 0 * p])
    xg, yg, t, U, V, P, _ = export_window.load_openfoam_samples(tmp_path)
    assert np.allclose(t, [10.0, 10.05, 10.1])
    assert np.allclose(P[1], X + 1) and np.allclose(U[2], 2 * (X + 2)) and np.allclose(V[0], 3 * X)


def test_force_coefficient_analysis(tmp_path):
    d = tmp_path / "postProcessing" / "forceCoeffs1" / "0"
    d.mkdir(parents=True)
    t = np.arange(0, 150, 0.004)
    cl = 0.33 * np.sin(2 * np.pi * 0.165 * t)
    cd = 1.33 + 0.01 * np.sin(4 * np.pi * 0.165 * t)
    body = "\n".join(f"{a:.5f} {b:.6f} {c:.6f}" for a, b, c in zip(t, cd, cl))
    (d / "forceCoeffs.dat").write_text("# Force coefficients\n# Time Cd Cl\n" + body + "\n")
    tt, cdd, cll = grid_independence.read_coefficients(grid_independence.coefficient_files(tmp_path))
    res = grid_independence.analyse(tt, cdd, cll, window=80.0)
    assert res["St"] == pytest.approx(0.165, abs=5e-4)
    assert res["Cd_mean"] == pytest.approx(1.33, abs=2e-3)
    assert res["Cl_amp"] == pytest.approx(0.33, abs=2e-3)


def test_richardson_recovers_second_order():
    f_exact, c = 1.0, 0.1
    vals = [f_exact + c * h**2 for h in (0.25, 0.5, 1.0)]     # fine, medium, coarse
    r = grid_independence.richardson(*vals)
    assert r["order"] == pytest.approx(2.0) and r["extrapolated"] == pytest.approx(f_exact)


def test_block_mesh_topology_is_consistent():
    for level in generate_cases.LEVELS:
        generate_cases.check_topology(level)                  # raises on any inconsistency
        assert level.n_cells() > 0
