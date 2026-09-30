# piv_pressure/tests/test_norms.py
"""Vector norms, normalised norms, scaled residuals and the refactored metrics."""
import numpy as np
import pytest

from analysis.metrics import eps_l2, eps_max, error_norms, standard_metrics
from analysis.norms import (l1, l1_normalised, l2, l2_normalised, linf, norm_summary,
                            residual_scale, scaled)


def test_lecture_example():
    """Lecture 11 slide: r = (0.0001, 0.0004, 0.0006, 0.002)."""
    r = np.array([0.0001, 0.0004, 0.0006, 0.002])
    assert l1(r) == pytest.approx(0.0031)
    assert l2(r) == pytest.approx(0.00213, abs=5e-6)              # sqrt(4.53e-6) = 0.002128
    assert l2_normalised(r) == pytest.approx(0.0011, abs=5e-5)    # 0.001064 ~ "representative" 0.0011
    assert linf(r) == pytest.approx(0.002)
    assert l1_normalised(r) == pytest.approx(0.0031 / 4)


def test_known_vector_and_shapes():
    r = np.array([[3.0, -4.0], [0.0, 0.0]])
    s = norm_summary(r)
    assert s == pytest.approx({"l1": 7.0, "l2": 5.0, "linf": 4.0, "l1n": 1.75, "l2n": 2.5})


@pytest.mark.parametrize("n", [4, 100, 10_000])
def test_normalised_norms_independent_of_length(n):
    c = 0.37
    r = np.full(n, c)
    assert l1_normalised(r) == pytest.approx(c) and l2_normalised(r) == pytest.approx(c)
    assert l1(r) == pytest.approx(n * c) and l2(r) == pytest.approx(np.sqrt(n) * c)


def test_norm_ordering():
    """For any vector: L2n <= Linf and L1n <= L2n (power-mean inequality)."""
    r = np.random.default_rng(1).standard_normal(500)
    assert l1_normalised(r) <= l2_normalised(r) <= linf(r)


def test_scaled_and_scales():
    b = np.array([2.0, -2.0, 2.0, -2.0])
    assert residual_scale(b, "rhs") == pytest.approx(2.0)
    assert residual_scale(kind="physical", rho=1.2, U=2.0, L=0.5) == pytest.approx(1.2 * 4 / 0.25)
    assert np.allclose(scaled(b, 2.0), b / 2.0)
    with pytest.raises(ValueError):
        scaled(b, 0.0)
    assert residual_scale(np.zeros(3), "rhs") == 1.0


def test_metrics_use_norms_consistently():
    rng = np.random.default_rng(3)
    p, q = rng.random((20, 30)), rng.random((20, 30))
    e = (p - p.mean()) - (q - q.mean())
    assert eps_l2(p, q) == pytest.approx(np.linalg.norm(e) / np.linalg.norm(q - q.mean()), rel=1e-15)
    assert eps_max(p, q, 1.0, 1.0) == pytest.approx(np.abs(e).max() / 0.5, rel=1e-15)
    en = error_norms(p, q)
    assert en["err_pct"] == pytest.approx(100 * eps_l2(p, q))
    assert en["err_l2n"] == pytest.approx(np.sqrt(np.mean(e**2)))
    assert set(standard_metrics(p, q)) >= {"eps_l2", "eps_max", "err_l1n", "err_l2n", "err_linf", "err_pct"}
