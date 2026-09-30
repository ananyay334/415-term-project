# piv_pressure/core/pipeline.py
"""One entry point for the whole reconstruction: velocity -> G -> pressure.

    field --(material derivative + viscous term)--> G = grad p
          --(Poisson | direct | omni integration)--> p

Example
-------
>>> from fields.analytic import taylor_green
>>> from core.pipeline import reconstruct
>>> rec = reconstruct(taylor_green(64), method="poisson")
>>> rec.p.shape
(64, 64)
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from analysis.norms import residual_scale as residual_scale_value
from core.direct_integration import direct_integration
from core.iterative import DEFAULT_MAX_ITER, DEFAULT_TOL
from core.material_derivative import pressure_gradient
from core.omni_integration import omni_integration
from core.poisson import bernoulli_pressure, poisson_pressure
from fields.field import FlowField

METHODS: tuple[str, ...] = ("poisson", "direct", "omni")
BCS: tuple[str, ...] = ("neumann", "bernoulli", "exact")


@dataclass
class Reconstruction:
    """Result of :func:`reconstruct`.

    Attributes
    ----------
    p : reconstructed pressure (ny, nx); defined up to a constant.
    p_ref : reference pressure at the same instant (or ``None``).
    Gx, Gy : pressure gradient used for the integration.
    time_index : snapshot index that was reconstructed.
    method : integration method name.
    diagnostics : solver/compatibility/timing information.
    extra : additional fields (e.g. the two direct-integration paths, or
        ``"history"`` - the IterationHistory of an iterative Poisson solve).
    """

    p: np.ndarray
    p_ref: np.ndarray | None
    Gx: np.ndarray
    Gy: np.ndarray
    time_index: int
    method: str
    diagnostics: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)


def reconstruct(
    fld: FlowField,
    method: str = "poisson",
    material_derivative: str = "eulerian",
    bc: str = "neumann",
    pin: str = "mean",
    time_index: int | None = None,
    source: str = "div_G",
    solver: str = "direct",
    dirichlet_edges: tuple[str, ...] = ("top", "bottom"),
    U_inf: float | None = None,
    p_inf: float = 0.0,
    n_directions: int = 32,
    direct_path: str = "average",
    corner: str = "SW",
    interpolation: str = "auto",
    include_viscous: bool = True,
    order: int = 2,
    tol: float = DEFAULT_TOL,
    max_iter: int = DEFAULT_MAX_ITER,
    omega: float | None = None,
    residual_scale: str | float = "rhs",
    L_ref: float = 1.0,
) -> Reconstruction:
    """Reconstruct the pressure field at one snapshot.

    Parameters
    ----------
    fld : velocity data (common interface).
    method : ``"poisson"`` (main), ``"direct"`` or ``"omni"``.
    material_derivative : ``"eulerian"`` or ``"lagrangian"``.
    bc : Poisson boundary treatment - ``"neumann"`` on all edges,
        ``"bernoulli"`` (Dirichlet from Bernoulli on ``dirichlet_edges``) or
        ``"exact"`` (Dirichlet from ``fld.p_ref``; an oracle used as a control
        in the BC study, standing in for e.g. wall pressure taps).
    pin : ``"mean"`` or ``"node"`` for the all-Neumann Poisson problem.
    time_index : snapshot to reconstruct (default: central snapshot).
    source : ``"div_G"`` or ``"velocity"`` (2-D incompressible identity).
    solver : Poisson solver - ``"direct"`` (spsolve, default), ``"jacobi"``,
        ``"gs"``, ``"sor"`` or ``"cg"`` (same system, residual history recorded).
    dirichlet_edges : edges used for the Dirichlet options.
    U_inf, p_inf : free-stream values for Bernoulli (U_inf default: fld.U_ref).
    n_directions : ray directions for the omni method.
    direct_path : ``"path1"``, ``"path2"`` or ``"average"`` for the direct method.
    corner : reference corner for the direct method.
    interpolation : interpolation for the Lagrangian scheme.
    include_viscous : keep the mu lap(u) term in G.
    order : spatial order of all finite differences and of the Neumann BCs
        (2 = default; 1 = first order, for the grid-convergence comparison).
    tol, max_iter : iterative stopping rule (scaled L2_n residual < tol).
    omega : relaxation factor (Jacobi default 0.95, SOR default omega_opt).
    residual_scale : ``"rhs"`` (S = L2_n of the right-hand side),
        ``"physical"`` (S = rho U_ref^2 / L_ref^2) or a number.
    L_ref : length scale for the physical residual scale.
    """
    t0 = time.perf_counter()
    grad = pressure_gradient(fld, time_index, material_derivative, include_viscous, interpolation,
                             order)
    n = grad.time_index
    t1 = time.perf_counter()
    diag: dict = {"material_derivative": material_derivative, **grad.info,
                  "time_gradient_s": t1 - t0}
    extra: dict = {}

    if method == "poisson":
        if bc == "neumann":
            edges: tuple[str, ...] = ()
            values = None
        elif bc == "bernoulli":
            edges = dirichlet_edges
            values = bernoulli_pressure(
                fld.u[n], fld.v[n], fld.rho, fld.U_ref if U_inf is None else U_inf, p_inf
            )
        elif bc == "exact":
            if fld.p_ref is None:
                raise ValueError("bc='exact' needs a reference pressure")
            edges, values = dirichlet_edges, fld.p_ref[n]
        else:
            raise ValueError(f"unknown bc {bc!r}; use one of {BCS}")
        if residual_scale == "rhs":
            S = None
        elif isinstance(residual_scale, str):
            S = residual_scale_value(kind=residual_scale, rho=fld.rho, U=fld.U_ref, L=L_ref)
        else:
            S = float(residual_scale)
        p, info = poisson_pressure(
            grad.Gx, grad.Gy, fld.x, fld.y, source=source, u=fld.u[n], v=fld.v[n],
            rho=fld.rho, dirichlet_edges=edges, dirichlet_values=values, pin=pin, solver=solver,
            tol=tol, max_iter=max_iter, order=order, omega=omega, residual_scale=S,
        )
        if "history" in info:
            extra["history"] = info.pop("history")
        diag.update(info)
    elif method == "direct":
        paths = direct_integration(grad.Gx, grad.Gy, fld.x, fld.y, corner)
        if direct_path not in paths:
            raise ValueError(f"direct_path must be one of {tuple(paths)}")
        p = paths[direct_path]
        extra.update(paths)
        d = (paths["path1"] - paths["path1"].mean()) - (paths["path2"] - paths["path2"].mean())
        diag["path_rms_difference"] = float(np.sqrt(np.mean(d**2)))
        diag["direct_path"] = direct_path
    elif method == "omni":
        p, info = omni_integration(grad.Gx, grad.Gy, fld.x, fld.y, n_directions)
        diag.update(info)
    else:
        raise ValueError(f"unknown method {method!r}; use one of {METHODS}")

    diag["time_integration_s"] = time.perf_counter() - t1
    diag["time_total_s"] = time.perf_counter() - t0
    diag["method"] = method
    p_ref = None if fld.p_ref is None else fld.p_ref[n]
    return Reconstruction(p=p, p_ref=p_ref, Gx=grad.Gx, Gy=grad.Gy, time_index=n,
                          method=method, diagnostics=diag, extra=extra)
