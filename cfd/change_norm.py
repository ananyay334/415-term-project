# piv_pressure/cfd/change_norm.py
"""Time-loop convergence criteria for the CFD runs (Lecture 11).

1. **Steady-state criterion** ||u^{n+1} - u^n|| < eps: the change between two
   saved snapshots, measured with the normalised L2 norm over the window
   (both velocity components, divided by U_inf).  For a steady flow it falls
   towards zero and "breaks the time loop" once below eps.
2. **Periodicity criterion** for the Re = 100 wake: the flow never becomes
   steady (vortices keep shedding), so ||u^{n+1} - u^n|| oscillates at a
   finite level and can never satisfy a steady-state tolerance.  Instead the
   *cycle-to-cycle* change of the lift signal is used: for each shedding cycle
   (between upward zero crossings of Cl - mean) the period T_k and amplitude
   A_k = (max - min)/2 are measured; the flow is periodic once
   |A_k - A_{k-1}|/A_k and |T_k - T_{k-1}|/T_k stay below ``tol`` for
   ``n_consecutive`` cycles in a row.

Usage:  python cfd/change_norm.py [--tol 1e-3] [--eps 1e-6]
Output: results/cfd/periodicity_<mesh>.csv, change_norm_window.csv,
        steady_re20_change.csv and results/cfd/time_loop_summary.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis.norms import l2_normalised  # noqa: E402
from analysis.plots import RESULTS_DIR  # noqa: E402
from cfd import export_window, grid_independence  # noqa: E402
from fields.cfd_loader import DEFAULT_CFD_NPZ, cfd_available, load_cfd  # noqa: E402

PERIODIC_LEVELS = ("coarse", "medium", "fine")
STEADY_CASE = "steady_re20"


def change_norms(u: np.ndarray, v: np.ndarray, U_ref: float = 1.0) -> np.ndarray:
    """||u^{n+1} - u^n|| (L2_n over both components) / U_ref for each pair of snapshots."""
    du = np.diff(u, axis=0)
    dv = np.diff(v, axis=0)
    return np.array([l2_normalised(np.stack([a, b])) for a, b in zip(du, dv)]) / U_ref


def first_below(t: np.ndarray, values: np.ndarray, eps: float) -> tuple[int, float]:
    """Index and time of the first value below ``eps`` (-1, NaN if never)."""
    idx = np.nonzero(np.asarray(values) < eps)[0]
    return (int(idx[0]), float(t[idx[0]])) if idx.size else (-1, float("nan"))


def cycle_statistics(t: np.ndarray, cl: np.ndarray) -> dict[str, np.ndarray]:
    """Per-cycle start time, period and amplitude of a (nearly) periodic signal.

    Cycles run between consecutive upward zero crossings of cl - mean, where
    the mean is taken over the second half of the record (after start-up).
    Crossing times are refined by linear interpolation.
    """
    t = np.asarray(t, float)
    cl = np.asarray(cl, float)
    s = cl - np.mean(cl[t >= t[0] + 0.5 * (t[-1] - t[0])])
    up = np.nonzero((s[:-1] < 0.0) & (s[1:] >= 0.0))[0]
    tc = t[up] - s[up] * (t[up + 1] - t[up]) / (s[up + 1] - s[up])
    starts, periods, amps = [], [], []
    for a, b, k0, k1 in zip(tc[:-1], tc[1:], up[:-1], up[1:]):
        seg = cl[k0 + 1:k1 + 1]
        if seg.size < 3:
            continue
        starts.append(a)
        periods.append(b - a)
        amps.append(0.5 * (seg.max() - seg.min()))
    return {"t_start": np.array(starts), "period": np.array(periods), "amplitude": np.array(amps)}


def periodicity_check(t: np.ndarray, cl: np.ndarray, tol: float = 1e-3,
                      n_consecutive: int = 3) -> dict:
    """Cycle-to-cycle periodicity test of the lift coefficient.

    Returns the per-cycle arrays, the relative changes, whether the criterion
    is met and the time from which it holds (start of the first cycle of the
    first run of ``n_consecutive`` cycles that all satisfy it).
    """
    c = cycle_statistics(t, cl)
    A, T = c["amplitude"], c["period"]
    dA = np.abs(np.diff(A)) / A[1:] if A.size > 1 else np.array([])
    dT = np.abs(np.diff(T)) / T[1:] if T.size > 1 else np.array([])
    ok = (dA < tol) & (dT < tol)
    t_periodic, satisfied = float("nan"), False
    for k in range(ok.size - n_consecutive + 1):
        if ok[k:k + n_consecutive].all():
            t_periodic, satisfied = float(c["t_start"][k + 1]), True
            break
    return {**c, "rel_change_amplitude": dA, "rel_change_period": dT, "satisfied": satisfied,
            "t_periodic": t_periodic, "tol": tol, "n_consecutive": n_consecutive,
            "n_cycles": int(A.size)}


def write_periodicity_csv(res: dict, path: Path) -> Path:
    """One row per shedding cycle."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("cycle,t_start,period,amplitude,rel_change_period,rel_change_amplitude\n")
        for k in range(res["n_cycles"]):
            dT = res["rel_change_period"][k - 1] if k > 0 else float("nan")
            dA = res["rel_change_amplitude"][k - 1] if k > 0 else float("nan")
            fh.write(f"{k},{res['t_start'][k]:.6g},{res['period'][k]:.6g},{res['amplitude'][k]:.6g},"
                     f"{dT:.6g},{dA:.6g}\n")
    return path


def write_change_csv(t: np.ndarray, change: np.ndarray, path: Path) -> Path:
    """Time (of the later snapshot) and change norm."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("t,change_l2n\n")
        for a, b in zip(t, change):
            fh.write(f"{a:.6g},{b:.6g}\n")
    return path


def main(argv: list[str] | None = None) -> dict:
    """Run every check for which data exist; write CSVs and a JSON summary."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", type=Path, default=ROOT / "cfd" / "cases")
    ap.add_argument("--outdir", type=Path, default=RESULTS_DIR / "cfd")
    ap.add_argument("--tol", type=float, default=1e-3, help="cycle-to-cycle periodicity tolerance")
    ap.add_argument("--n-consecutive", type=int, default=3)
    ap.add_argument("--eps", type=float, default=1e-6, help="steady-state tolerance on the change norm")
    args = ap.parse_args(argv)
    summary: dict = {"periodicity": {}, "tol": args.tol, "eps": args.eps}
    for level in PERIODIC_LEVELS:
        files = grid_independence.coefficient_files(args.cases / level)
        if not files:
            continue
        t, _, cl = grid_independence.read_coefficients(files)
        res = periodicity_check(t, cl, args.tol, args.n_consecutive)
        write_periodicity_csv(res, args.outdir / f"periodicity_{level}.csv")
        summary["periodicity"][level] = {k: res[k] for k in ("satisfied", "t_periodic", "n_cycles")}
        print(f"[periodicity] {level}: {res['n_cycles']} cycles, periodic from t = {res['t_periodic']:.4g} "
              f"(satisfied={res['satisfied']})")
    if cfd_available():
        fld = load_cfd(DEFAULT_CFD_NPZ)
        ch = change_norms(fld.u, fld.v, fld.U_ref)
        write_change_csv(fld.t[1:], ch, args.outdir / "change_norm_window.csv")
        summary["window_change"] = {"min": float(ch.min()), "max": float(ch.max()),
                                    "mean": float(ch.mean())}
        print(f"[change] periodic wake window: ||u^(n+1)-u^n|| between {ch.min():.3g} and {ch.max():.3g}")
    steady = args.cases / STEADY_CASE
    if (steady / "postProcessing" / "sampleWindow").is_dir():
        x, y, t, U, V, P, info = export_window.load_openfoam_samples(steady)
        ch = change_norms(U, V, float(info.get("U_inf", 1.0)))
        write_change_csv(t[1:], ch, args.outdir / "steady_re20_change.csv")
        k, t_eps = first_below(t[1:], ch, args.eps)
        summary["steady"] = {"t_below_eps": t_eps, "t_stop": float(t[-1]), "final_change": float(ch[-1]),
                             "eps": args.eps, "sample_dt": float(t[1] - t[0])}
        print(f"[steady] change < {args.eps:g} first at t = {t_eps:.4g}; OpenFOAM stopped at t = {t[-1]:.4g}")
    if not (summary["periodicity"] or "window_change" in summary or "steady" in summary):
        print("[skip] no CFD force coefficients, window export or steady-case samples found.")
        return summary
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "time_loop_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
