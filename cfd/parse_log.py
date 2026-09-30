# piv_pressure/cfd/parse_log.py
"""Extract the solver residuals from OpenFOAM run logs.

For every time step, OpenFOAM prints lines such as::

    Time = 0.0125
    Courant Number mean: 0.0123 max: 0.487
    smoothSolver:  Solving for Ux, Initial residual = 0.0021, Final residual = 1.2e-09, No Iterations 3
    smoothSolver:  Solving for Uy, Initial residual = 0.0154, Final residual = 4.1e-09, No Iterations 3
    GAMG:  Solving for p, Initial residual = 0.0347, Final residual = 0.00031, No Iterations 4
    ...

The *initial* residual is the normalised imbalance of the equation before
the linear solver starts, i.e. how far the previous solution is from
satisfying the new equations.  In a steady problem it falls towards zero as
the time loop converges; in the periodic cylinder wake it levels off,
because the solution keeps changing from step to step.

Per time step the first "Solving for" line of each field is kept (``p`` is
solved several times per step; the first solve is the one foamLog calls
``p_0``).  pimpleFoam prints the Courant number (and deltaT) *before* the
"Time =" line of the step it belongs to, so it is held until that line.  Both the ESI format (``Time = 0.0125``) and the Foundation format
(``Time = 0.0125s``) are accepted [verify for your version].  Output of the
``foamLog`` utility (``logs/Ux_0`` etc., two columns: time, value) can be
read with :func:`read_foamlog_dir`.

Usage:  python cfd/parse_log.py [--cases cfd/cases]   -> results/cfd/residuals_<case>.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis.plots import RESULTS_DIR  # noqa: E402

FIELDS: tuple[str, ...] = ("Ux", "Uy", "p")
_TIME = re.compile(r"^Time = ([-+0-9.eE]+)s?\s*$")
_SOLVE = re.compile(r"Solving for (\w+), Initial residual = ([-+0-9.eE]+), "
                    r"Final residual = ([-+0-9.eE]+), No Iterations (\d+)")
_COURANT = re.compile(r"^Courant Number mean: ([-+0-9.eE]+) max: ([-+0-9.eE]+)")


def parse_log(path: str | Path, fields: tuple[str, ...] = FIELDS) -> dict[str, np.ndarray]:
    """Initial residuals per time step (first solve of each field) and max Co.

    Returns a dict with ``t`` and one array per field (NaN where a field was
    not solved in a step) plus ``Co_max``.
    """
    t: list[float] = []
    data: dict[str, list[float]] = {f: [] for f in fields}
    co: list[float] = []
    current: dict[str, float] | None = None
    co_now = np.nan
    co_pending = np.nan

    def flush() -> None:
        if current is not None:
            for f in fields:
                data[f].append(current.get(f, np.nan))
            co.append(co_now)

    with Path(path).open(errors="replace") as fh:
        for line in fh:
            m = _COURANT.match(line.strip())
            if m:                       # belongs to the next "Time =" step
                co_pending = float(m.group(2))
                continue
            m = _TIME.match(line.strip())
            if m:
                flush()
                t.append(float(m.group(1)))
                current, co_now, co_pending = {}, co_pending, np.nan
                continue
            if current is None:
                continue
            m = _SOLVE.search(line)
            if m and m.group(1) in fields and m.group(1) not in current:
                current[m.group(1)] = float(m.group(2))
    flush()
    out = {"t": np.array(t), "Co_max": np.array(co)}
    out.update({f: np.array(v) for f, v in data.items()})
    return out


def read_foamlog_dir(logs: str | Path, fields: tuple[str, ...] = FIELDS) -> dict[str, np.ndarray]:
    """Read ``foamLog`` output (files ``<field>_0`` with columns time, value)."""
    logs = Path(logs)
    out: dict[str, np.ndarray] = {}
    for f in fields:
        path = logs / f"{f}_0"
        if path.is_file():
            arr = np.loadtxt(path, ndmin=2)
            out.setdefault("t", arr[:, 0])
            out[f] = np.interp(out["t"], arr[:, 0], arr[:, 1])
    if "t" not in out:
        raise FileNotFoundError(f"no foamLog files (e.g. Ux_0) in {logs}")
    return out


def case_residuals(case: Path) -> dict[str, np.ndarray] | None:
    """Residuals of one case from ``log.pimpleFoam`` (or ``logs/`` from foamLog)."""
    log = case / "log.pimpleFoam"
    if log.is_file():
        return parse_log(log)
    if (case / "logs").is_dir():
        return read_foamlog_dir(case / "logs")
    return None


def write_residual_csv(res: dict[str, np.ndarray], path: Path, max_rows: int = 20000) -> Path:
    """CSV of time and residuals, thinned to at most ``max_rows`` rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    n = res["t"].size
    idx = np.unique(np.linspace(0, n - 1, min(n, max_rows)).astype(int)) if n else np.array([], int)
    keys = [k for k in ("t",) + FIELDS + ("Co_max",) if k in res]
    with path.open("w") as fh:
        fh.write(",".join(keys) + "\n")
        for i in idx:
            fh.write(",".join(f"{res[k][i]:.6g}" for k in keys) + "\n")
    return path


def main(argv: list[str] | None = None) -> dict[str, dict]:
    """Parse the logs of every case directory that has one."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", type=Path, default=ROOT / "cfd" / "cases")
    ap.add_argument("--outdir", type=Path, default=RESULTS_DIR / "cfd")
    args = ap.parse_args(argv)
    found = {}
    for case in sorted(p for p in args.cases.iterdir() if p.is_dir()) if args.cases.is_dir() else []:
        res = case_residuals(case)
        if res is None or res["t"].size == 0:
            print(f"[skip] no solver log in {case}")
            continue
        found[case.name] = res
        path = write_residual_csv(res, args.outdir / f"residuals_{case.name}.csv")
        print(f"[residuals] {case.name}: {res['t'].size} time steps -> {path}")
    return found


if __name__ == "__main__":
    main()
