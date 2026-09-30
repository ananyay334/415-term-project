# piv_pressure/run_all.py
"""Run every study, regenerate every figure and write results/SUMMARY.md.

Order matters only for the combined figures: run_bc before run_method
(figure 8), run_window before run_denoise (figure 9) and cfd_grid before
cfd_convergence (figure 18).  CFD parts - including the optional steady
Re = 20 case - are skipped automatically if the OpenFOAM results are missing.

The convergence / error-analysis studies (Lecture 11) are ``residuals``
(iterative solvers, norms, SOR stability; figures 11-15), ``order``
(first- vs second-order grid convergence; figure 16) and
``cfd_convergence`` (solver residuals, periodicity, steady-state criterion,
Richardson/GCI; figures 17-18).

Usage:
    python run_all.py              # full run (10 seeds, grids up to 256^2)
    python run_all.py --quick      # fast check (3 seeds, coarser sweeps)
    python run_all.py --only noise dt
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from analysis import summary  # noqa: E402
from cfd import grid_independence  # noqa: E402
from studies import (run_bc, run_cfd_convergence, run_cfd_validation, run_convergence,  # noqa: E402
                     run_denoise, run_dt, run_method, run_noise, run_order, run_residuals,
                     run_window)

STUDIES = {
    "convergence": run_convergence.main,
    "noise": run_noise.main,
    "dt": run_dt.main,
    "bc": run_bc.main,
    "method": run_method.main,
    "window": run_window.main,
    "denoise": run_denoise.main,
    "cfd_grid": grid_independence.main,
    "cfd_validation": run_cfd_validation.main,
    "residuals": run_residuals.main,
    "order": run_order.main,
    "cfd_convergence": run_cfd_convergence.main,
}


def main(argv: list[str] | None = None) -> None:
    """Run the selected studies and the summary."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="fewer seeds and coarser sweeps")
    ap.add_argument("--seeds", type=int, default=10, help="noise realisations (full run)")
    ap.add_argument("--only", nargs="+", choices=list(STUDIES), help="run only these studies")
    ap.add_argument("--no-animation", action="store_true", help="skip the CFD GIF")
    args = ap.parse_args(argv)
    selected = args.only or list(STUDIES)
    t_all = time.perf_counter()
    for name in selected:
        t0 = time.perf_counter()
        print(f"\n===== {name} =====")
        if name == "cfd_grid":
            argv_s: list[str] = []
        else:
            argv_s = ["--seeds", str(args.seeds)] + (["--quick"] if args.quick else [])
            if name == "cfd_validation" and not args.no_animation:
                argv_s.append("--animate")
        STUDIES[name](argv_s)
        print(f"----- {name}: {time.perf_counter() - t0:.1f} s")
    summary.main()
    print(f"\nAll done in {time.perf_counter() - t_all:.1f} s. Figures: results/figures, "
          "tables: results/*.csv, summary: results/SUMMARY.md")


if __name__ == "__main__":
    main()
