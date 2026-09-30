# piv_pressure/cfd/grid_independence.py
"""Grid-independence table for the cylinder wake (Figure 4).

For each mesh (coarse, medium, fine) the force-coefficient history written by
the ``forceCoeffs1`` function object is analysed over the last ``--window``
time units (after shedding is periodic):

* Strouhal number St = f D / U from the FFT of Cl(t) (Hann window, peak
  refined by parabolic interpolation; Cl is first resampled on a uniform
  time grid because the adaptive time step is not constant);
* mean drag coefficient and the Cl amplitude (max - min)/2;
* a periodicity check: Cl amplitude in the first vs second half of the window;
* the % change between successive meshes;
* Richardson extrapolation (observed order p, extrapolated value) and the
  Grid Convergence Index of the fine and medium meshes (Roache; procedure of
  Celik et al. 2008, J. Fluids Eng. 130, 078001 [verify]), with refinement
  ratios from the representative cell size h = sqrt(domain area / cells).

The header line of the .dat file is parsed, so both the ESI
(``coefficient.dat``: "# Time Cd Cd(f) Cd(r) Cl ...") and Foundation
(``forceCoeffs.dat``: "# Time Cm Cd Cl ...") layouts work [verify for your
version].

Usage:  python cfd/grid_independence.py [--window 80]
Output: results/cfd/grid_independence.csv / .md, results/cfd/cfd_summary.json,
        figure fig04_cfd_grid_independence.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt  # noqa: E402

from analysis.plots import PALETTE, RESULTS_DIR, save_figure, set_style  # noqa: E402

LEVELS = ("coarse", "medium", "fine")
#: 2-D fluid area of the domain: 30 D x 20 D minus the cylinder (pi D^2 / 4).
DOMAIN_AREA = 30.0 * 20.0 - np.pi / 4.0
SAFETY_FACTOR = 1.25
#: Published Re = 100 values for comparison.  All marked [verify]: check them
#: against the original papers before quoting; unbounded-domain values, whereas
#: this domain has +-10 D symmetry sides (5 % blockage), which raises Cd and St.
LITERATURE = [
    {"source": "Spec range (St 0.16-0.17) [verify]", "St": "0.16-0.17", "Cd_mean": "", "Cl_amp": ""},
    {"source": "Williamson 1996, experiments [verify]", "St": "~0.164", "Cd_mean": "", "Cl_amp": ""},
    {"source": "Park, Kwon & Choi 1998, 2-D DNS [verify]", "St": "0.165", "Cd_mean": "1.33",
     "Cl_amp": "0.33"},
]


def coefficient_files(case: Path, fo_name: str = "forceCoeffs1") -> list[Path]:
    """All force-coefficient files of a case (one per restart time)."""
    base = case / "postProcessing" / fo_name
    if not base.is_dir():
        return []
    return sorted([p for p in base.glob("*/*.dat") if p.name in ("coefficient.dat", "forceCoeffs.dat")],
                  key=lambda p: float(p.parent.name))


def read_coefficients(files: list[Path]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Concatenate (t, Cd, Cl) from the files, using the header to find columns."""
    t_all, cd_all, cl_all = [], [], []
    for f in files:
        header, rows = None, []
        for line in f.read_text().splitlines():
            if line.startswith("#"):
                tokens = line.lstrip("#").split()
                if tokens and tokens[0] == "Time":
                    header = tokens
                continue
            if line.strip():
                rows.append([float(v) for v in line.split()])
        if header is None or not rows:
            continue
        data = np.array(rows)
        t_all.append(data[:, header.index("Time")])
        cd_all.append(data[:, header.index("Cd")])
        cl_all.append(data[:, header.index("Cl")])
    if not t_all:
        raise FileNotFoundError("no readable force-coefficient data")
    t, cd, cl = np.concatenate(t_all), np.concatenate(cd_all), np.concatenate(cl_all)
    order = np.argsort(t, kind="stable")
    t, cd, cl = t[order], cd[order], cl[order]
    keep = np.concatenate([[True], np.diff(t) > 0])       # drop restart overlaps
    return t[keep], cd[keep], cl[keep]


def uniform_resample(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Linear resampling on a uniform grid with the median time step."""
    dt = float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, y)


def strouhal(t: np.ndarray, cl: np.ndarray, D: float = 1.0, U: float = 1.0) -> tuple[float, np.ndarray, np.ndarray]:
    """St from the dominant FFT peak of Cl; returns (St, frequencies, amplitude)."""
    tu, y = uniform_resample(t, cl)
    y = (y - y.mean()) * np.hanning(y.size)
    nfft = 8 * int(2 ** np.ceil(np.log2(y.size)))          # zero padding for a smooth spectrum
    spec = np.abs(np.fft.rfft(y, nfft))
    freq = np.fft.rfftfreq(nfft, tu[1] - tu[0])
    k = int(np.argmax(spec[1:]) + 1)
    if 1 <= k < spec.size - 1:                               # parabolic peak refinement
        a, b, c = np.log(spec[k - 1:k + 2] + 1e-300)
        k_ref = k + 0.5 * (a - c) / (a - 2 * b + c)
    else:
        k_ref = k
    f_peak = k_ref * (freq[1] - freq[0])
    return f_peak * D / U, freq, spec


def analyse(t: np.ndarray, cd: np.ndarray, cl: np.ndarray, window: float) -> dict:
    """St, mean Cd, Cl amplitude and a periodicity check over the last ``window``."""
    m = t >= t[-1] - window
    tw, cdw, clw = t[m], cd[m], cl[m]
    St, _, _ = strouhal(tw, clw)
    half = tw < tw[0] + 0.5 * (tw[-1] - tw[0])
    amp1 = 0.5 * (clw[half].max() - clw[half].min())
    amp2 = 0.5 * (clw[~half].max() - clw[~half].min())
    return {
        "St": St, "Cd_mean": float(np.mean(cdw)), "Cd_amp": float(0.5 * (cdw.max() - cdw.min())),
        "Cl_amp": float(0.5 * (clw.max() - clw.min())), "Cl_rms": float(np.std(clw)),
        "t_start": float(tw[0]), "t_end": float(tw[-1]),
        "periods_analysed": float((tw[-1] - tw[0]) * St),
        "amp_drift": float(abs(amp2 - amp1) / max(amp2, 1e-12)),
    }


def representative_h(cells: int, area: float = DOMAIN_AREA) -> float:
    """Representative 2-D cell size h = sqrt(area / cells)."""
    return float(np.sqrt(area / cells))


def richardson(f_fine: float, f_med: float, f_coarse: float, r21: float = 2.0,
               r32: float | None = None, max_iter: int = 100) -> dict:
    """Richardson extrapolation and GCI from three meshes (fine, medium, coarse).

    With e21 = f2 - f1 and e32 = f3 - f2 (1 = fine) the observed order solves
    (Celik et al. 2008)::

        p = | ln|e32/e21| + q(p) | / ln r21,   q(p) = ln((r21^p - s)/(r32^p - s)),
        s = sign(e32/e21)

    (for r21 = r32 = r this is simply p = ln(e32/e21)/ln r).  Then

        f_ext  = (r21^p f1 - f2) / (r21^p - 1)
        GCI_21 = 1.25 |e21/f1| / (r21^p - 1)      (fine mesh)
        GCI_32 = 1.25 |e32/f2| / (r32^p - 1)      (medium mesh)

    and GCI_32 / (r21^p GCI_21) close to 1 indicates the asymptotic range.
    Oscillatory convergence (s < 0) is flagged; the order is then unreliable.
    """
    r32 = r21 if r32 is None else r32
    e21, e32 = f_med - f_fine, f_coarse - f_med
    pct21 = 100.0 * abs(e21 / f_med) if f_med != 0 else float("nan")
    pct32 = 100.0 * abs(e32 / f_coarse) if f_coarse != 0 else float("nan")
    out = {"order": float("nan"), "extrapolated": float("nan"), "GCI_fine": float("nan"),
           "GCI_medium": float("nan"), "asymptotic_ratio": float("nan"),
           "pct_change_coarse_medium": pct32, "pct_change_medium_fine": pct21,
           "r21": r21, "r32": r32, "oscillatory": False, "monotone": False}
    if e21 == 0 or e32 == 0:
        return out
    s = float(np.sign(e32 / e21))
    out["oscillatory"] = s < 0
    p = np.log(abs(e32 / e21)) / np.log(r21)
    for _ in range(max_iter):                 # fixed-point iteration for unequal ratios
        q = np.log((r21**p - s) / (r32**p - s))
        p_new = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
        if abs(p_new - p) < 1e-12:
            break
        p = p_new
    if s < 0 or not np.isfinite(p) or p <= 0:
        return out
    ext = (r21**p * f_fine - f_med) / (r21**p - 1.0)
    gci21 = SAFETY_FACTOR * abs(e21 / f_fine) / (r21**p - 1.0)
    gci32 = SAFETY_FACTOR * abs(e32 / f_med) / (r32**p - 1.0)
    out.update(order=float(p), extrapolated=float(ext), GCI_fine=float(gci21), GCI_medium=float(gci32),
               asymptotic_ratio=float(gci32 / (r21**p * gci21)), monotone=True)
    return out


def table_markdown(rows: list[dict], rich: dict) -> str:
    """Markdown version of the grid-independence table."""
    out = ["| Mesh | Cells | St | mean Cd | Cl amplitude | periods analysed | amplitude drift |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        out.append(f"| {r['level']} | {r['cells']} | {r['St']:.4f} | {r['Cd_mean']:.4f} | "
                   f"{r['Cl_amp']:.4f} | {r['periods_analysed']:.1f} | {100 * r['amp_drift']:.2f} % |")
    if rich:
        out += ["", "| Quantity | % change coarse->medium | % change medium->fine | observed order p "
                "| Richardson extrapolated | GCI medium | GCI fine | asymptotic ratio |",
                "|---|---:|---:|---:|---:|---:|---:|---:|"]
        for q in ("St", "Cd_mean", "Cl_amp"):
            if q in rich:
                R = rich[q]
                out.append(f"| {q} | {R['pct_change_coarse_medium']:.2f} % | {R['pct_change_medium_fine']:.2f} % "
                           f"| {R['order']:.2f} | {R['extrapolated']:.4f} | {100 * R['GCI_medium']:.2f} % "
                           f"| {100 * R['GCI_fine']:.2f} % | {R['asymptotic_ratio']:.3f} |")
    out += ["", "| Literature (Re = 100) | St | mean Cd | Cl amplitude |", "|---|---:|---:|---:|"]
    for L in LITERATURE:
        out.append(f"| {L['source']} | {L['St']} | {L['Cd_mean']} | {L['Cl_amp']} |")
    return "\n".join(out) + "\n"


def plot(rows: list[dict], series: dict, rich: dict) -> None:
    """Figure 4: table (top), Cl(t) (bottom left) and Cl spectra (bottom right)."""
    fig = plt.figure(figsize=(9.0, 6.4))
    ax_t = fig.add_axes([0.03, 0.62, 0.94, 0.34])
    ax_t.axis("off")
    cells = [[r["level"], f"{r['cells']:d}", f"{r['St']:.4f}", f"{r['Cd_mean']:.4f}", f"{r['Cl_amp']:.4f}"]
             for r in rows]
    if all(q in rich for q in ("St", "Cd_mean", "Cl_amp")):
        cells.append(["Richardson extrap.", "", *(f"{rich[q]['extrapolated']:.4f}"
                                                   for q in ("St", "Cd_mean", "Cl_amp"))])
    for L in LITERATURE:
        cells.append([L["source"], "", L["St"], L["Cd_mean"], L["Cl_amp"]])
    tab = ax_t.table(cellText=cells, colLabels=["Mesh / source", "Cells", "St", "mean Cd", "Cl amplitude"],
                     loc="center", cellLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8)
    tab.scale(1.0, 1.25)
    ax_t.set_title("Grid independence, cylinder wake Re = 100 (literature values: [verify])")
    ax1 = fig.add_axes([0.08, 0.08, 0.40, 0.44])
    ax2 = fig.add_axes([0.57, 0.08, 0.40, 0.44])
    markers = {"coarse": "v", "medium": "s", "fine": "o"}
    for k, (level, (t, cl, f, spec)) in enumerate(series.items()):
        m = t >= t[-1] - 30.0
        ax1.plot(t[m], cl[m], color=PALETTE[k], lw=1.0, label=level, marker=markers[level],
                 markevery=max(1, int(m.sum() / 12)), ms=4)
        ax2.plot(f, spec / spec.max(), color=PALETTE[k], lw=1.2, label=level, marker=markers[level],
                 markevery=max(1, f.size // 40), ms=4)
    ax1.set_xlabel(r"time $t U/D$ (-)")
    ax1.set_ylabel(r"$C_L$ (-)")
    ax1.set_title("Lift coefficient, last 30 time units")
    ax1.legend(loc="upper right", ncol=3)
    ax2.set_xlim(0, 0.6)
    ax2.set_xlabel(r"frequency $f D/U$ (-)")
    ax2.set_ylabel("normalised amplitude (-)")
    ax2.set_title(r"FFT of $C_L$ (Hann window)")
    ax2.legend(loc="upper right")
    save_figure(fig, "fig04_cfd_grid_independence")


def main(argv: list[str] | None = None) -> list[dict]:
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", type=Path, default=ROOT / "cfd" / "cases")
    ap.add_argument("--window", type=float, default=80.0, help="analysis window (time units)")
    ap.add_argument("--outdir", type=Path, default=RESULTS_DIR / "cfd")
    args = ap.parse_args(argv)
    set_style()
    rows, series = [], {}
    for level in LEVELS:
        case = args.cases / level
        files = coefficient_files(case)
        if not files:
            print(f"[skip] no forceCoeffs output in {case}")
            continue
        t, cd, cl = read_coefficients(files)
        res = analyse(t, cd, cl, args.window)
        info = json.loads((case / "mesh_info.json").read_text()) if (case / "mesh_info.json").is_file() else {}
        res.update(level=level, cells=int(info.get("cells", 0)))
        rows.append(res)
        _, f, spec = strouhal(t[t >= t[-1] - args.window], cl[t >= t[-1] - args.window])
        series[level] = (t, cl, f, spec)
        if res["amp_drift"] > 0.02:
            print(f"[warn] {level}: Cl amplitude changes by {100 * res['amp_drift']:.1f} % within "
                  "the window - shedding may not be fully periodic yet")
    if not rows:
        print("No CFD results found - run cfd/run_cases.sh first.")
        return rows
    rich = {}
    for r in rows:
        r["h_rep"] = representative_h(r["cells"]) if r["cells"] > 0 else float("nan")
    if [r["level"] for r in rows] == list(LEVELS):
        h1, h2, h3 = rows[2]["h_rep"], rows[1]["h_rep"], rows[0]["h_rep"]
        r21, r32 = h2 / h1, h3 / h2
        if not (np.isfinite(r21) and np.isfinite(r32) and r21 > 1.01 and r32 > 1.01):
            print("[warn] cell counts missing or not refined - assuming refinement ratio 2")
            r21 = r32 = 2.0
        for q in ("St", "Cd_mean", "Cl_amp"):
            rich[q] = richardson(rows[2][q], rows[1][q], rows[0][q], r21, r32)
        args.outdir.mkdir(parents=True, exist_ok=True)
        with (args.outdir / "richardson.csv").open("w") as fh:
            keys = ["order", "extrapolated", "GCI_fine", "GCI_medium", "asymptotic_ratio",
                    "pct_change_coarse_medium", "pct_change_medium_fine", "r21", "r32", "monotone",
                    "oscillatory"]
            fh.write("quantity," + ",".join(keys) + "\n")
            for q, R in rich.items():
                fh.write(q + "," + ",".join(f"{R[k]:.6g}" if isinstance(R[k], float) else str(R[k])
                                             for k in keys) + "\n")
    args.outdir.mkdir(parents=True, exist_ok=True)
    keys = ["level", "cells", "h_rep", "St", "Cd_mean", "Cd_amp", "Cl_amp", "Cl_rms", "t_start",
            "t_end", "periods_analysed", "amp_drift"]
    with (args.outdir / "grid_independence.csv").open("w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    (args.outdir / "grid_independence.md").write_text(table_markdown(rows, rich))
    finest = rows[-1]
    (args.outdir / "cfd_summary.json").write_text(json.dumps(
        {"St": finest["St"], "period": 1.0 / finest["St"], "level": finest["level"],
         "richardson": rich}, indent=2))
    plot(rows, series, rich)
    print(table_markdown(rows, rich))
    return rows


if __name__ == "__main__":
    main()
