import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle, Polygon, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap

PAL = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
plt.rcParams.update({
    "font.family": "DejaVu Serif", "mathtext.fontset": "dejavuserif",
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12, "legend.fontsize": 10.5,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 250, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    "legend.frameon": False,
})
DIV = LinearSegmentedColormap.from_list("div", ["#1d4f91", "#5b9bd5", "#f2f2f2", "#e8826b", "#a3241a"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#ffffff", "#9ec5ea", "#2a78d6", "#0b2f66"])
M = {"poisson": (PAL[0], "o", "Poisson"), "direct": (PAL[1], "s", "Direct (path avg.)"),
     "omni": (PAL[2], "^", "Omni-directional")}
CASES = {"tg": "Taylor–Green", "kov20": "Kovasznay Re = 20", "kov40": "Kovasznay Re = 40"}

def rows(name):
    with open(ROOT / "results" / f"{name}.csv") as f:
        return list(csv.DictReader(f))

def stat(rs, val, **m):
    v = [float(r[val]) for r in rs if all(str(r.get(k)) == str(x) for k, x in m.items()) and r[val] not in ("", "nan")]
    v = np.array(v)
    return (v.mean(), v.std(ddof=1) if v.size > 1 else 0.0) if v.size else (np.nan, np.nan)

def save(fig, name):
    fig.savefig(OUT / f"{name}.png"); plt.close(fig); print("saved", name)

