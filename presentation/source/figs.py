"""Presentation figures for the PIV-pressure deck (generated from the project's code and CSVs)."""
import sys, os, csv, math
from pathlib import Path
ROOT = Path(sys.argv[1])            # piv_pressure dir
OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT))
os.environ["PIV_RESULTS_DIR"] = str(OUT / "_tmp")

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

# ---------------------------------------------------------------- 1. grid convergence (clean / noisy)
conv = rows("convergence")
def conv_fig(condition, name):
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.7), sharey=(condition == "noisy"))
    for ax, case in zip(axs, CASES):
        best = None; slopes = []
        for meth, (c, mk, lab) in M.items():
            hs = sorted({float(r["h"]) for r in conv if r["case"] == case and r["method"] == meth and r["condition"] == condition})
            mu, sd = zip(*[stat(conv, "eps_l2", case=case, method=meth, condition=condition, h=r_h) for r_h in
                           [next(r["h"] for r in conv if r["case"] == case and float(r["h"]) == h) for h in hs]]) if hs else ([], [])
            mu, sd = np.array(mu), np.array(sd)
            ok = np.isfinite(mu)
            hs = np.array(hs)[ok]; mu = mu[ok]; sd = sd[ok]
            if condition == "clean":
                p = np.polyfit(np.log(hs), np.log(mu), 1)[0]
                slopes.append(f"{lab.split()[0]} {p:.2f}")
            ax.plot(hs, mu, color=c, marker=mk, lw=2, ms=6, label=lab)
            if condition == "noisy" and np.any(sd > 0):
                ax.fill_between(hs, np.maximum(mu - sd, 1e-4), mu + sd, color=c, alpha=0.15, lw=0)
            if condition == "noisy" and meth == "poisson":
                i = int(np.argmin(mu)); best = (hs[i], mu[i])
        if condition == "clean":
            h0 = np.array([0.008, 0.05]); ref = 0.5 * (mu[-1] if False else 1)
            yref = 0.6 * ax.get_lines()[0].get_ydata()[0] * (h0 / 0.15) ** 2
            ax.plot(h0, yref, "--", color=INK2, lw=1.2, label="reference slope 2")
            ax.text(0.97, 0.04, "fitted slopes:\n" + "\n".join(slopes), transform=ax.transAxes, va="bottom", ha="right", fontsize=9.5, color=INK2)
        else:
            ax.plot(*best, marker="*", ms=17, color=INK, ls="none", label="Poisson optimum", zorder=5)
            ax.annotate(f"h = {best[0]:.3f}", xy=best, xytext=(best[0] * 0.75, best[1] * 0.55), fontsize=10, ha="right")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(CASES[case]); ax.set_xlabel("PIV grid spacing h")
    axs[0].set_ylabel(r"error $\varepsilon_{L2}$")
    hd, lb = axs[0].get_legend_handles_labels()
    fig.legend(hd, lb, loc="lower center", ncol=len(lb), bbox_to_anchor=(0.5, -0.14), fontsize=10.5)
    save(fig, name)

conv_fig("clean", "f_conv_clean")
conv_fig("noisy", "f_conv_noisy")

# ---------------------------------------------------------------- 2. noise (log-log)
noise = rows("noise")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
for ax, case in zip(axs, ("tg", "kov20")):
    for meth, (c, mk, lab) in M.items():
        sig = sorted({float(r["sigma"]) for r in noise if r["case"] == case and float(r["sigma"]) > 0})
        mu = np.array([stat([r for r in noise if float(r["sigma"]) == s], "eps_l2", case=case, method=meth)[0] for s in sig])
        sd = np.array([stat([r for r in noise if float(r["sigma"]) == s], "eps_l2", case=case, method=meth)[1] for s in sig])
        s = np.array(sig) * 100
        slope = np.polyfit(np.log(s), np.log(mu), 1)[0]
        ax.plot(s, mu, color=c, marker=mk, lw=2, ms=6, label=f"{lab} (slope {slope:.2f})")
        ax.fill_between(s, mu - sd, mu + sd, color=c, alpha=0.15, lw=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(CASES[case] + r"  (n = 64, W = h)")
    ax.set_xlabel(r"velocity noise $\sigma$ (% of $U_{ref}$)")
    ax.legend(loc="upper left", fontsize=9.5)
axs[0].set_ylabel(r"error $\varepsilon_{L2}$")
save(fig, "f_noise")

# ---------------------------------------------------------------- 3. BC + method bars (noisy)
bc = rows("bc"); meth = rows("method")
BCL = [("neumann", "All-Neumann, pin mean"), ("neumann_node", "All-Neumann, pin node"),
       ("exact_1", "Exact p on 1 edge"), ("exact_2", "Exact p on 2 edges"),
       ("bernoulli_1", "Bernoulli on 1 edge"), ("bernoulli_2", "Bernoulli on 2 edges")]
MTL = [("direct_path1", "Direct, path 1"), ("direct_path2", "Direct, path 2"), ("direct", "Direct, path avg."),
       ("poisson_velocity", "Poisson, velocity source"), ("poisson", "Poisson, div G source"), ("omni", "Omni-directional")]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
for ax, data, labels, title in ((axs[0], meth, MTL, "(a) Integration method"), (axs[1], bc, BCL, "(b) Boundary condition (Poisson)")):
    y = np.arange(len(labels)); w = 0.38
    for j, (case, col) in enumerate((("tg", PAL[0]), ("kov20", PAL[1]))):
        mu = [stat(data, "eps_l2", case=case, condition="noisy", variant=k)[0] for k, _ in labels]
        sd = [stat(data, "eps_l2", case=case, condition="noisy", variant=k)[1] for k, _ in labels]
        ax.barh(y + (j - 0.5) * w, mu, w, xerr=sd, color=col, alpha=0.9, label=CASES[case],
                error_kw=dict(ecolor=INK2, lw=1, capsize=2))
        for yy, m, e in zip(y + (j - 0.5) * w, mu, sd):
            ax.text(m + e + 0.015, yy, f"{m:.3f}", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in labels]); ax.invert_yaxis()
    ax.set_xlabel(r"error $\varepsilon_{L2}$ (σ = 1 %, mean ± std, 10 seeds)")
    ax.set_title(title, loc="left"); ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 1.05)
axs[0].legend(loc="lower right", fontsize=9.5)
fig.tight_layout(w_pad=2.5)
save(fig, "f_bc_method")

# ---------------------------------------------------------------- 4. time step
dt = rows("dt")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
for ax, (case, cond, title) in zip(axs, (("tgc", "clean", "Convected vortex, clean data"),
                                          ("tg", "noisy", "Stationary vortex, σ = 1 % noise"))):
    for sch, c, mk, lab in (("eulerian", PAL[0], "o", "Eulerian"), ("lagrangian", PAL[1], "s", "Lagrangian")):
        d = sorted({float(r["dt"]) for r in dt if r["case"] == case})
        mu = np.array([stat([r for r in dt if float(r["dt"]) == x], "eps_l2", case=case, condition=cond, scheme=sch)[0] for x in d])
        sd = np.array([stat([r for r in dt if float(r["dt"]) == x], "eps_l2", case=case, condition=cond, scheme=sch)[1] for x in d])
        ax.plot(d, mu, color=c, marker=mk, lw=2, ms=6, label=lab)
        if cond == "noisy": ax.fill_between(d, mu - sd, mu + sd, color=c, alpha=0.15, lw=0)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(title)
    ax.set_xlabel(r"time between snapshots $\Delta t$"); ax.legend(loc="upper left" if cond == "clean" else "upper right")
axs[0].set_ylabel(r"error $\varepsilon_{L2}$")
axs[0].annotate("truncation error\n" + r"$\propto \Delta t^2$ grows", xy=(0.2, 0.018), xytext=(0.035, 0.008),
                fontsize=9.5, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2))
axs[1].annotate("noise amplified\n" + r"as $\sigma/\Delta t$", xy=(0.012, 0.25), xytext=(0.03, 0.2),
                fontsize=9.5, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2))
save(fig, "f_dt")

# ---------------------------------------------------------------- 5. denoising
den = rows("denoise")
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
for k, c, mk, lab in (("none", PAL[0], "o", "No denoising"), ("gaussian", PAL[1], "s", "Gaussian filter"),
                      ("savgol", PAL[2], "^", "Savitzky–Golay")):
    sig = sorted({float(r["sigma"]) for r in den if r["case"] == "tg_core"})
    mu = np.array([stat([r for r in den if float(r["sigma"]) == s], "core_ratio", case="tg_core", denoise=k)[0] for s in sig])
    e = np.array([stat([r for r in den if float(r["sigma"]) == s], "eps_l2", case="tg_core", denoise=k)[0] for s in sig])
    s = np.array(sig) * 100
    axs[0].plot(s, mu, color=c, marker=mk, lw=2, ms=6, label=lab)
    axs[1].plot(s[1:], e[1:], color=c, marker=mk, lw=2, ms=6, label=lab)
    axs[1].plot([0.2], [e[0]], marker=mk, color=c, ms=7, mfc="white", mew=1.6)
axs[0].axhline(1, ls="--", color=INK2, lw=1.2)
axs[0].set_title("Vortex-core suction (reconstructed / exact)"); axs[0].set_xlabel(r"noise $\sigma$ (% of $U_{ref}$)")
axs[0].set_ylabel(r"$p_{core,rec}\,/\,p_{core,ref}$"); axs[0].legend(loc="upper left")
axs[1].set_xscale("log"); axs[1].set_yscale("log"); axs[1].set_title(r"Field error $\varepsilon_{L2}$")
axs[1].set_xlabel(r"noise $\sigma$ (% of $U_{ref}$); hollow = clean data")
axs[1].set_xticks([0.2, 0.5, 1, 2, 5, 10]); axs[1].set_xticklabels(["0", "0.5", "1", "2", "5", "10"])
save(fig, "f_denoise")

# ---------------------------------------------------------------- 6. field maps (run the code)
from studies.common import make_case
from core.pipeline import reconstruct
def fields_fig(case, name, title):
    fld = make_case(case, n=64, window_h=1.0, sigma=0.01, seed=0)
    rec = reconstruct(fld, method="poisson")
    p = rec.p - rec.p.mean(); pr = rec.p_ref - rec.p_ref.mean()
    err = np.abs(p - pr) / (0.5 * fld.rho * fld.U_ref ** 2)
    X, Y = np.meshgrid(fld.x, fld.y)
    vmax = np.abs(pr).max()
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    for ax, f, t in ((axs[0], pr, "exact pressure"), (axs[1], p, "reconstructed (σ = 1 %, Poisson)")):
        im = ax.pcolormesh(X, Y, f, cmap=DIV, vmin=-vmax, vmax=vmax, shading="auto")
        n = fld.u.shape[0] // 2; s = 6
        ax.quiver(X[::s, ::s], Y[::s, ::s], fld.u[n][::s, ::s], fld.v[n][::s, ::s], color="#222", alpha=0.55, width=0.004)
        ax.set_title(t); ax.set_aspect("equal"); ax.grid(False)
    fig.colorbar(im, ax=axs[:2], shrink=0.9, label=r"$p-\bar p$  ($\rho U^2$)")
    im2 = axs[2].pcolormesh(X, Y, err, cmap=SEQ, shading="auto")
    axs[2].set_title(r"error $|p-p_{ref}|\,/\,\frac{1}{2}\rho U^2$"); axs[2].set_aspect("equal"); axs[2].grid(False)
    fig.colorbar(im2, ax=axs[2], shrink=0.9)
    for ax in axs: ax.set_xlabel("x")
    axs[0].set_ylabel("y")
    eps = np.linalg.norm(p - pr) / np.linalg.norm(pr)
    print(name, "eps_l2 =", eps)
    save(fig, name)
    return eps
eps_tg = fields_fig("tg_core", "f_fields_tg", "Taylor–Green")
eps_kv = fields_fig("kov20", "f_fields_kov", "Kovasznay")

# ---------------------------------------------------------------- 7. test-case thumbnails
from fields.analytic import taylor_green_uvp, kovasznay_uvp
fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
x = np.linspace(0.3, 2.5, 120); X, Y = np.meshgrid(x, x)
u, v, p = taylor_green_uvp(X, Y, 1.0)
axs[0].contourf(X, Y, p, 30, cmap=DIV); axs[0].streamplot(x, x, u, v, color="#222", density=1.0, linewidth=0.7, arrowsize=0.8)
axs[0].set_title("Taylor–Green vortex (unsteady)")
xk = np.linspace(-0.5, 1.5, 160); Xk, Yk = np.meshgrid(xk, xk)
uk, vk, pk = kovasznay_uvp(Xk, Yk, 20.0)
axs[1].contourf(Xk, Yk, pk, 30, cmap=DIV); axs[1].streamplot(xk, xk, uk, vk, color="#222", density=1.1, linewidth=0.7, arrowsize=0.8)
axs[1].set_title("Kovasznay flow, Re = 20 (steady)")
for ax in axs: ax.set_aspect("equal"); ax.grid(False); ax.set_xlabel("x"); ax.set_ylabel("y")
save(fig, "f_testcases")

# ---------------------------------------------------------------- 8. PIV schematic (intro)
rng = np.random.default_rng(3)
fig, ax = plt.subplots(figsize=(6.4, 4.0)); ax.set_xlim(0, 10); ax.set_ylim(-0.7, 6.3); ax.axis("off")
ax.add_patch(Polygon([[0.3, 0.5], [9.6, 0.5], [9.6, 3.6], [0.3, 3.6]], closed=True, fc="#c9f0d6", ec="#1baf7a", lw=1.2, alpha=0.6))
ax.text(9.45, 0.62, "laser light sheet", color="#0d7a52", fontsize=10.5, ha="right", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.add_patch(Circle((2.6, 2.05), 0.62, fc="#7c7c7c", ec=INK, lw=1.2, zorder=3))
ax.text(2.6, 2.05, "body", ha="center", va="center", color="white", fontsize=10, zorder=4)
for k in range(5):
    ax.annotate("", xy=(1.3, 0.9 + 0.62 * k), xytext=(0.35, 0.9 + 0.62 * k), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.2))
ax.text(0.4, 3.75, r"$U_\infty$", fontsize=12, color=INK2)
px = rng.uniform(0.4, 9.5, 260); py = rng.uniform(0.55, 3.55, 260)
keep = (px - 2.6) ** 2 + (py - 2.05) ** 2 > 0.7 ** 2
ax.scatter(px[keep], py[keep], s=3.5, color="#155e3f", zorder=2)
for cx, cy, s in ((4.6, 2.6, 1), (6.2, 1.5, -1), (7.8, 2.6, 1)):
    t = np.linspace(0, 2 * np.pi, 60)
    ax.plot(cx + 0.55 * np.cos(t), cy + 0.42 * np.sin(t), color=PAL[0] if s > 0 else PAL[7], lw=1.3, alpha=0.9)
    ax.annotate("", xy=(cx + 0.55 * s * 0.0 + 0.01, cy + 0.42), xytext=(cx - 0.05 * s, cy + 0.42),
                arrowprops=dict(arrowstyle="->", color=PAL[0] if s > 0 else PAL[7], lw=1.3))
# interrogation grid
for i in range(6):
    ax.plot([5.2 + 0.4 * i] * 2, [0.9, 2.9], color="#444", lw=0.5, alpha=0.6)
for j in range(6):
    ax.plot([5.2, 7.2], [0.9 + 0.4 * j] * 2, color="#444", lw=0.5, alpha=0.6)
ax.add_patch(Rectangle((5.6, 1.3), 0.4, 0.4, fc="none", ec=PAL[1], lw=2))
ax.annotate("interrogation window W\n→ one velocity vector", xy=(6.0, 1.5), xytext=(7.0, -0.55), fontsize=9.5,
            color=PAL[1], arrowprops=dict(arrowstyle="->", color=PAL[1]))
# camera
ax.add_patch(FancyBboxPatch((4.3, 4.85), 1.6, 0.9, boxstyle="round,pad=0.05", fc="#333", ec=INK))
ax.add_patch(Rectangle((4.85, 4.45), 0.5, 0.42, fc="#555", ec=INK))
ax.text(6.1, 5.25, "camera: two frames, Δt apart", fontsize=10, color=INK, va="center")
ax.plot([5.1, 2.0], [4.45, 3.6], color=INK2, lw=0.8, ls=":"); ax.plot([5.1, 8.6], [4.45, 3.6], color=INK2, lw=0.8, ls=":")
ax.text(5.0, 4.0, "seeding particles", fontsize=9.5, color="#155e3f", ha="center", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
save(fig, "f_piv_schematic")

# ---------------------------------------------------------------- 9. methods schematic
fig, axs = plt.subplots(1, 3, figsize=(12, 3.9))
for ax in axs:
    ax.set_xlim(-0.5, 6.5); ax.set_ylim(-0.5, 6.5); ax.set_aspect("equal"); ax.axis("off")
    for i in range(7):
        ax.plot([i, i], [0, 6], color=GRID, lw=0.8); ax.plot([0, 6], [i, i], color=GRID, lw=0.8)
    gx, gy = np.meshgrid(range(7), range(7)); ax.scatter(gx, gy, s=6, color="#9a9a9a")
# (a) direct
a = axs[0]; a.set_title("(a) Direct integration", fontsize=13)
a.plot([0, 4], [0, 0], color=PAL[3], lw=3); a.plot([4, 4], [0, 4], color=PAL[3], lw=3, label="path 1: x then y")
a.plot([0, 0], [0, 4], color=PAL[4], lw=3); a.plot([0, 4], [4, 4], color=PAL[4], lw=3, label="path 2: y then x")
a.plot(0, 0, "s", color=INK, ms=9); a.text(0.15, -0.45, r"$p_0$ (reference)", fontsize=10)
a.plot(4, 4, "o", color=PAL[1], ms=10); a.text(4.2, 4.25, r"$p(x,y)$", fontsize=11)
a.text(0.1, 5.3, r"$p = p_0 + \int_C \mathbf{G}\cdot d\mathbf{l}$", fontsize=12)
a.legend(loc="lower right", fontsize=9, bbox_to_anchor=(1.02, 0.12))
# (b) omni
b = axs[1]; b.set_title("(b) Omni-directional", fontsize=13)
cx, cy = 3.6, 2.6
for th in np.linspace(0, 2 * np.pi, 16, endpoint=False):
    d = np.array([np.cos(th), np.sin(th)]); ts = []
    for k, lim in ((0, (0, 6)), (1, (0, 6))):
        if abs(d[k]) > 1e-9:
            for L in lim:
                t = (L - (cx, cy)[k]) / d[k]
                if t > 0: ts.append(t)
    tmax = min(t for t in ts if all(-1e-9 <= (cx, cy)[kk] + t * d[kk] <= 6 + 1e-9 for kk in (0, 1)))
    b.plot([cx, cx + tmax * d[0]], [cy, cy + tmax * d[1]], color=PAL[2], lw=1.3, alpha=0.85)
    b.plot(cx + tmax * d[0], cy + tmax * d[1], "o", color=PAL[2], ms=4)
b.plot(cx, cy, "o", color=PAL[1], ms=10, zorder=5)
b.text(0.05, 6.15, "average over rays from all boundary points", fontsize=9.5)
# (c) Poisson stencil
c = axs[2]; c.set_title("(c) Pressure-Poisson", fontsize=13)
for (i, j) in ((3, 3), (2, 3), (4, 3), (3, 2), (3, 4)):
    c.plot(i, j, "o", color=PAL[0], ms=11 if (i, j) == (3, 3) else 8)
c.plot([2, 4], [3, 3], color=PAL[0], lw=1.5); c.plot([3, 3], [2, 4], color=PAL[0], lw=1.5)
c.text(4.25, 1.55, "5-point\nLaplacian", fontsize=9.5, color=PAL[0])
c.plot(0, 4, "o", color=PAL[0], ms=8); c.plot(-0.5, 4, "o", mfc="white", mec=PAL[7], ms=8, mew=1.6)
c.plot([0, 1], [4, 4], color=PAL[0], lw=1.2); c.plot([-0.5, 0], [4, 4], color=PAL[7], lw=1.2, ls="--")
c.annotate("ghost node:\n" + r"$\partial p/\partial n = \mathbf{n}\cdot\mathbf{G}$", xy=(-0.5, 4.0), xytext=(0.1, 5.05),
           fontsize=9.5, color=PAL[7], arrowprops=dict(arrowstyle="->", color=PAL[7]))
c.text(0.6, 0.6, r"$\nabla^2 p = \nabla\cdot\mathbf{G}$", fontsize=13)
c.set_xlim(-1.0, 6.5)
fig.tight_layout()
save(fig, "f_methods")

# ---------------------------------------------------------------- 10. noise-vs-resolution cartoon (why optimum h)
h = np.logspace(-2.2, -0.8, 100)
trunc = 3.0 * h ** 2; noise_t = 4e-4 / h
fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.plot(h, trunc, color=PAL[3], lw=2, label=r"truncation $\propto h^2$")
ax.plot(h, noise_t, color=PAL[7], lw=2, label=r"noise $\propto \sigma/h$")
ax.plot(h, trunc + noise_t, color=INK, lw=2.4, label="total")
i = np.argmin(trunc + noise_t); ax.plot(h[i], (trunc + noise_t)[i], "*", ms=16, color=INK)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("grid spacing h"); ax.set_ylabel("error")
ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_xticks([], minor=True); ax.set_yticks([], minor=True)
ax.text(h[8], trunc[8] * 1.6, r"truncation $\propto h^2$", color=PAL[3], fontsize=10.5)
ax.text(h[-30], noise_t[-30] * 0.38, r"noise $\propto \sigma/h$", color=PAL[7], fontsize=10.5)
ax.text(h[i], (trunc + noise_t)[i] * 1.45, "total", ha="center", fontsize=10.5); ax.set_title("Why an optimum spacing exists", fontsize=12)
save(fig, "f_tradeoff")
print("TG eps", eps_tg, "KOV eps", eps_kv)
