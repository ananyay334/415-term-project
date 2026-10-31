"""Extra presentation figures: equations, hero field, split bar charts, CFD domain schematic."""
import sys, os, csv, json
from pathlib import Path
ROOT = Path(sys.argv[1]); OUT = Path(sys.argv[2])
sys.path.insert(0, str(ROOT)); os.environ["PIV_RESULTS_DIR"] = str(OUT / "_tmp")
exec(open(Path(__file__).with_name("figs_style.py")).read())

# ---------------- equations (transparent PNG)
EQ = {
    "eq_G": r"$\nabla p \;=\; \mathbf{G} \;=\; -\rho\,\dfrac{D\mathbf{u}}{Dt} \;+\; \mu\,\nabla^{2}\mathbf{u}$",
    "eq_euler": r"$\left(\dfrac{D\mathbf{u}}{Dt}\right)_{E} \approx \dfrac{\mathbf{u}^{n+1}-\mathbf{u}^{n-1}}{2\,\Delta t} + (\mathbf{u}^{n}\!\cdot\!\nabla)\,\mathbf{u}^{n}$",
    "eq_lagr": r"$\left(\dfrac{D\mathbf{u}}{Dt}\right)_{L} \approx \dfrac{\mathbf{u}(\mathbf{x}+\mathbf{u}\Delta t,\,t^{n+1})-\mathbf{u}(\mathbf{x}-\mathbf{u}\Delta t,\,t^{n-1})}{2\,\Delta t}$",
    "eq_poisson": r"$\nabla^{2} p = \nabla\!\cdot\!\mathbf{G}\;\;\mathrm{in}\;\Omega, \qquad \dfrac{\partial p}{\partial n} = \mathbf{n}\!\cdot\!\mathbf{G}\;\;\mathrm{on}\;\partial\Omega$",
    "eq_compat": r"$\oint_{\partial\Omega} \mathbf{n}\!\cdot\!\mathbf{G}\;ds \;=\; \iint_{\Omega} \nabla\!\cdot\!\mathbf{G}\;dA$",
    "eq_eps": r"$\varepsilon_{L2} = \dfrac{\| (p-\bar p) - (p_{ref}-\bar p_{ref}) \|_{2}}{\| p_{ref}-\bar p_{ref} \|_{2}}$",
    "eq_direct": r"$p(x,y) = p_0 + \int_C \mathbf{G}\cdot d\mathbf{l}$",
    "eq_lap": r"$\dfrac{p_{i+1,j}+p_{i-1,j}+p_{i,j+1}+p_{i,j-1}-4p_{i,j}}{h^{2}} = (\nabla\!\cdot\!\mathbf{G})_{i,j}$",
    "eq_ns": r"$\nabla\!\cdot\!\mathbf{u}=0,\qquad \rho\left(\dfrac{\partial \mathbf{u}}{\partial t}+(\mathbf{u}\!\cdot\!\nabla)\mathbf{u}\right) = -\nabla p + \mu\nabla^{2}\mathbf{u}$",
}
sizes = {}
for k, tex in EQ.items():
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, tex, fontsize=22, color="#111111")
    fig.savefig(OUT / f"{k}.png", dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

# ---------------- split bars
bc = rows("bc"); meth = rows("method")
BCL = [("neumann", "All-Neumann, pin mean"), ("neumann_node", "All-Neumann, pin node"),
       ("exact_1", "Exact p on 1 edge"), ("exact_2", "Exact p on 2 edges"),
       ("bernoulli_1", "Bernoulli on 1 edge"), ("bernoulli_2", "Bernoulli on 2 edges")]
MTL = [("direct_path1", "Direct, path 1"), ("direct_path2", "Direct, path 2"), ("direct", "Direct, path avg."),
       ("poisson_velocity", "Poisson, velocity source"), ("poisson", "Poisson, div G source"), ("omni", "Omni-directional")]
def bars(data, labels, name, legend):
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    y = np.arange(len(labels)); w = 0.38
    for j, (case, col) in enumerate((("tg", PAL[0]), ("kov20", PAL[1]))):
        mu = [stat(data, "eps_l2", case=case, condition="noisy", variant=k)[0] for k, _ in labels]
        sd = [stat(data, "eps_l2", case=case, condition="noisy", variant=k)[1] for k, _ in labels]
        ax.barh(y + (j - 0.5) * w, mu, w, xerr=sd, color=col, alpha=0.9, label=CASES[case],
                error_kw=dict(ecolor=INK2, lw=1, capsize=2))
        for yy, m, e in zip(y + (j - 0.5) * w, mu, sd):
            ax.text(m + e + 0.015, yy, f"{m:.3f}", va="center", fontsize=9, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in labels]); ax.invert_yaxis()
    ax.set_xlabel(r"error $\varepsilon_{L2}$  (σ = 1 %, mean ± std of 10 seeds)")
    ax.grid(axis="y", visible=False); ax.set_xlim(0, 1.05)
    if legend: ax.legend(loc=legend if isinstance(legend, str) else "lower right")
    save(fig, name)
bars(meth, MTL, "f_method_bar", True)
bars(bc, BCL, "f_bc_bar", "upper right")

# ---------------- hero: reconstructed TG core pressure
from studies.common import make_case
from core.pipeline import reconstruct
fld = make_case("tg_core", n=64, window_h=1.0, sigma=0.01, seed=0)
rec = reconstruct(fld, method="poisson")
p = rec.p - rec.p.mean(); X, Y = np.meshgrid(fld.x, fld.y); vmax = np.abs(rec.p_ref - rec.p_ref.mean()).max()
fig, ax = plt.subplots(figsize=(4.4, 3.9))
im = ax.pcolormesh(X, Y, p, cmap=DIV, vmin=-vmax, vmax=vmax, shading="auto")
n = fld.u.shape[0] // 2; s = 5
ax.quiver(X[::s, ::s], Y[::s, ::s], fld.u[n][::s, ::s], fld.v[n][::s, ::s], color="#222", alpha=0.6, width=0.005)
ax.set_aspect("equal"); ax.grid(False); ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=ax, shrink=0.85, label=r"$p-\bar p$ ($\rho U^2$)")
ax.set_title("Pressure rebuilt from noisy\nPIV-like velocity (arrows)", fontsize=12)
save(fig, "f_hero")

# ---------------- CFD domain schematic
from matplotlib.patches import Rectangle, Circle
fig, ax = plt.subplots(figsize=(7.2, 3.6)); ax.set_xlim(-11.5, 21.5); ax.set_ylim(-13.6, 11.5); ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(Rectangle((-10, -10), 30, 20, fc="#f4f8fb", ec=INK, lw=1.4))
ax.add_patch(Circle((0, 0), 0.5, fc="#777", ec=INK, lw=1, zorder=4))
ax.add_patch(Rectangle((1, -2), 4, 4, fc="none", ec=PAL[1], lw=2, zorder=5))
ax.text(5.3, 2.4, "PIV window\n81 × 81 points\nx/D ∈ [1, 5], y/D ∈ [−2, 2]", color=PAL[1], fontsize=9.5, va="bottom")
for k in range(5):
    ax.annotate("", xy=(-8.0, -7 + 3.5 * k), xytext=(-10.0, -7 + 3.5 * k), arrowprops=dict(arrowstyle="->", color=PAL[0], lw=1.4))
ax.text(-10.3, 0, "inlet\nU = (1, 0)\n∂p/∂n = 0", ha="right", va="center", fontsize=9.5, color=PAL[0])
ax.text(20.3, 0, "outlet\np = 0\n∂U/∂n = 0", ha="left", va="center", fontsize=9.5)
ax.text(5, 10.4, "symmetry plane (±10 D)", ha="center", fontsize=9.5); ax.text(5, -12.9, "symmetry plane", ha="center", fontsize=9.5)
ax.text(1.0, -3.2, "cylinder D = 1, no-slip", fontsize=9.5, ha="center", va="top", bbox=dict(fc="#f4f8fb", ec="none", pad=1))
ax.annotate("", xy=(0, -10.6), xytext=(-10, -10.6), arrowprops=dict(arrowstyle="<->", color=INK2)); ax.text(-5, -9.8, "10 D", ha="center", fontsize=9, color=INK2)
ax.annotate("", xy=(20, -10.6), xytext=(0, -10.6), arrowprops=dict(arrowstyle="<->", color=INK2)); ax.text(10, -9.8, "20 D", ha="center", fontsize=9, color=INK2)
for sgn, c in ((1, PAL[0]), (-1, PAL[7])):
    for xc in np.arange(7.5, 19, 5.0) + (0 if sgn > 0 else 2.5):
        t = np.linspace(0, 2 * np.pi, 50); ax.plot(xc + 0.9 * np.cos(t), sgn * 0.9 + 0.7 * np.sin(t), color=c, lw=1, alpha=0.6)
save(fig, "f_cfd_domain")

# sizes for layout
for f in sorted(OUT.glob("*.png")):
    from PIL import Image
    w, h = Image.open(f).size; sizes[f.stem] = [w, h]
json.dump(sizes, open(OUT / "sizes.json", "w"), indent=1)
print(len(sizes), "images")
