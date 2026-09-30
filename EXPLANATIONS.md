# EXPLANATIONS — Pressure from PIV

This document explains **every file and every function** of the project, the mathematics behind
them, and the questions you are likely to be asked in a viva, with short answers.

* **Part A** — derivations: the Navier–Stokes pressure gradient, material derivative, Poisson
  problem, ghost-point Neumann stencil, singular system, 2-D source identity, direct and
  omni-directional integration, PIV error model, metrics, pipeline pseudo-code.
* **Part B** — places where the original specification was wrong or would not work, and what
  was done instead.
* **Part C** — file-by-file: what each function does, which equation it implements, and viva
  questions.
* **Part D** — non-obvious results found while running the studies, and why they happen.
* **Part E** — the convergence and error-analysis add-on (Lecture 11): norms, scaled
  residuals, iterative solvers and their pseudo-time view, SOR stability, the steady-state
  criterion, first- vs second-order grid convergence, CFD residuals and periodicity,
  Richardson/GCI. It lists its own spec corrections and viva questions.

Notation: h is the grid spacing, Δt the time between snapshots, ρ the density, μ = ρν the
dynamic viscosity, **u** = (u, v) the velocity, and **G** the pressure gradient estimated
from the velocity. Arrays are stored as `[n, j, i]` = (time, y, x), so axis −1 is x and axis −2 is y.

---

# Part A — Theory and derivations

## A1. Pressure gradient from the Navier–Stokes equations

Incompressible, constant-property Navier–Stokes in 2-D:

    ρ ( ∂u/∂t + (u·∇)u ) = −∇p + μ ∇²u ,      ∇·u = 0

Solving for the pressure gradient gives the quantity that is computed from PIV data:

    ∇p = G ≡ −ρ Du/Dt + μ ∇²u ,      Du/Dt = ∂u/∂t + (u·∇)u

Component form (what `core/material_derivative.py` computes):

    Gx = −ρ ( u_t + u u_x + v u_y ) + μ ( u_xx + u_yy )
    Gy = −ρ ( v_t + u v_x + v v_y ) + μ ( v_xx + v_yy )

* The **material acceleration** Du/Dt is the acceleration of a fluid particle. It is the
  dominant term at moderate and high Reynolds number.
* The **viscous term** μ∇²u matters at low Re. In Kovasznay flow at Re = 20 or 40,
  dropping it leaves an O(1) error (test `test_viscous_term_is_needed`).
* In the non-dimensional CFD set-up ρ = 1, so p is the *kinematic* pressure p/ρ, exactly as
  in OpenFOAM's incompressible solvers.

## A2. The material derivative

### A2.1 Eulerian estimate

At a fixed node, with snapshots n−1, n and n+1:

    ∂u/∂t ≈ (uⁿ⁺¹ − uⁿ⁻¹) / (2Δt)                         error −(Δt²/6) ∂³u/∂t³
    ∂u/∂x ≈ (u_{i+1} − u_{i−1}) / (2h)  (interior)         error −(h²/6)  ∂³u/∂x³
    ∂u/∂x ≈ (−3u₀ + 4u₁ − u₂) / (2h)   (left edge)         error +(h²/3)  ∂³u/∂x³

Derivation of the one-sided stencil: write u₁ and u₂ as Taylor series about x₀.
- u₁ = u₀ + h u′ + h²/2 u″ + h³/6 u‴
- u₂ = u₀ + 2h u′ + 2h² u″ + 4h³/3 u‴

Choose a, b, c in (a u₀ + b u₁ + c u₂)/h so that the u₀ and u″ terms vanish and the u′ term
equals 1. This gives a = −3/2, b = 2, c = −1/2.

**Noise amplification.** If every velocity sample carries independent noise of standard
deviation s = σU:
- the time difference has noise s/(√2 Δt);
- a spatial first derivative has noise s/(√2 h);
- the Laplacian has noise ≈ √20 s/h² (five-point stencil: 1 + 1 + 1 + 1 + 16 = 20).

So the viscous noise μ√20 s/h² grows fastest as h → 0 (see Part D).

### A2.2 Lagrangian (pseudo-tracking) estimate

Follow a straight path through the node with the local velocity **a** = **u**(x, tⁿ):

    x_f = x + a Δt ,   x_b = x − a Δt
    Du/Dt ≈ [ u(x_f, tⁿ⁺¹) − u(x_b, tⁿ⁻¹) ] / (2Δt)

**Accuracy.** Let L = a·∇ + ∂/∂t, with a frozen. Then u(x ± aΔt, t ± Δt) = e^{±ΔtL} u, so

    estimate = sinh(Δt L) u / Δt = L u + (Δt²/6) L³u + O(Δt⁴)

L u = u_t + a·∇u is *exactly* Du/Dt, so the scheme is second-order in Δt.

**Why it can beat the Eulerian scheme.** Consider a pattern w carried at speed c:
u(x, t) = c + w(x − ct). Then:
- the Eulerian error ∝ ∂³u/∂t³ = −(c·∇)³w, which scales like **c³**;
- the Lagrangian error ∝ L³u = (w·∇)³w, which scales like **w³**, where w is the velocity
  *relative to the pattern*.

In a wake, the vortex-induced velocity w is smaller than the convection speed c, so the
Lagrangian scheme wins. For the stationary Taylor–Green vortex (c = 0) the Eulerian time
error is ≈ 0, so there it loses (Part B, item 4).

**Implementation details.**
* Velocities at x_f and x_b come from `scipy.interpolate.RegularGridInterpolator`: cubic if
  available, otherwise linear. Interpolation error is divided by 2Δt:
  - cubic: O(h⁴/Δt);
  - linear: O(h²/Δt), which *grows* as Δt → 0.
* **Pitfall found and fixed.** Since SciPy 1.13 the cubic spline coefficients are computed by
  an *iterative* solver with atol = 1e-6. On SciPy 1.15.3 the interpolation error stalled
  around 1e-5 whatever the grid (checked while building the project), and the Lagrangian
  error then *grew* under refinement. The code passes `solver=spsolve` (direct), which
  restores 4th-order interpolation.
* If x_f or x_b leaves the window, the interpolator returns NaN and the Eulerian value is
  used at that node. The fraction of such nodes is reported as `fallback_fraction`.
* Steady fields have a single snapshot, so the scheme tracks through that one snapshot with a
  pseudo-step Δt = h / max|u|.

## A3. The pressure-Poisson problem

Take the divergence of ∇p = G:

    ∇²p = ∇·G  in Ω ,        ∂p/∂n = n·G  on ∂Ω  (Neumann, all four edges)

**Compatibility.** Integrate ∇²p over Ω and apply Gauss's theorem: ∬ ∇²p dA = ∮ ∂p/∂n ds.
A solution can only exist if

    R = ∬ f dA − ∮ n·G ds = 0 .

For an exact gradient field R = 0. For discretised or noisy G, R ≠ 0; the code reports it
(`compatibility_residual`, and its relative value).

**Least-squares interpretation (important for the viva).** Minimising

    J[p] = ∬ |∇p − G|² dA

gives, through the calculus of variations, ∇²p = ∇·G in Ω and ∂p/∂n = n·G on ∂Ω. The Poisson
problem is therefore the *projection* of the measured G onto the space of gradient fields: it
keeps the curl-free part of G and throws away the curl, which can only be error.
- This is why Poisson integration is path-independent and filters much of the noise.
- It is also why the source term and the Neumann data must come *from the same G* (Part D).

## A4. Discretisation: 5-point Laplacian and ghost-point Neumann conditions

Unknowns sit at the grid nodes, k = j·nx + i. In the interior:

    (p_{i+1,j} − 2p_{i,j} + p_{i−1,j})/dx² + (p_{i,j+1} − 2p_{i,j} + p_{i,j−1})/dy² = f_{i,j}

with truncation error −(dx²/12) p_xxxx − (dy²/12) p_yyyy, i.e. second order.

**Ghost-point elimination** on the left edge (i = 0). Add a ghost node at i = −1. The Neumann
condition ∂p/∂x = Gx (outward normal −x, so ∂p/∂n = −Gx = n·G), discretised centrally, is

    (p₁ − p₋₁)/(2dx) = Gx,₀   ⇒   p₋₁ = p₁ − 2 dx Gx,₀

Substituting into the Laplacian row gives

    (2p₁ − 2p₀)/dx² + [y-part] = f₀ + 2Gx,₀/dx

Likewise:
- right edge: p_N = p_{N−2} + 2dx Gx, so b −= 2Gx/dx;
- bottom edge: b += 2Gy/dy;
- top edge: b −= 2Gy/dy.

**Corners.** A corner node has two ghost neighbours, one in x and one in y. Each is eliminated
with *its own* edge condition (∂p/∂x = Gx and ∂p/∂y = Gy). Both conditions are therefore used
exactly and no averaging is needed (Part B, item 2).

**Order of accuracy.** The boundary rows have a local truncation error of O(h), since
(2p₁ − 2p₀)/h² − 2g/h = p″ + (h/3)p‴ + …. The global error is nevertheless O(h²), because:
- the O(h) error lives only on a set of O(h) measure;
- multiplying edge rows by ½ and corner rows by ¼ turns the scheme into a vertex-centred
  finite-volume scheme on half and quarter cells, whose flux errors are O(h²).

The tests confirm second order: slopes 1.8–2.2 in `tests/test_poisson_exact.py`, and 1.96–2.03
between the two finest grids for the analytic cases.

**Symmetric form.** Let w = trapezoidal weights (1 interior, ½ on edges, ¼ at corners) and
W = diag(w). Then W·A is symmetric, which the test `test_matrix_properties_and_compatibility`
checks. S = −W A is symmetric positive semi-definite, which is what conjugate gradients needs.

## A5. The singular all-Neumann system

A·**1** = 0, because adding a constant to p changes nothing. A is therefore singular: the
pressure is defined only up to a constant.

The left null vector is w: wᵀA = 0. So A p = b is solvable only if wᵀb = 0, and

    dx·dy·wᵀb = trapezoidal ∬ f dA − trapezoidal ∮ n·G ds

This is the **discrete compatibility residual**, and the code reports exactly this number.

Two fixes are implemented:

1. **`pin="node"`** — replace row k by p_k = 0.
   - The remaining equations become solvable.
   - The row that was thrown away was the only place the incompatibility R could "go", so R
     acts as a point source at the pinned node. This gives a logarithmic disturbance
     ~ (R/2π) ln r around that node.
   - That is why pin-node errors are larger (BC study: TG clean 5.8e-3 vs 7.2e-4 for
     pin-mean).
2. **`pin="mean"`** — solve the augmented system

        [ A   1 ] [ p ]   [ b ]
        [ 1ᵀ  0 ] [ λ ] = [ 0 ]

   - This means A p + λ**1** = b and Σp = 0.
   - Multiplying the first row by wᵀ gives λ = wᵀb / wᵀ**1**.
   - The data are therefore corrected by a **uniform** source −λ, which spreads the
     incompatibility evenly instead of dumping it at one node.
   - The matrix is non-singular because the constraint row is not orthogonal to the null
     vector **1**, and the extra column **1** is not in the range of A (wᵀ**1** ≠ 0).

The conjugate-gradient path uses S = −WA:
- pinned or Dirichlet unknowns are *eliminated* (symmetric reduction);
- for the mean option, the right-hand side is made compatible by subtracting the same uniform
  correction.

The test `test_direct_and_cg_agree` checks that both routes give the same answer.

**Dirichlet edges.** Rows of nodes on a Dirichlet edge are replaced by p_k = value. Such an
edge anchors the constant, so no pin is needed.

## A6. The 2-D incompressible source-term identity

Take the divergence of G:
- ∇·(∂u/∂t) = ∂(∇·u)/∂t = 0;
- ∇·∇²u = ∇²(∇·u) = 0.

What remains is the convective term. With the summation convention,

    ∇·((u·∇)u) = ∂_i (u_j ∂_j u_i) = (∂_i u_j)(∂_j u_i) + u_j ∂_j(∂_i u_i)

The last term vanishes when ∇·u = 0. Written out in 2-D:

    (∂_i u_j)(∂_j u_i) = u_x² + u_y v_x + v_x u_y + v_y² = u_x² + 2 u_y v_x + v_y²

so

    ∇²p = −ρ [ u_x² + 2 u_y v_x + v_y² ]

Using ∇·u = 0 again (v_y = −u_x) gives the equivalent form 2ρ(u_x v_y − u_y v_x).

For Taylor–Green: u_x = −sin x sin y F, v_y = +sin x sin y F, u_y = cos x cos y F,
v_x = −cos x cos y F. Then

    −ρ[2 sin²x sin²y − 2cos²x cos²y]F² = ρ(cos 2x + cos 2y)F² = ∇²p   ✓

**Caveats.**
1. The identity assumes ∇·u = 0 exactly. Noisy PIV velocity is not solenoidal.
2. The Neumann data still come from G.
3. The quadratic terms bias noisy data by about +σ²/h² (this uniform part is removed by the
   compatibility correction).

## A7. Direct spatial integration

p(B) − p(A) = ∫_A^B G·dl, integrated with the cumulative trapezoidal rule, which is second
order:
- **path 1**: along the first row in x, then up each column in y;
- **path 2**: along the first column in y, then along each row in x.

For an exact gradient field both paths give the same answer (Stokes: ∮ G·dl = ∬ curl G dA = 0).
Measured G has a curl, so the two paths differ. The difference measures the **path
dependence**. Each path also accumulates its errors like a random walk, which is why direct
integration is by far the most noise-sensitive method (e.g. ε_L2 0.42 vs 0.082 for Poisson
at σ = 1 %, TG, n = 64).

## A8. Omni-directional integration (Liu & Katz 2006)

The idea is to average line integrals of G over many straight paths coming from all
directions, so that path-dependent errors cancel. The implementation (`core/omni_integration.py`)
is a *rotating-ray* variant:

1. Directions θ_k = 2πk/N, k = 0…N−1 (N = 32 by default).
2. From every node P, trace the ray backwards along −e_k to the boundary point E_k(P).
   Integrate G·e_k from E_k to P with the trapezoidal rule, step h, bilinear interpolation of
   G (both second order).
3. **Boundary stage.** For the M boundary nodes, the entry points are on the boundary too.
   p(E) is interpolated linearly along the perimeter between the two neighbouring boundary
   nodes. Averaging over the inward directions gives p_i = mean_k[Σ_m A_im p_m + I_i]: a
   linear system (I − A)p_b = c. Its null space is the constants, and because curl G ≠ 0 it is
   not exactly compatible. It is solved with the augmented zero-mean system (as in A5), which
   is the direct equivalent of the iteration used in the original method.
4. **Interior stage.** p(P) = mean over all N directions of [p_b(E_k) + ∫ G·e_k ds].

Liu & Katz (2006) use a "virtual boundary" outside the domain with parallel rays and iterate
the boundary values. This implementation keeps the idea (many directions, boundary values
consistent with the rays) but starts the rays at each node. Compare the details with the paper
before quoting it as the exact original algorithm **[verify]**.

Cost is O(N × nodes × ray length). Measured wall-clock times on an Apple M2 were about 0.4 s
at 64² (`noise.csv`), ~3 s at 128² and ~13 s at 256² (`convergence.csv`), against about
15 ms, 0.16 s and 0.75 s for Poisson. That is why the grid study runs omni on noisy data only
up to 64².

## A9. Bernoulli Dirichlet edge

For **steady, irrotational** flow, p + ½ρ|u|² is constant, so on a free-stream edge

    p = p∞ + ½ρ (U∞² − |u|²)

For unsteady *irrotational* flow the correct form is p + ρ ∂φ/∂t + ½ρ|u|² = const. The
steady form therefore misses ρ∂φ/∂t, and in rotational flow it is simply wrong.
- Neither Taylor–Green nor Kovasznay has an irrotational edge (both carry vorticity
  everywhere), so a Bernoulli edge *must* increase their error. The BC study shows it does.
- For the cylinder wake the |y| = 2D edges are nearly irrotational but unsteady.

## A10. Synthetic PIV model and denoising

* **Box filter, width W.** A continuous box average has transfer function
  sin(kW/2)/(kW/2) ≈ 1 − (kW)²/24. It therefore biases the velocity by −(W²/24)∇²u and
  attenuates peaks. Pressure is quadratic in u, so the vortex-core suction is under-predicted.
  - For an even number of points the plain moving average is off-centre by half a cell. The
    code averages the two half-shifted windows, which is exactly a centred trapezoidal box of
    width W.
* **Subsampling** keeps every r-th node. The truth is generated r = 4 times finer than the
  PIV grid, with a margin, so that the box never reaches past the data.
* **Noise** is N(0, (σU)²), independent for u, v, every node and every snapshot, and seeded.
* **Gaussian filter** (std s grid points): transfer function exp(−(ksh)²/2).
* **Savitzky–Golay 2-D** (window m × m, order q): at each node, a least-squares fit of a
  polynomial of total degree q, evaluated at the node. It preserves polynomials of degree ≤ q
  exactly, even at the edges, where the window is shifted inwards and the fit is evaluated off
  centre. It therefore removes noise with less peak bias than the Gaussian.

## A11. Error metrics

Pressure is known only up to a constant, so both fields are **demeaned** first (p′ = p − p̄):

    ε_L2  = ‖p′ − p′_ref‖₂ / ‖p′_ref‖₂
    ε_max = max|p′ − p′_ref| / (½ρU²)

Additional metrics:
- **Error map**: |p′ − p′_ref|.
- **Line profile**: the grid row nearest to the vortex core.
- **Core ratio**: (reconstructed minimum within 3h of the reference core) / (reference core
  pressure). A value below 1 means the suction is under-predicted.
- **Convergence slope**: least-squares slope of log ε vs log h.

## A12. Pipeline pseudo-code

```
input : FlowField F (x, y, t, u[nt,ny,nx], v, rho, mu, p_ref?)
        options: method, material_derivative, bc, pin, source, solver, ...

1.  (synthetic PIV only) F <- degrade(truth, box W, subsample r, noise sigma, seed)
2.  (optional)           F <- denoise(F, gaussian | savgol)
3.  n <- chosen snapshot (default: middle)
4.  Du/Dt <- EULERIAN : (u[n+1]-u[n-1])/(2dt) + u*du/dx + v*du/dy    (2nd-order FD)
              LAGRANGIAN: [u(x+u dt, t+dt) - u(x-u dt, t-dt)]/(2dt)  (cubic interp.,
                          Eulerian fall-back outside the window)
5.  G <- -rho*Du/Dt + mu*lap(u)
6.  if method == POISSON:
        f <- div G   (or -rho[u_x^2 + 2u_y v_x + v_y^2])
        b <- f + ghost-point Neumann terms from G on the edges
        R <- trapz(f) - contour(n.G)                      # report
        if Dirichlet edges: replace those rows by p = value
        elif pin == node : replace row k by p_k = 0
        else             : augmented [[A,1],[1^T,0]] (zero mean)
        p <- spsolve(...)  (or CG on -W A)
    elif method == DIRECT: p <- cumulative trapezoid along path 1 / path 2 / average
    elif method == OMNI  : boundary system from N ray directions, then interior ray average
7.  metrics: demean both; eps_L2, eps_max, error map, core ratio, line profile
```

---

# Part B — Corrections to the specification (Rule 4)

| # | Spec says | Problem | What was done |
|---|---|---|---|
| 1 | "Taylor–Green: reconstructed G matches the analytic **−∇p**" | By definition G = −ρDu/Dt + μ∇²u = **+∇p** (Section 3 of the spec itself). | Tests compare G with +∇p. |
| 2 | "Corners take the average of the two edge conditions" | With ghost-point elimination a corner has *two* ghost nodes, and each needs its own condition (∂p/∂x = Gx and ∂p/∂y = Gy). Averaging would throw away information and lose accuracy. | Each ghost node is eliminated with its own edge condition, so both conditions are used exactly. |
| 3 | TG window [0.3, 2.5]² for vortex-core metrics | The window contains **no vortex core**. The cores (pressure minima) are at (0, 0), (π, π), (0, π), (π, 0). The window's centre (π/2, π/2) is a saddle point (pressure maximum). | The spec window is kept for everything else. Core metrics (window and denoising studies) use a second, still non-periodic window [π−1.1, π+1.1]² centred on the core at (π, π) (`tg_core`). |
| 4 | Δt study: "Eulerian degrades faster than Lagrangian" | For the *stationary* TG, ∂u/∂t = −2νu is tiny, so the Eulerian time error is negligible, while the Lagrangian error still grows as Δt². The expectation holds only for *advected* structures, with Lagrangian winning when the vortex velocity w is below the convection speed c (A2.2). | Added a **convected Taylor–Green** case: an exact NS solution by Galilean invariance, amplitude A = 0.4, c = 1. The stationary TG is kept as a control case. |
| 5 | BC study: "Dirichlet edge reduces error" (Bernoulli) | Bernoulli is valid only on steady irrotational edges; TG and Kovasznay have none. Also, even an *exact* Dirichlet value on a single edge can raise the error (Part D). | Added control variants with **exact** p on one edge and on two opposite edges, and Bernoulli on one and on two edges. The results are reported as found. |
| 6 | Box filter of width W via `uniform_filter` | For even widths the moving average is off-centre by half a cell (a phase error). | Even widths average the two half-shifted windows, giving a centred trapezoidal box of exact width W (still `uniform_filter`). |
| 7 | "`RegularGridInterpolator`: cubic if available" | On SciPy ≥ 1.13 the default cubic build uses an iterative solver with atol = 1e-6. The accuracy stalls near 1e-5 and the Lagrangian error then grows under refinement (seen on SciPy 1.15.3). | The direct solver (`solver=spsolve`) is requested when supported. Old SciPy versions fall back automatically. |
| 8 | pin = node "replace one row" | Correct, but the incompatibility R is then dumped at that node as a point source. | Both options are implemented as specified. **pin = mean is the default**, and R is always reported. |
| 9 | CFD export with OpenFOAM `sample` | Writing ~500 full field dumps for post-processing would take gigabytes. | A run-time `sets` function object (the same sampling engine as `sample`) writes only the 81 × 81 window every 0.05. |
| 10 | CFD start-up | A perfectly symmetric mesh and inflow can delay vortex shedding for a long time. | The inlet has v = 0.1 for t < 5, ramped to 0 at t = 6. From t = 6 on it is exactly U = (1, 0), as specified. |
| 11 | Code layout "exactly this" | The common interface, the `reconstruct()` function, shared study code, the summary generator and the case generator need a home. | Added `fields/field.py`, `core/pipeline.py`, `studies/common.py`, `analysis/summary.py`, `cfd/generate_cases.py`, `cfd/run_cases.sh`, `cfd/export_window.py`, `conftest.py`, `tests/test_piv_and_io.py`. Nothing in the listed layout was removed. |
| 12 | Dependencies "NumPy, SciPy, Matplotlib, pytest only" | Running the notebook needs Jupyter; `.h5` JHTDB files need h5py. | Both are optional; the core code imports neither. |
| 13 | "Deliver every file in full in the response" | This project was written straight to disk, where every file exists in full. Pasting ~6 000 lines into chat would add nothing. | All files are complete on disk, with no placeholders or TODOs. |
| 14 | Omni-directional "Liu & Katz 2006" | The original uses a virtual boundary and iteration; this is a rotating-ray variant (A8). | Documented as a variant **[verify]**. |

---

# Part C — File-by-file explanations and viva questions

## `fields/field.py`

**Purpose.** Defines `FlowField`, the common interface that every data source returns and
every algorithm consumes.

| Function / member | What it does |
|---|---|
| `FlowField` (dataclass) | Holds x, y, t, u, v (shape [nt, ny, nx]), rho, mu, p_ref, U_ref, name, meta. |
| `__post_init__` | Converts inputs to float arrays and calls `validate`. |
| `validate` | Checks the array shapes, and that x, y, t are uniform and increasing (the finite differences assume uniform spacing). |
| `nt, ny, nx, dx, dy, dt, steady, nu` | Convenience properties. `dt` is NaN for a single steady snapshot; `nu = mu/rho`. |
| `copy_with` | `dataclasses.replace`: a new object with some fields changed (arrays not copied). |
| `subset_time` | Keeps selected snapshots (e.g. n−k, n, n+k for the Δt study). |
| `grids` | 2-D coordinate arrays via `meshgrid(indexing="xy")`, shape (ny, nx). |

**Viva.**
1. *Why a common interface?* The degradation, integration and metric code need not know
   whether the data are analytic, CFD or JHTDB. Adding a new source means filling one
   dataclass.
2. *Why check uniform spacing?* The stencils (e.g. (f_{i+1} − f_{i−1})/2h) assume a constant h.
   On a non-uniform grid they lose accuracy or are simply wrong.
3. *What does `indexing="xy"` do?* `X[j, i] = x[i]` and `Y[j, i] = y[j]`, which matches the
   [ny, nx] layout.

## `fields/analytic.py`

**Purpose.** The exact solutions used for validation.

| Function | Equation |
|---|---|
| `taylor_green_uvp` | u = c_x + A cos ξ sin η F, v = c_y − A sin ξ cos η F, p = −(ρA²/4)(cos 2ξ + cos 2η)F². Here F = e^{−2νt}, ξ = x − c_x t, η = y − c_y t (Galilean shift c, amplitude A; the spec case is c = 0, A = 1). |
| `taylor_green_grad_p` | ∂p/∂x = (ρA²/2) sin 2ξ F², ∂p/∂y = (ρA²/2) sin 2η F² (used to test G). |
| `taylor_green_on_grid`, `taylor_green` | Sample on given grids / on an n × n window with nt snapshots centred on t_center. |
| `kovasznay_lambda` | λ = Re/2 − √(Re²/4 + 4π²). |
| `kovasznay_uvp` | u = 1 − e^{λx}cos 2πy, v = (λ/2π)e^{λx} sin 2πy, p = (ρ/2)(1 − e^{2λx}). |
| `kovasznay_grad_p` | ∂p/∂x = −ρλe^{2λx}, ∂p/∂y = 0. |
| `kovasznay_on_grid`, `kovasznay` | Single steady snapshot with μ = ρ/Re. |

**Viva.**
1. *Why is Taylor–Green an exact NS solution?* Each velocity component is a Laplacian
   eigenfunction (∇²u = −2u). So ∂u/∂t = ν∇²u when F = e^{−2νt}, and the convective term is a
   gradient, balanced exactly by −∇p/ρ.
2. *Why is the convected version exact?* Navier–Stokes is Galilean invariant. With
   u′(x, t) = c + u(x − ct, t), the terms ∂u′/∂t = u_t − c·∇u and (u′·∇)u′ = ((c + u)·∇)u add
   up to the original material derivative.
3. *Why a non-periodic window?* On a full period the boundary terms cancel and the Neumann
   treatment would never be tested.
4. *What sets λ in Kovasznay flow?* Substituting the ansatz into steady NS gives
   λ² − Re λ − 4π² = 0. The negative root gives a decaying wake.
5. *Why does Kovasznay's pressure depend on x only?* ∇p = −ρ(u·∇)u + μ∇²u collapses to
   −ρλe^{2λx} x̂ for this special solution.

## `fields/cfd_loader.py`

| Function | What it does |
|---|---|
| `DEFAULT_CFD_NPZ` | `results/cfd/cylinder_re100_window.npz`; the environment variable `PIV_CFD_NPZ` overrides it. |
| `cfd_available` | True if the export exists (studies skip CFD otherwise). |
| `save_cfd_npz` | Writes x, y, t, u, v, p, rho, nu, U_inf, D and a JSON metadata string (compressed). |
| `load_cfd` | Reads that file into a `FlowField` with μ = ρν and U_ref = U∞. |

**Viva.** *Why store the kinematic pressure?* OpenFOAM's incompressible solvers solve for
p/ρ, and with ρ = 1 the two are identical. *Why `allow_pickle=False`?* To load only plain
arrays; pickled objects could run arbitrary code.

## `fields/jhtdb_loader.py` (optional)

| Function | What it does |
|---|---|
| `load_jhtdb_npz` | Loads a converted slice (x, y, t, u, v, p [, w]). |
| `load_jhtdb_h5` | Reads plane `z_index` of a cutout HDF5 file (needs h5py; dataset names **[verify]**). |
| `load_jhtdb_slice` | Chooses the reader by file suffix. |
| `out_of_plane_terms` | For planes z−Δz, z, z+Δz: M_x = −ρ w ∂u/∂z + μ ∂²u/∂z² and the same for v, which are the terms a planar reconstruction misses. |

**Viva.**
1. *What is the out-of-plane error?* Planar PIV has no w or ∂/∂z. The 2-D G therefore misses
   −ρw∂u/∂z + μ∂²u/∂z². In 3-D turbulence these terms are as large as the in-plane ones, so the
   error is O(1) and does not go away with refinement.
2. *How would you quantify it?* Take a three-plane cutout, compute M with
   `out_of_plane_terms`, and compare ‖M‖ with ‖G_3D‖, or add M to G and reconstruct again.

## `piv/degrade.py`

| Function | What it does |
|---|---|
| `box_filter` | Centred 2-D moving average over `width` points: `uniform_filter` for odd widths; for even widths the average of the four half-shifted windows, giving a centred trapezoidal box. |
| `subsample` | Every `stride`-th node, then `crop` nodes removed from each side. |
| `add_noise` | u, v += N(0, (σ U_ref)²), seeded with `default_rng(seed)`. |
| `degrade` | Box filter → subsample → noise; p_ref is only subsampled. |
| `degrade_analytic` | Builds the truth r = 4× finer with a margin, degrades it, and crops back to exactly `linspace(x0, x1, n)`. |

**Viva.**
1. *Why generate the truth finer than the PIV grid?* A PIV vector is an average over the
   interrogation window, so the underlying field must be resolved below h.
2. *Why is the reference pressure never filtered?* We want the error with respect to the
   *true* pressure, including the bias that filtering causes.
3. *Why fixed seeds?* Reproducibility. The same seeds are reused across parameter values
   ("common random numbers"), so differences between parameter values are less noisy.
4. *Is fixed σ realistic across window sizes?* Not fully. In real PIV the random error drops
   for larger windows (more particles). This model follows the spec (fixed σ), so W only adds
   bias here.

## `piv/denoise.py`

| Function | What it does |
|---|---|
| `gaussian_denoise` | `gaussian_filter(f, sigma=size, mode="nearest")`. |
| `savgol_exponents` | Monomials xᵃyᵇ with a + b ≤ order. |
| `savgol_coefficients` | Vandermonde matrix V of the window offsets; least-squares operator C = pinv(V). |
| `savitzky_golay_2d` | Polynomial coefficients for every in-array window (`sliding_window_view` @ Cᵀ). Each node uses the nearest window and evaluates the polynomial at its own offset, so edges are handled properly. |
| `denoise_field` | Applies a denoiser to u and v of every snapshot. |

**Viva.**
1. *Why does denoising hurt at σ = 0?* It also smooths the true signal, adding bias and no
   benefit. SUMMARY §7 shows higher error at 0 % for both filters.
2. *Why is Savitzky–Golay better than Gaussian on Kovasznay at low σ?* SG reproduces
   quadratics exactly, so its bias on smooth peaks is smaller. The Gaussian always smooths
   curvature.
3. *What is `pinv(V)`?* The Moore–Penrose pseudo-inverse. It gives the least-squares
   polynomial coefficients from the window values.

## `core/derivatives.py`

| Function | Stencil |
|---|---|
| `_moveaxis_last` | Moves the differentiation axis last so one code path serves both x and y. |
| `first_derivative` | Interior (f_{i+1} − f_{i−1})/2h; ends (∓3f₀ ± 4f₁ ∓ f₂)/2h. |
| `second_derivative` | Interior (f_{i+1} − 2f_i + f_{i−1})/h²; ends (2f₀ − 5f₁ + 4f₂ − f₃)/h². |
| `gradient`, `laplacian`, `divergence` | Combinations of the above. |
| `time_derivative` | Central (fⁿ⁺¹ − fⁿ⁻¹)/2Δt; one-sided second-order at the first and last snapshot. |

**Viva.**
1. *Derive the one-sided second-derivative stencil.* Write f₁, f₂, f₃ as Taylor series about
   x₀ and demand that f₀, f′, f‴ cancel and the f″ coefficient equals 1. The solution is
   (2, −5, 4, −1)/h², with error (11/12)h² f⁗.
2. *Why must boundary stencils be second order?* The Neumann data n·G and the source term on
   the boundary rows come from them. A first-order stencil there would drop the global order
   to 1.
3. *Why not `np.gradient`?* `np.gradient(edge_order=2)` would give the same first derivative,
   but writing it out makes every coefficient visible. The second derivative isn't available
   in NumPy anyway.

## `core/material_derivative.py`

| Function | What it does |
|---|---|
| `PressureGradient` | Container: Gx, Gy, Dux, Duy, time index, info. |
| `_check_time_index` | Defaults to the middle snapshot; 0 for steady data. |
| `convective_acceleration` | (u u_x + v u_y, u v_x + v v_y). |
| `eulerian_material_derivative` | ∂u/∂t (central, or one-sided at the ends) + convective term; ∂u/∂t = 0 if steady. |
| `_make_interpolator` | RGI over (y, x). Tries cubic with `solver=spsolve`, then cubic without it (old SciPy), then linear. |
| `lagrangian_material_derivative` | Pseudo-tracking (A2.2) with Eulerian fall-back where the tracked points leave the window. |
| `pressure_gradient` | G = −ρ Du/Dt + μ∇²u (the viscous term can be switched off for tests). |

**Viva.**
1. *Why is the interpolator built on (y, x) and not (x, y)?* The arrays are [ny, nx], so the
   first axis is y.
2. *What limits the Lagrangian scheme as Δt → 0?* Interpolation error and noise divided by
   2Δt.
3. *Why fall back to Eulerian near the edges?* The tracked points leave the measured region,
   and extrapolation is unreliable.
4. *Which terms carry the most noise at fine grids?* The viscous term, ~μσ/h² (Part D).

## `core/poisson.py`

| Function | What it does |
|---|---|
| `neumann_second_difference_1d` | 1-D matrix [1, −2, 1]/h² with first row [−2, 2] and last row [2, −2] (ghost points eliminated). |
| `trapezoid_weights_1d` | ½, 1, …, 1, ½. |
| `assemble_neumann_laplacian` | A = I_y ⊗ D_xx + D_yy ⊗ I_x (Kronecker products), plus the weights w = w_y ⊗ w_x. |
| `edge_mask` | Boolean mask of an edge's nodes. |
| `source_div_G` | f = ∇·G. |
| `source_velocity` | f = −ρ(u_x² + 2u_y v_x + v_y²). |
| `neumann_rhs` | b = f ± 2G/h on the edges (both at corners). |
| `compatibility_residual` | Trapezoidal ∬f − ∮n·G, and that value divided by (∬\|f\| + ∮\|n·G\|). |
| `bernoulli_pressure` | p∞ + ½ρ(U∞² − u² − v²). |
| `_cg_solve` | Jacobi-preconditioned CG (handles SciPy's `rtol`/`tol` rename) and counts iterations. |
| `solve_poisson` | Dirichlet rows, pin-node row replacement, or the augmented zero-mean system (spsolve); or the symmetric reduced system for CG. |
| `poisson_pressure` | Builds the source term and calls `solve_poisson`. |

**Viva.**
1. *Why is the all-Neumann system singular, and how do you fix it?* A·1 = 0: pin a node, or
   add a zero-mean constraint (A5).
2. *What does the compatibility residual mean physically?* It is the net "source" of the
   Neumann problem that cannot be balanced by boundary flux. It is non-zero when G is not an
   exact gradient.
3. *Why is W·A symmetric?* It is the finite-volume form with half and quarter cells on the
   boundary: every interior flux appears once in each of its two cells, with the same
   coefficient.
4. *Why CG only on −W A?* CG needs a symmetric positive (semi-)definite matrix. The raw
   ghost-point matrix A is not symmetric.
5. *What does λ in the augmented system represent?* A uniform source correction that makes
   the data compatible.

## `core/direct_integration.py`

| Function | What it does |
|---|---|
| `_integrate_from_first_corner` | Paths 1 and 2 from node [0, 0] with `cumulative_trapezoid`. |
| `direct_integration` | Flips the arrays so that the chosen corner becomes [0, 0]. The reversed coordinates make dx negative, so the integration direction stays correct. Returns path 1, path 2 and their average. |

**Viva.** *Why do the paths differ?* curl G ≠ 0 for measured data. *Why does error grow along
the path?* The trapezoidal sum of noisy G is a random walk: error ∝ σ_G h √(number of steps).
*Why can path 1 and path 2 have different mean errors?* See Part D: the error of a path
depends on the flow direction relative to that path. Mirroring the data swaps the two errors
exactly.

## `core/omni_integration.py`

| Function | What it does |
|---|---|
| `bilinear` | Bilinear interpolation on the uniform grid; points outside by round-off are clamped. |
| `boundary_nodes` | Boundary indices counter-clockwise from the SW corner, plus arc length s. |
| `arc_length` | Perimeter coordinate of a boundary point (edge chosen by the nearest side). |
| `_perimeter_weights` | Periodic linear-interpolation weights along the perimeter. |
| `_exit_distance` | Distance along a direction to the rectangle edge. |
| `_ray_integrals` | Trapezoidal ∫G·e ds over rays, with rays sorted by length and processed in memory-bounded chunks. |
| `omni_integration` | Boundary system (augmented, zero mean), then the interior average over N directions. Reports the correction μ, the RMS spread of the single-direction estimates, and the time. |

**Viva.** *Why is omni less noise-sensitive than one path?* Averaging N nearly independent
path errors reduces the variance. *Why a linear system for the boundary?* Each boundary value
depends on other boundary values through the rays. *Cost vs Poisson?* About 30× slower at 64²
(0.4 s vs 15 ms in the noise study), with similar accuracy.

## `core/pipeline.py`

| Function | What it does |
|---|---|
| `Reconstruction` | Result container: p, p_ref, G, time index, method, diagnostics (compatibility, solver, timings), extra (the direct paths). |
| `reconstruct` | Computes G, then Poisson / direct / omni. `bc`: neumann, bernoulli (Dirichlet from Bernoulli on `dirichlet_edges`) or exact (Dirichlet from p_ref, an oracle control). Records wall-clock times. |

**Viva.** *Why an "exact" BC option?* It is a control: it shows what a *correct* Dirichlet
edge (e.g. from pressure taps) would achieve, separately from the Bernoulli modelling error.

## `analysis/metrics.py`

| Function | What it does |
|---|---|
| `demean` | p − mean(p). |
| `eps_l2`, `eps_max`, `error_map` | A11. |
| `locate_minimum` | Minimum inside a region, skipping a 2-node band at the edges. |
| `line_profile` | Row nearest to y = y_line. |
| `core_pressure` | Reference core = minimum of demeaned p_ref; reconstructed = minimum within 3h of it. Returns the ratio and the normalised error. |
| `convergence_slope` | `polyfit(log h, log ε, 1)`. |
| `standard_metrics` | ε_L2 and ε_max in one dict. |

**Viva.** *Why search within 3h for the reconstructed minimum?* A small displacement of the
core should not be counted as an amplitude error. *Can the noisy core ratio exceed 1?* Yes.
Taking the minimum of a noisy field is biased low (a selection bias), which is seen in the
window study (1.02).

## `analysis/plots.py`

| Function | What it does |
|---|---|
| `PALETTE`, `STYLE`, `style_of` | Colour-blind-checked categorical palette in fixed order. Each entity (method, scheme, BC, case) has a fixed colour *and* marker. |
| `CMAP_DIVERGING`, `CMAP_SEQUENTIAL` | Blue–grey–red map for signed pressure; one-hue blue map for error magnitude. |
| `set_style` | Shared rcParams (fonts, thin lines, light grid, 300-dpi export). |
| `save_figure` | PNG (300 dpi) + PDF into `results/figures`; `PIV_RESULTS_DIR` redirects the output. |
| `add_slope_triangle` | Reference line ∝ h^order. |
| `plot_series` | Line + markers ± std band. |
| `field_row` | Several fields with one shared colour scale and one colour bar. |
| `error_panel` | Error map with its own sequential colour bar. |

**Viva.** *Why shared colour limits?* Otherwise two different fields can look identical.
*Why markers as well as colours?* So that identity does not depend on colour vision alone.

## `analysis/summary.py`

Builds `results/SUMMARY.md` from the CSVs. Each `section_*` function turns one study into
tables (mean ± std) plus automatically generated, checkable sentences: fitted slopes, optimum
h, which scheme is lower at the largest Δt, the σ values at which denoising helps, and so on.

**Viva.** *How do you guarantee there are no invented numbers?* Every number in the summary is
computed from the CSVs by this script. Re-running the studies regenerates everything.

## `studies/common.py`

Case factory (`analytic_truth`, `make_case`, `cfd_piv`, `get_cfd`), the evaluation wrapper
(`evaluate`: reconstruct → metrics → one CSV row), statistics (`mean_std`), CSV input/output
(`write_csv`, `read_csv`), and the shared command line (`base_parser`, `seeds_for`).

Baseline: n = 64, W = h, σ = 1 %, Δt = 0.05, refine = 4.

**Viva.** *What are the case keys?* See the module docstring: tg, tg_core, tgc, kov20, kov40,
cfd.

## Study scripts (`studies/run_*.py`)

Each script parses its arguments, loops over cases × parameter × seeds, writes one CSV row per
reconstruction, and saves its figures.

| Script | Sweep | Figures |
|---|---|---|
| `run_convergence.py` | n = 16…256; clean (point samples) vs noisy (W = h, σ = 1 %); Poisson / direct / omni | fig01 (TG), fig01b (all cases), fig02 (Kovasznay fields) |
| `run_noise.py` | σ = 0…10 % | fig03 |
| `run_dt.py` | Δt = 0.01…0.2 (CFD: k × 0.05), Eulerian vs Lagrangian | fig07 |
| `run_bc.py` | neumann (pin mean or node), Bernoulli / exact on 1 or 2 edges | fig08b |
| `run_method.py` | direct paths 1 / 2 / average, Poisson (div G or velocity source), omni | fig08a, fig08 (with BC) |
| `run_window.py` | W = 0, 1, 2, 4 h | fig09a |
| `run_denoise.py` | none / Gaussian / Savitzky–Golay × σ | fig09b, fig09 (with window) |
| `run_cfd_validation.py` | one CFD snapshot, error over one period, GIF | fig05, fig05b, fig06, fig10 |

**Viva for the studies.**
1. *How is the optimum h found?* argmin of the mean noisy ε_L2 over the tested spacings.
2. *Why ≥ 10 seeds?* To estimate the mean and standard deviation of a random error. The
   standard error of the mean falls as 1/√N.
3. *Why vary one parameter at a time?* To attribute each effect cleanly (a local sensitivity
   study around the baseline).

## `run_all.py`

Runs every study in dependency order (BC before method, window before denoise), then the CFD
grid-independence and validation steps (skipped without data), then the summary.
`--quick` uses 3 seeds and coarse sweeps; `--only` selects studies.

## `cfd/generate_cases.py`

| Function | What it does |
|---|---|
| `MeshLevel`, `LEVELS` | Cell counts and gradings per mesh level (coarse / medium / fine = 6 690 / 26 760 / 107 040 cells; ~2× refinement per direction). |
| `blocks` | 12 hex blocks: base quad, cells, simpleGrading. |
| `check_topology` | Asserts counter-clockwise blocks, equal cell counts and grading on shared edges, every outer edge in exactly one patch, and outward patch normals. |
| `first_last_cell` | Geometric grading: first = L(k−1)/(kⁿ−1), with k = ratio^{1/(n−1)}. |
| `foam_header`, `block_mesh_dict`, `sample_window_dict`, `control_dict` | Writers for the OpenFOAM dictionaries. |
| `write_case`, `main` | Write each case directory plus `mesh_info.json`. |

**Viva.**
1. *Why an O-grid?* It gives body-fitted, nearly orthogonal cells, clustered radially at the
   wall, for the boundary layer.
2. *Why one cell in z with `empty` patches?* That is OpenFOAM's way of making a 2-D case.
3. *Why Aref = 0.1?* Aref = D × span, and the span is the cell thickness 0.1.
4. *How do you check the Courant number?* `maxCo 0.5` with `adjustTimeStep yes` lets
   pimpleFoam pick Δt; the logs report Co.

## `cfd/run_cases.sh`, `cfd/cases/*/Allrun`

Generate → blockMesh → checkMesh → pimpleFoam (optionally in parallel with decomposePar and
mpirun) for the three levels, then the grid-independence and export scripts.

## `cfd/export_window.py`

| Function | What it does |
|---|---|
| `window_grid` | The 81 × 81 grid from `mesh_info.json`. |
| `fields_in_name`, `read_set_file` | Handle both the separate (`window_p.xy`, `window_U.xy`) and combined (`window_p_U.xy`) raw layouts. |
| `to_grid_indices` | Map samples to grid indices by their coordinates (so point order doesn't matter). |
| `read_time_dir`, `load_openfoam_samples` | Every sample time → arrays [nt, ny, nx]. |
| `uniform_in_time` | Checks that dt is uniform; otherwise interpolates in time. |
| `load_fluent_exports` | Fluent ASCII fallback, linear `griddata` onto the window. |

## `cfd/grid_independence.py`

| Function | What it does |
|---|---|
| `coefficient_files`, `read_coefficients` | Header-based parsing of the ESI / Foundation force-coefficient files; restarts concatenated. |
| `uniform_resample`, `strouhal` | Uniform resampling, Hann window, zero-padded FFT, parabolic peak refinement → St = fD/U. |
| `analyse` | St, mean Cd, Cl amplitude and RMS, and a periodicity check (first vs second half). |
| `richardson` | p = ln(e₃₂/e₂₁)/ln 2, extrapolated value, GCI = 1.25\|e₂₁/f₁\|/(2ᵖ − 1). |
| `table_markdown`, `plot`, `main` | Table (markdown + figure), Cl(t), spectra, `cfd_summary.json`. |

**Viva.**
1. *Why resample before the FFT?* The adaptive time step makes the samples non-uniform.
2. *Why a Hann window?* To reduce spectral leakage from a non-integer number of periods.
3. *What is the GCI?* A conservative error band on the fine-grid value, based on Richardson
   extrapolation with a safety factor of 1.25.
4. *Why could your St be higher than the unbounded-flow value?* The ±10 D symmetry sides give
   5 % blockage, which accelerates the flow past the cylinder.

## `cfd/FLUENT_SETTINGS.md`

Equivalent Fluent set-up (laminar transient pressure-based solver, ρ = 1, μ = 0.01, the same
boundary conditions, second-order schemes, Co ≤ 0.5, force reports, ASCII export) plus how to
convert the exports.

## Tests (`tests/`)

| File | What it proves |
|---|---|
| `test_derivatives.py` | Slopes 1.85–2.3 in the max norm, *including the boundary nodes*; exact for quadratics and cubics; time derivative exact for quadratics. |
| `test_poisson_exact.py` | p = sin x cos y with exact Neumann data converges at order ≈ 2 for pin mean or node, direct or CG. Also: the null space; wᵀA = 0; W·A symmetric; the trapezoidal residual equals dx·dy·wᵀb; the residual is O(h²); Dirichlet edges; linear p reproduced exactly. |
| `test_taylor_green.py` | G = +∇p to second order (stationary, convected, amplitude 0.4); full reconstruction second order for all three methods; the convected TG satisfies NS; Lagrangian is second order in Δt. |
| `test_kovasznay.py` | G and the reconstruction converge at second order (both source forms); without μ∇²u the error does not converge; steady Lagrangian option works. |
| `test_integration_methods.py` | All methods exact for linear p; exact ∇p integrates correctly; the methods agree on clean data; noisy G shows path dependence, and omni / Poisson beat a single path; Bernoulli edge values applied correctly. |
| `test_piv_and_io.py` | Box filter centred, normalised and preserves linear fields; noise seeding and std; degrade grid, crop and reference; SG preserves quadratics; the denoisers reduce noise; npz round trip; export parser (both layouts, shuffled points); St, Cd and Cl from a synthetic force file; Richardson recovers order 2; mesh topology checks. |

**Viva.** *How do you test second order?* Refine h by 2 and fit the log–log slope. The tests
accept 1.8–2.3 because of pre-asymptotic effects on coarse grids.

## `results.ipynb`

A walkthrough of the validation ladder that only imports package functions. Each step has a
markdown cell before it (what is tested, which equation) and after it (what the result shows).
`QUICK = True` re-runs fast sweeps into `results/notebook_quick/`; `QUICK = False` loads the
full results.

---

# Part D — Non-obvious findings (with explanations)

All numbers are in `results/SUMMARY.md` (seeds 0–9, fixed; the runs are deterministic).

1. **The noisy error has an optimum grid spacing.** For TG and Kovasznay the Poisson noisy
   error is lowest at n = 64 (h ≈ 0.03) and rises again for finer grids. The reason is noise
   in the derivatives:
   - the viscous term contributes ~μ√20 σU/h²;
   - the convective term contributes ~U σU/h;
   - the time derivative contributes σU/(√2Δt), independent of h.

   At h = 0.035 the viscous noise is already comparable to the convective noise, and it grows
   as h⁻².
2. **Direct integration: path 2 is worse than path 1 for Taylor–Green.** This is not a bug:
   mirroring the data (x ↔ y) swaps the two errors exactly. Noise in the convective and
   viscous terms comes from the *same* snapshot, so the two are correlated. Whether they add
   or cancel when integrated along a path depends on the flow direction relative to the path.
   Reversing the velocity flips the sign of the viscous and time-derivative noise but not of
   the convective noise, and swaps which path is better.
3. **The velocity-source Poisson form is much worse on noisy data** (TG ε ≈ 0.42 vs 0.08 for
   the div-G form). Its interior source has no time-derivative noise, but its Neumann data
   (from G) do, so source and boundary data no longer come from one vector field. The
   least-squares structure of A3 is lost. With noise-free boundary G the velocity form was
   actually *better* (checked during development). Its weakness is the mismatch, not the
   identity.
4. **One exact Dirichlet edge can hurt, two opposite edges help.** With a single Dirichlet
   edge the Neumann-data error is anchored on one side and accumulates across the domain. Two
   opposite edges bracket it (TG noisy: 1 edge 0.093, 2 edges 0.044, all-Neumann 0.082).
   Bernoulli edges are much worse on both analytic cases, as expected for rotational flow.
5. **Eulerian vs Lagrangian depends on w/c (A2.2).** Stationary TG: Eulerian is flat and
   Lagrangian grows. Convected TG with w/c = 0.4: at Δt = 0.2 Lagrangian is ~3× lower. With
   noise, both errors *fall* with Δt over this range, because the ~1/Δt noise amplification
   dominates the truncation error.
6. **Window averaging has a small effect on the Taylor–Green core** (core ratio 0.997 at W = 4h),
   because the vortex (size ~π) is much larger than W = 0.14. The effect grows with W/(core
   size), so it is expected to be much stronger for the cylinder-wake vortices (core radius
   ~0.5 D, W up to 0.4 D).

---

# Part E — Convergence and error analysis (Lecture 11 add-on)

All defaults are unchanged. `spsolve` stays the default Poisson solver and `order=2` the
default scheme. The existing results were re-run and their existing CSV columns compared
value by value with the previous run; they are identical. Everything below is opt-in through
`reconstruct(..., solver=, order=, tol=, max_iter=, omega=, residual_scale=, L_ref=)` or runs
in the new study scripts.

## E0. Corrections to the add-on specification

| # | Spec says | Problem | What was done |
|---|---|---|---|
| 1 | "Jacobi" | On the all-Neumann Poisson problem the checkerboard vector (−1)^(i+j) is an exact eigenvector of D⁻¹A with eigenvalue 2, so classic Jacobi's iteration matrix has an eigenvalue exactly −1. That mode never decays and the residual **stalls** (e.g. at 5.0e-7 on 32², above tol = 1e-8). In pseudo-time terms, Δτ = h²/4 is exactly the explicit stability limit (E3). | Jacobi takes a relaxation ω_J. The default **ω_J = 0.95** (Δτ just inside the limit) damps the mode at a 5 % cost. Classic ω_J = 1 is still available, and Figure 11 shows its stall. |
| 2 | ω_opt ≈ 2/(1 + sin(πh/L)) | Young's formula for the **Dirichlet** model problem. With all-Neumann walls the slowest non-constant mode is cos(πx/L), with half the eigenvalue, so ρ_J is closer to 1 and the optimum is larger (e.g. 1.932 instead of 1.905 at 64²). | `omega_opt(..., neumann=True)` computes Young's ω from the correct ρ_J; the SOR default uses the Neumann estimate. Figure 15 marks both formulas and the measured best ω. |
| 3 | "Same matrix A for all solvers" | CG needs a symmetric positive (semi-)definite matrix, and the ghost-point A is not symmetric. | CG runs on S = −WA: the same equations, each row multiplied by its trapezoidal weight, so the same solution. Its recorded residual is that of the original A p = b. |
| 4 | "Same pinning" | The zero-mean augmented system [[A, 1], [1ᵀ, 0]] has a zero on the diagonal, so it cannot be iterated with Jacobi, Gauss–Seidel or SOR. | The iterative solvers solve the equivalent problem: b is made compatible with the same uniform correction λ = wᵀb/wᵀ1 that the augmented system applies, and the mean is removed after every iteration. The solution is identical to the direct one (tests: agreement < 1e-6). |
| 5 | "Gauss–Seidel" | A lexicographic sweep needs a Python loop over nodes, which is far too slow. | **Red-black** ordering: all red nodes (i + j even), then all black ones, each half-sweep vectorised. For the 5-point stencil red nodes couple only to black nodes, and the ordering is consistently ordered, so ρ_GS = ρ_J² and Young's SOR theory hold. |
| 6 | "`order=1`: first-order Neumann BCs" | The corner "average of the two edge conditions" is again unnecessary. | The ghost node is eliminated with the one-sided condition p₋₁ = p₀ − h g. The matrix stays symmetric (CG still works), its left null vector becomes **1**, and corners use both conditions. Interior derivatives are forward differences. The time difference stays central, because the study varies Δx at fixed Δt. |
| 7 | "Unnormalised vs normalised L2 across grids" | The ghost-point boundary rows of b contain flux terms 2G/h that grow as 1/h, so even L2ₙ of the full b grows (about N^0.2). | Figure 12b uses the interior rows, the rows of the 5-point Laplacian, where L1ₙ and L2ₙ are flat. The full-b L2ₙ is plotted too, dotted and labelled. |
| 8 | "Iterations vs N, 32²–256², each solver" | Jacobi and Gauss–Seidel at 256² need ~2–4×10⁵ iterations (Gauss–Seidel took 7.5 min on the M2; Jacobi about twice as long). | Default `--slow-max-n 128` for Jacobi and Gauss–Seidel; `--slow-max-n 256` gives the full range. SOR and CG always run up to 256². |
| 9 | "Stop the steady run when ‖uⁿ⁺¹ − uⁿ‖ < ε" | OpenFOAM cannot stop on a field-change norm without custom code. | The optional Re = 20 case stops through `runTimeControl` on the initial residuals **[verify syntax]**. The lecture's change-norm criterion is evaluated in Python on the sampled window, and Figure 17(f) shows both. |
| 10 | Periodic wake and ‖uⁿ⁺¹ − uⁿ‖ | The Re = 100 wake is periodic, so the change never goes to zero. | A cycle-to-cycle periodicity criterion (E7). |
| 11 | Richardson with three meshes | The refinement ratios come from cell counts and need not be exactly equal. | The Celik et al. (2008) fixed-point procedure for p, which reduces to ln(e₃₂/e₂₁)/ln r for equal ratios **[verify reference]**. |
| 12 | "Existing results must not change" | New metrics add CSV columns. | Existing columns are value-identical. New columns: err_l1n, err_l2n, err_linf, err_pct. The Re = 100 OpenFOAM dictionaries regenerate byte-identically. `mesh_info.json` gains two descriptive keys (`kick`, `residual_stop`), and its `sample_start` for the coarse and medium meshes changes from null to 175 (unused, since those meshes do not sample). |

## E1. Norms — `analysis/norms.py`

| Function | Equation | Meaning |
|---|---|---|
| `l1(r)` | Σ\|rᵢ\| | total absolute imbalance; all cells weigh the same |
| `l2(r)` | √(Σ rᵢ²) | Euclidean length; large components dominate |
| `linf(r)` | max\|rᵢ\| | worst cell only |
| `l1_normalised(r)` | (1/n) Σ\|rᵢ\| | mean absolute imbalance per cell |
| `l2_normalised(r)` | √((1/n) Σ rᵢ²) | RMS imbalance per cell |
| `scaled(r, S)` | r/S | relative ("percentage") imbalance |
| `residual_scale(b, kind, ρ, U, L)` | `"rhs"`: S = L2ₙ(b); `"physical"`: S = ρU²/L² | the characteristic scale |
| `norm_summary(r)` | all five | convenience |

**Why normalise?** For a residual of constant size c on n cells, L1 = nc and L2 = √n·c, while
L1ₙ = L2ₙ = c. The same physical imbalance therefore gives larger L1 and L2 on finer grids, so
a fixed tolerance is harder to meet the finer the grid.

**Ordering.** L1ₙ ≤ L2ₙ ≤ L∞ for any vector (power-mean inequality).

**Lecture example** (`tests/test_norms.py`). For r = (0.0001, 0.0004, 0.0006, 0.002):
- L1 = 0.0031;
- L2 = √4.53e-6 = 0.00213;
- L2ₙ = √(4.53e-6/4) = 0.00106, the "representative" ≈ 0.0011 on the slide.

**Scales.**
- `"rhs"`: with p⁰ = 0 the first residual is b, so the scaled residual starts at exactly 1 and
  "tol = 1e-8" means "reduce the imbalance by 10⁸". This is the Poisson analogue of scaling
  finite-volume residuals by the face flux.
- `"physical"`: ∇²p ~ ρU²/L² (dynamic pressure over length squared), independent of the
  numerics.

**Metrics refactor.** `eps_l2 = L2(e)/L2(p′_ref)` and `eps_max = L∞(e)/(½ρU²)` now call these
functions (bit-identical results). `error_norms` adds L1ₙ, L2ₙ and L∞ of the error, plus
err % = 100·eps_L2.

**Viva.**
1. *Why does L2 favour large components?* Squaring gives large entries more weight, so one bad
   cell dominates. L1 weights all cells equally.
2. *Which norm would you use for a convergence check?* A normalised one (L2ₙ), scaled by a
   physical or source scale. Add L∞ if you care about the worst cell.
3. *What does a scaled residual of 1e-3 mean?* The average imbalance per cell is 0.1 % of the
   characteristic source magnitude.
4. *Is L∞ normalised?* It needs no normalisation: the maximum does not grow with n.

## E2. Scaled residuals and "levelling off"

The lecture's point: the exact value of the scaled residual matters less than the shape of
its history. It should fall steadily and **level off at a low value**. A curve that levels off
high means a problem, such as an incompatible or inconsistent system or an undamped mode.
Classic Jacobi on the all-Neumann problem (E0 #1) is a real example, where the residual
levels off at an iteration-independent floor.

## E3. Iterative solvers — `core/iterative.py`

**The system.** `LinearSystem` holds:
- A, with identity rows for prescribed unknowns;
- b, made compatible for the zero-mean problem;
- the fixed-node mask and the initial guess;
- `zero_mean`: whether to remove the mean every iteration;
- the symmetrising weights w;
- the grid shape.

`core/poisson.py:poisson_system` and `linear_system` build it from exactly the data the direct
solver uses.

**Jacobi** (`jacobi`): p^{k+1} = p^k + ω D⁻¹(b − A p^k), with D = diag(A).

*Pseudo-time view.* For dx = dy = h every row of A has diagonal −4/h², so D⁻¹ = −h²/4 and

    p^{k+1} = p^k + ω (h²/4) (A p^k − b) = p^k + Δτ (∇²_h p^k − f),   Δτ = ω h²/4

This is exactly explicit (forward-Euler) time marching of the pseudo-transient equation
∂p/∂τ = ∇²p − f toward its steady state ∇²p = f. The steady-state criterion ‖p^{k+1} − p^k‖ < ε
is the lecture's rule for breaking the (pseudo-)time loop.

*Stability.* Explicit diffusion amplifies each eigenmode by 1 + Δτλ, where λ ∈ [−8/h², 0] are
the eigenvalues of the 5-point Laplacian. Stability requires |1 − 8Δτ/h²| ≤ 1, i.e.
Δτ ≤ h²/4.
- At Δτ = h²/4 (classic Jacobi, ω = 1) the mode with λ = −8/h², the checkerboard (−1)^(i+j),
  has amplification −1: it flips sign every step and never decays.
- With Dirichlet walls the most negative eigenvalue is slightly above −8/h², so Jacobi still
  converges, very slowly for that mode.
- With all-Neumann walls the checkerboard is an exact eigenvector (each ghost-point boundary
  row also gives −8/h²), so Jacobi stalls at a floor set by the checkerboard content of the
  solution.
- ω = 0.95 (Δτ = 0.95h²/4) gives amplification −0.9.

*Convergence rate.* The slowest (smoothest) mode has amplification 1 − ωμ_min, with
μ_min ≈ (πh/L)²/4 for the Neumann problem. The iteration count to reduce the residual by a
factor F is about ln F/μ_min ∝ (L/h)² ∝ N. That is why Jacobi's iteration count grows like
N¹ (Figure 14).

**Gauss–Seidel and SOR** (`gauss_seidel`, `sor`): p_i ← p_i + ω(bᵢ − (Ap)ᵢ)/A_ii, using the newest
values, in red-black order. For consistently ordered matrices (Young):
- ρ_GS = ρ_J², so Gauss–Seidel needs half the Jacobi iterations (measured: 12,328 vs 25,366 at
  64²);
- SOR's spectral radius is minimised at

      ω_opt = 2 / (1 + √(1 − ρ_J²)),   giving ρ_SOR = ω_opt − 1,

  which makes the iteration count scale like √N instead of N.

Assumptions: 5-point Laplacian, consistent ordering, real Jacobi eigenvalues with
|μ| ≤ ρ_J < 1 (excluding the null-space pair).
- Dirichlet model problem on a square: ρ_J = cos(πh/L), giving ω_opt = 2/(1 + sin(πh/L)), the
  spec formula (`omega_opt_dirichlet_formula`).
- All-Neumann walls (dx = dy): the lowest non-constant mode varies in one direction only, so
  μ_min = sin²(πh/2L) and ρ_J = cos²(πh/2L) (`jacobi_spectral_radius(neumann=True)`). Then
  √(1 − ρ_J²) ≈ πh/(√2 L) and ω_opt ≈ 2/(1 + πh/(√2 L)), which is larger than the Dirichlet
  value.
- **Kahan's theorem**: ρ_SOR ≥ |ω − 1|, so SOR cannot converge for ω ≤ 0 or ω ≥ 2. Figure 15b
  shows ω = 2 (no decay) and ω = 2.1 (divergence).

**Conjugate gradient** (`conjugate_gradient`). Jacobi-preconditioned CG on S = −WA (symmetric
positive semi-definite); prescribed unknowns are eliminated. Per iteration:

    α = (r·z)/(d·Sd);  x += α d;  r −= α Sd;  z = M⁻¹r;  β = (r·z)_new/(r·z)_old;  d = z + β d

- CG minimises the error in the S-norm over a growing Krylov space. The iteration count scales
  like √κ(S) ~ L/h ∝ √N.
- For the singular zero-mean problem, CG from p⁰ = 0 converges because the right-hand side is
  compatible (orthogonal to the null space).

**History and stopping** (`IterationHistory`, `_Recorder`). Every iteration records L1ₙ, L2ₙ,
L∞ and the scaled L2ₙ of r = b − Ap (free rows), plus the change L2ₙ(p^{k+1} − p^k). The loop
stops when:
- the scaled L2ₙ < tol (**converged**);
- k = max_iter;
- the residual grows by 10¹⁰ or becomes non-finite (**diverged**).

`to_rows` writes the history to CSV. An optional `callback(p, k)` lets a study compute extra
quantities, such as the true error in Figure 13.

**Measured at 64²** (TG, tol 1e-8):

| Solver | Iterations | Time |
|---|---:|---:|
| Jacobi (ω = 0.95) | 25,366 | 1.3 s |
| Gauss–Seidel | 12,328 | 0.9 s |
| SOR (ω = 1.932) | 314 | 0.02 s |
| CG | 180 | 0.01 s |
| direct `spsolve` | – | ~0.015 s |

The spec formula ω = 1.905 needs 628 SOR iterations.

**Viva.**
1. *Show that Jacobi is explicit pseudo-time marching.* See the derivation above: Δτ = ωh²/4.
2. *Why does classic Jacobi stall here but not on a Dirichlet problem?* The Neumann
   checkerboard mode has amplification exactly −1 (at the stability limit). Dirichlet walls
   shift the extreme eigenvalue so that |amplification| < 1.
3. *Why is SOR so much faster?* Over-relaxation pushes the spectral radius from 1 − O(h²) to
   1 − O(h). Iterations scale like √N instead of N.
4. *Why does SOR diverge for ω ≥ 2?* Kahan: ρ_SOR ≥ |ω − 1| ≥ 1.
5. *Why not always use CG?* Here you could (it is the fastest iterative method). It requires a
   symmetric positive (semi-)definite matrix, which is why it runs on −WA and not on A. The
   direct `spsolve` is still fastest at these sizes and stays the default.

## E4. Residual vs change norm; the steady-state criterion

For Jacobi, p^{k+1} − p^k = ωD⁻¹r^k, so

    ‖p^{k+1} − p^k‖ = Δτ ‖r^k‖,   with Δτ = ω h²/4 ≈ 3·10⁻⁴ at 64²

The change is always about 10³–10⁴ times smaller than the residual. Worse, the error satisfies
‖p^k − p*‖ ≈ ‖p^{k+1} − p^k‖/(1 − ρ), and with ρ = 1 − O(h²) the factor 1/(1 − ρ) is huge. A
small change per iteration therefore does **not** mean a converged solution when convergence
is slow. Figure 13 marks the iteration where the change first drops below ε = 1e-6: the scaled
residual and the true error are still orders of magnitude larger there. That is why the solvers
stop on the scaled residual, not on the change.

For physical time-stepping (CFD) the same criterion ‖uⁿ⁺¹ − uⁿ‖ < ε ends a *steady* run. A
periodic flow never satisfies it (E7).

**Viva.** *Is ‖uⁿ⁺¹ − uⁿ‖ < ε enough to declare steady state?* Only together with a residual
check. When the rate of change is small because the time step (or pseudo-time step) is small,
the change can be tiny while the solution is still far from steady.

## E5. First-order scheme — `order=1`

**Derivatives** (`core/derivatives.py`):
- forward first difference (f_{i+1} − fᵢ)/h = f′ + (h/2)f″ + …, so O(h);
- forward second difference (f_{i+2} − 2f_{i+1} + fᵢ)/h² = f″ + h f‴ + …, so O(h).

`gradient`, `laplacian` and `divergence` take `order`. `core/material_derivative.py` passes it
to the convective and viscous terms; the time difference stays central.

**Neumann BC** (`core/poisson.py`). With the one-sided condition (p₀ − p₋₁)/h = g, the boundary
row becomes (p₁ − p₀)/h² − g/h. A Taylor expansion shows this row approximates f with an
**O(1)** local error (it equals p″/2 + …). Because that error sits on a strip of width O(h), the
global error is **O(h)**. The matrix rows [−1, 1]/h² make A symmetric with left null vector 1,
so the compatibility condition becomes Σb = 0, applied through the same augmented or
uniform-correction machinery.

**Why the global slope is 1.** The largest truncation error is O(h) in both G and the boundary
treatment, and the Poisson solve preserves the order of its data. The test
`tests/test_scheme_order.py` measures slope ≈ 1 for `order=1` and ≈ 2 for `order=2` on clean TG.

**Viva.**
1. *Why keep the time difference central in `order=1`?* The comparison is about Δx at fixed Δt.
   A first-order time error would add an h-independent floor and flatten the curve.
2. *Why does an O(1) boundary truncation error still give O(h) convergence?* It acts only on a
   layer of O(h) width. The same argument makes the second-order scheme's O(h) boundary rows
   globally O(h²).
3. *Why are fitted slopes over all grids lower than the asymptotic order?* The coarsest grids
   are pre-asymptotic (higher-order terms still matter). Fit over the finest grids.

## E6. Grid-convergence plot (Figure 16) — `studies/run_order.py`

- Error (%) = 100·ε_L2, plotted against Δx on log–log axes with Δx **reversed** (finer to the
  right).
- Straight lines of slope p mean error ≈ CΔxᵖ.
- Slope triangles: legs Δx → 2Δx and error → 2ᵖ·error.
- Slopes are fitted over the finest five grids (`slope_fine`); `slope_all` is also reported.
- **1 % target:** `crossing()` returns the first tested Δx with error < 1 % and the log–log
  interpolated crossing between the bracketing grids:

      Δx* = exp(ln Δx_{k−1} + t (ln Δx_k − ln Δx_{k−1})),   t = (ln 1 − ln e_{k−1}) / (ln e_k − ln e_{k−1})

  If a scheme never reaches 1 %, the crossing is extrapolated with its asymptotic slope,
  Δx* = Δx_min (1/e_min)^{1/p}, and marked "(extrapolated)".
- **Noisy overlay** (σ = 1 %, W = h, 10 seeds): the error stops falling and rises again on fine
  grids, from noise amplified as 1/Δx (convective term) and 1/Δx² (viscous term).

**Viva.**
1. *What resolution do you need for 1 % error?* Read the guide lines; the order-2 scheme
   reaches it on a much coarser grid.
2. *How would you halve the error with each scheme?* Order 1: refine Δx by 2. Order 2: by √2.
3. *Why not always refine?* With measurement noise the error rises again (see the dashed
   lines). The noisy optimum is reported in SUMMARY §10.

## E7. CFD solver residuals, change norm and periodicity — `cfd/parse_log.py`, `cfd/change_norm.py`

**`parse_log(path)`** reads `Time = …` lines (ESI or Foundation form, `…s`) and, for each time
step, the first `Solving for Ux/Uy/p, Initial residual = …`. For p, the first of the several
pressure solves per step is kept (foamLog's `p_0`).
- The Courant number is printed *before* its `Time =` line, so it is held until that line.
- `read_foamlog_dir` reads the `foamLog` output instead.
- The **initial residual** is the normalised imbalance of the equation at the start of the time
  step, before the linear solver acts. In a steady flow it falls toward zero; in the periodic
  wake it levels off.

**`change_norms(u, v)`**: L2ₙ of (uⁿ⁺¹ − uⁿ, vⁿ⁺¹ − vⁿ)/U∞ between saved snapshots.
**`first_below`**: the first time the change drops below ε.

**`cycle_statistics` / `periodicity_check`.**
- Cycles are cut at the upward zero crossings of Cl − mean, with crossing times refined by
  linear interpolation.
- Each cycle gives a period T_k and an amplitude A_k = (max − min)/2.
- The flow counts as periodic from the first cycle after which |A_k − A_{k−1}|/A_k and
  |T_k − T_{k−1}|/T_k both stay below `tol` (default 1e-3) for `n_consecutive` (3) cycles.
- This replaces the steady-state criterion, which a periodic flow can never meet.

**Optional steady case `steady_re20`** (`cfd/generate_cases.py`, `OPTIONAL_LEVELS`):
- medium mesh, ν = 0.05 (Re = 20, below the shedding onset near Re ≈ 47 **[verify]**), no kick;
- window sampled every 0.5 D/U from t = 0;
- a `runTimeControl` function object ends the run when every initial residual is below 1e-5
  **[verify syntax]**;
- run it with `sh cfd/run_cases.sh N --steady`.

`change_norm.py` then reports when ‖uⁿ⁺¹ − uⁿ‖ first drops below ε and when OpenFOAM stopped.

**Viva.**
1. *Why do the Re = 100 residuals not fall to machine zero?* The solution changes every time
   step (periodic shedding), so each new step starts with a finite imbalance.
2. *How do you know the wake is fully developed?* The cycle-to-cycle changes of Cl amplitude and
   period stay below tolerance for several cycles (Figure 17e).
3. *Initial vs final residual?* Initial = before the linear solve (measures how much the
   solution must change). Final = after the solve (measures the linear-solver tolerance).

## E8. Richardson extrapolation and GCI — `cfd/grid_independence.py`

- `representative_h(cells)` = √(fluid area/cells).
- The refinement ratios are r₂₁ = h₂/h₁ and r₃₂ = h₃/h₂ (1 = fine).

`richardson(f1, f2, f3, r21, r32)` (Celik et al. 2008 **[verify]**):

    e21 = f2 − f1, e32 = f3 − f2, s = sign(e32/e21)
    p = | ln|e32/e21| + ln((r21^p − s)/(r32^p − s)) | / ln r21        (fixed-point iteration)
    f_ext = (r21^p f1 − f2)/(r21^p − 1)
    GCI_fine = 1.25 |e21/f1| / (r21^p − 1),    GCI_medium = 1.25 |e32/f2| / (r32^p − 1)
    asymptotic ratio = GCI_medium / (r21^p GCI_fine)   (≈ 1 in the asymptotic range)

It also reports the % change coarse → medium and medium → fine. Oscillatory convergence
(s < 0) is flagged, and no order is reported. The test builds f(h) = f₀ + Chᵖ with p = 2, 1 and
1.5 (unequal ratios) and recovers p and f₀ to 1e-8.

Figure 18: St and mean C_D vs cells, with the extrapolated value dashed and ±GCI·f error bars
on the medium and fine meshes. Literature values stay marked [verify].

**Viva.**
1. *What does GCI = 0.5 % on the fine mesh mean?* The fine-mesh value is estimated to lie within
   0.5 % of the grid-converged value, with safety factor 1.25.
2. *Why a safety factor?* The observed order comes from only three meshes and may not be
   asymptotic. 1.25 is Roache's value for three-grid studies.
3. *What if p comes out very different from 2?* The meshes are probably not in the asymptotic
   range, or the quantity has other error sources (e.g. time step, FFT resolution for St).
   Refine or check.

## E9. New study scripts, figures and tests

| File | What it does |
|---|---|
| `studies/run_residuals.py` | Figures 11–15. `baseline_system` builds the TG system and its direct solution. Parts: `part_history` (all solvers + classic Jacobi), `part_norms` (SOR norms; interior vs full norms across grids), `part_change` (Jacobi residual/change/true error through a callback), `part_sor` (ω scan, divergence runs), `part_grid` (iterations vs N with fitted exponents). |
| `studies/run_order.py` | Figure 16 (TG, Kovasznay): `run`, `crossing`, `summarise`, `plot_case`. |
| `studies/run_cfd_convergence.py` | Figures 17–18 from `parse_log`, `change_norm` and `grid_independence` outputs; every missing part is skipped. |
| `analysis/plots.py` | New entity styles (solvers, norms, orders, residual fields, meshes) and helpers: `draw_slope_triangle`, `threshold_guides`, `plot_history`. |
| `analysis/summary.py` | New sections 9 (solvers, norms, SOR) and 10 (order), plus the CFD time-loop criteria and the Richardson table in §8. |
| Tests | `test_norms.py`, `test_iterative.py`, `test_scheme_order.py`, `test_cfd_convergence.py` (37 new tests; 87 in total). |
