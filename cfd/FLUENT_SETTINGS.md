# ANSYS Fluent settings (fallback for the OpenFOAM cylinder case)

These settings reproduce the OpenFOAM set-up in `cfd/generate_cases.py` as closely as
Fluent allows. Menu names follow recent Fluent versions (2022 R1 and later); check them
against your version, marked **[verify]** where they are version-specific.

## 1. Geometry and mesh

| Item | Value |
|---|---|
| Dimension | 2-D planar |
| Domain | x ∈ [−10, 20] D, y ∈ [−10, 10] D; cylinder D = 1 m centred at (0, 0) |
| Mesh | same blocking as the OpenFOAM mesh: an O-grid from the cylinder to the square [−2, 2]², then rectangular blocks. Use the cell counts and gradings in `cfd/cases/<level>/mesh_info.json` (coarse 6 690, medium 26 760, fine 107 040 cells). |
| How | ANSYS Meshing (edge sizing + bias factor = the expansion ratios) or ICEM-CFD blocking. Alternatively convert the OpenFOAM mesh with `foamMeshToFluent`. That utility writes a 3-D mesh, so for a 2-D Fluent run rebuild the mesh in 2-D **[verify]**. |

## 2. General and models

* Solver: **Pressure-based**, **Transient**, 2D Space **Planar**, absolute velocity formulation.
* Models → Viscous: **Laminar**.
* Materials → create a fluid with **density = 1 kg/m³** and **viscosity = 0.01 kg/(m·s)**.
  With D = 1 m and U = 1 m/s this gives Re = ρUD/μ = 100, the same as OpenFOAM's ν = 0.01.

## 3. Boundary conditions

| Zone | Fluent type | Setting | OpenFOAM equivalent |
|---|---|---|---|
| inlet (x = −10) | velocity-inlet | 1 m/s normal to boundary (see §6 for the start-up kick) | fixedValue (1 0 0), p zeroGradient |
| outlet (x = 20) | pressure-outlet | gauge pressure 0 Pa | zeroGradient U, p = 0 |
| top / bottom (y = ±10) | symmetry | – | symmetryPlane |
| cylinder | wall | stationary, no slip | noSlip |

Operating pressure: leave at the default. The reconstruction compares pressures only up to a constant.

## 4. Solution methods and controls

| Item | Setting |
|---|---|
| Pressure–velocity coupling | **PISO** (skewness and neighbour correction 1) or SIMPLEC |
| Gradient | Least-squares cell based |
| Pressure | Second order |
| Momentum | Second-order upwind (like `linearUpwind`) or bounded central differencing |
| Transient formulation | (Bounded) second-order implicit (like OpenFOAM `backward`) |
| Iterations per time step | 20–30, with scaled residuals ≤ 1e-5 for continuity and momentum |
| Time step | fixed, chosen so that the Courant number is ≤ 0.5 everywhere (check the Courant-number contour). An exact divisor of 0.05 s makes exporting easy, e.g. 0.0025 s (fine), 0.005 s (medium), 0.01 s (coarse) **[verify Co]** |
| End time | 200 s (non-dimensional t = 200, about 33 shedding periods) |

## 5. Force monitors (Cd, Cl)

* Reference Values: area = 1 m² (D × 1 m depth for 2-D), length = 1 m, velocity = 1 m/s,
  density = 1 kg/m³.
* Report Definitions → Force → **Drag coefficient** on `cylinder`, force vector (1, 0).
  Add a **Lift coefficient** with force vector (0, 1).
* Report File: write both every time step together with the flow time.
* To reuse `cfd/grid_independence.py`, convert the report file into a whitespace table with
  the header line `# Time Cd Cl`. Save it as
  `cfd/cases/<level>/postProcessing/forceCoeffs1/0/forceCoeffs.dat`; the script reads the
  columns by name.

## 6. Starting the shedding

The OpenFOAM case adds a small transverse inlet velocity v = 0.1 for t < 5, ramped to 0 by
t = 6. In Fluent, use an expression for the inlet y-velocity, e.g.
`if(t < 5 [s], 0.1 [m/s], if(t < 6 [s], 0.1 [m/s] * (6 [s] - t) / 1 [s], 0 [m/s]))`
**[verify expression syntax]**. Alternatively patch v = 0.1 m/s into a small region behind the
cylinder at t = 0.

## 7. Export of the PIV window

* Calculation Activities → Automatic Export → **ASCII**, every *k* time steps so that the
  interval is 0.05 s, starting at t = 175 s and running to 200 s.
* Quantities: x-velocity, y-velocity, pressure. Location: node values, space or comma
  delimited, one file per export.
* Convert to the project format (linear interpolation onto the 81 × 81 window
  x ∈ [1, 5], y ∈ [−2, 2]):

```bash
python cfd/export_window.py --fluent-dir path/to/exports --fluent-dt 0.05
```

The parser looks for the header names `x-coordinate`, `y-coordinate`, `x-velocity`,
`y-velocity` and `pressure` **[verify against your export file]**. Files are ordered by
name, so use zero-padded time-step numbers.
