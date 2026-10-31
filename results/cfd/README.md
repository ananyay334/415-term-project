# CFD results (cylinder wake, Re = 100)

From the OpenFOAM v2606 runs (MinGW build, native Windows, 4 cores) on the coarse, medium and fine
meshes, post-processed with `cfd/grid_independence.py`, `cfd/export_window.py`,
`cfd/parse_log.py`, `cfd/change_norm.py` and `run_all.py`.

| File | Content |
|---|---|
| `grid_independence.md` | St, mean Cd, Cl amplitude per mesh; Richardson extrapolation and GCI |
| `figures/fig04_cfd_grid_independence.png` | grid table, Cl(t) and its spectrum |
| `figures/fig05_cfd_wake.png` | CFD pressure, reconstruction and error map (t = 187.5) |
| `figures/fig05b_cfd_error_history.png` | reconstruction error over one shedding period |
| `figures/fig06_vortex_core_line.png` | pressure along a line through a vortex core |
| `figures/fig07_dt_with_cfd.png` | error vs time step, including the CFD wake |
| `figures/fig17_cfd_residuals.png` | solver residuals and periodicity check |
| `figures/fig18_cfd_grid_convergence.png` | St and Cd vs mesh with GCI |
| `figures/paraview_fine_t200.png` | ParaView 6.2 view of the fine-mesh pressure at t = 200 |

Key numbers: fine-mesh St = 0.1695 (extrapolated 0.1697, observed order 2.48, GCI 0.13 %),
mean Cd = 1.388, Cl amplitude = 0.345; periodic from t = 76.3.

The raw `cylinder_re100_window.npz` export and the case time directories are not committed
(regenerate with `sh cfd/cases/<level>/Allrun N` and `python cfd/export_window.py`).
