# Results summary (generated automatically)

Generated 2026-09-30 08:41 by `analysis/summary.py` from the CSV files in `results/`. All numbers below come from those files; re-running the studies regenerates them.

## 1. Grid spacing h (run_convergence.py)

Clean = exact point samples; noisy = W = 1 h, sigma = 1 % (mean ± std over seeds).

### Taylor-Green

| h | n | Direct (avg.) clean | Omni clean | Poisson clean | Direct (avg.) noisy | Omni noisy | Poisson noisy |
|---|---|---|---|---|---|---|---|
| 0.1467 | 16 | 0.00702 | 0.0147 | 0.00941 | 0.255 ± 0.043 | 0.0895 ± 0.012 | 0.135 ± 0.029 |
| 0.07097 | 32 | 0.00157 | 0.0036 | 0.0027 | 0.294 ± 0.071 | 0.0599 ± 0.0089 | 0.101 ± 0.036 |
| 0.03492 | 64 | 0.00036 | 0.000931 | 0.000716 | 0.424 ± 0.11 | 0.0627 ± 0.0066 | 0.0823 ± 0.019 |
| 0.01732 | 128 | 8.54e-05 | 0.000216 | 0.000184 | 0.901 ± 0.13 | - | 0.162 ± 0.061 |
| 0.008627 | 256 | 2.07e-05 | 5.62e-05 | 4.66e-05 | 2.22 ± 0.6 | - | 0.202 ± 0.046 |

- Direct (avg.): fitted clean slope = **2.06**; noisy optimum h = **0.1467** (eps_L2 = 0.255); the noisy error rises again for smaller h.
- Omni: fitted clean slope = **1.97**; noisy optimum h = **0.07097** (eps_L2 = 0.0599); the noisy error rises again for smaller h.
- Poisson: fitted clean slope = **1.88**; noisy optimum h = **0.03492** (eps_L2 = 0.0823); the noisy error rises again for smaller h.

### Kovasznay Re=20

| h | n | Direct (avg.) clean | Omni clean | Poisson clean | Direct (avg.) noisy | Omni noisy | Poisson noisy |
|---|---|---|---|---|---|---|---|
| 0.1333 | 16 | 0.224 | 0.0671 | 0.164 | 0.263 ± 0.033 | 0.122 ± 0.0038 | 0.208 ± 0.011 |
| 0.06452 | 32 | 0.0809 | 0.0154 | 0.0535 | 0.197 ± 0.021 | 0.0406 ± 0.0053 | 0.0787 ± 0.015 |
| 0.03175 | 64 | 0.0205 | 0.00324 | 0.0148 | 0.404 ± 0.098 | 0.0666 ± 0.0071 | 0.059 ± 0.017 |
| 0.01575 | 128 | 0.00499 | 0.000708 | 0.00387 | 1.32 ± 0.44 | - | 0.133 ± 0.042 |
| 0.007843 | 256 | 0.00122 | 0.000166 | 0.00099 | 3.48 ± 0.68 | - | 0.216 ± 0.046 |

- Direct (avg.): fitted clean slope = **1.87**; noisy optimum h = **0.06452** (eps_L2 = 0.197); the noisy error rises again for smaller h.
- Omni: fitted clean slope = **2.13**; noisy optimum h = **0.06452** (eps_L2 = 0.0406); the noisy error rises again for smaller h.
- Poisson: fitted clean slope = **1.82**; noisy optimum h = **0.03175** (eps_L2 = 0.059); the noisy error rises again for smaller h.

### Kovasznay Re=40

| h | n | Direct (avg.) clean | Omni clean | Poisson clean | Direct (avg.) noisy | Omni noisy | Poisson noisy |
|---|---|---|---|---|---|---|---|
| 0.1333 | 16 | 0.288 | 0.096 | 0.098 | 0.336 ± 0.043 | 0.153 ± 0.0055 | 0.15 ± 0.02 |
| 0.06452 | 32 | 0.0979 | 0.023 | 0.0235 | 0.252 ± 0.034 | 0.0514 ± 0.0073 | 0.0721 ± 0.021 |
| 0.03175 | 64 | 0.0245 | 0.00512 | 0.00536 | 0.455 ± 0.11 | 0.0715 ± 0.0065 | 0.0671 ± 0.021 |
| 0.01575 | 128 | 0.00595 | 0.00118 | 0.00127 | 1.42 ± 0.48 | - | 0.15 ± 0.051 |
| 0.007843 | 256 | 0.00146 | 0.000282 | 0.000308 | 3.63 ± 0.7 | - | 0.236 ± 0.054 |

- Direct (avg.): fitted clean slope = **1.89**; noisy optimum h = **0.06452** (eps_L2 = 0.252); the noisy error rises again for smaller h.
- Omni: fitted clean slope = **2.07**; noisy optimum h = **0.06452** (eps_L2 = 0.0514); the noisy error rises again for smaller h.
- Poisson: fitted clean slope = **2.04**; noisy optimum h = **0.03175** (eps_L2 = 0.0671); the noisy error rises again for smaller h.

## 2. Noise level sigma (run_noise.py)

n = 64, W = 1 h, dt = 0.05; eps_L2 mean ± std over seeds.

### Taylor-Green

| sigma (%) | Direct (avg.) | Omni | Poisson |
|---|---|---|---|
| 0 | 0.000584 | 0.00116 | 0.000841 |
| 0.5 | 0.212 ± 0.056 | 0.0314 ± 0.0033 | 0.0411 ± 0.0097 |
| 1 | 0.424 ± 0.11 | 0.0627 ± 0.0066 | 0.0823 ± 0.019 |
| 2 | 0.849 ± 0.23 | 0.125 ± 0.013 | 0.165 ± 0.038 |
| 5 | 2.13 ± 0.57 | 0.313 ± 0.033 | 0.412 ± 0.095 |
| 10 | 4.29 ± 1.2 | 0.628 ± 0.068 | 0.841 ± 0.19 |

- Direct (avg.): log-log slope of eps_L2 vs sigma (sigma > 0) = **1.00** (1 = linear growth); mean wall-clock time 0.3 ms per reconstruction.
- Omni: log-log slope of eps_L2 vs sigma (sigma > 0) = **1.00** (1 = linear growth); mean wall-clock time 410.8 ms per reconstruction.
- Poisson: log-log slope of eps_L2 vs sigma (sigma > 0) = **1.01** (1 = linear growth); mean wall-clock time 15.0 ms per reconstruction.
- Lowest error at sigma = 10 %: **Omni**.

### Kovasznay Re=20

| sigma (%) | Direct (avg.) | Omni | Poisson |
|---|---|---|---|
| 0 | 0.0234 | 0.00653 | 0.0179 |
| 0.5 | 0.202 ± 0.048 | 0.0336 ± 0.0034 | 0.0323 ± 0.011 |
| 1 | 0.404 ± 0.098 | 0.0666 ± 0.0071 | 0.059 ± 0.017 |
| 2 | 0.807 ± 0.2 | 0.133 ± 0.015 | 0.116 ± 0.031 |
| 5 | 2.02 ± 0.5 | 0.333 ± 0.037 | 0.289 ± 0.073 |
| 10 | 4.03 ± 1 | 0.666 ± 0.074 | 0.581 ± 0.15 |

- Direct (avg.): log-log slope of eps_L2 vs sigma (sigma > 0) = **1.00** (1 = linear growth); mean wall-clock time 0.3 ms per reconstruction.
- Omni: log-log slope of eps_L2 vs sigma (sigma > 0) = **1.00** (1 = linear growth); mean wall-clock time 415.5 ms per reconstruction.
- Poisson: log-log slope of eps_L2 vs sigma (sigma > 0) = **0.97** (1 = linear growth); mean wall-clock time 15.1 ms per reconstruction.
- Lowest error at sigma = 10 %: **Poisson**.

## 3. Time step dt (run_dt.py)

### Taylor-Green

| dt | Eulerian clean | Lagrangian clean | Eulerian noisy | Lagrangian noisy |
|---|---|---|---|---|
| 0.01 | 0.000716 | 0.000684 | 0.267 ± 0.071 | 0.258 ± 0.056 |
| 0.02 | 0.000716 | 0.000686 | 0.142 ± 0.039 | 0.142 ± 0.043 |
| 0.05 | 0.000716 | 0.000738 | 0.0823 ± 0.019 | 0.0961 ± 0.027 |
| 0.1 | 0.000716 | 0.00147 | 0.0719 ± 0.014 | 0.0724 ± 0.018 |
| 0.2 | 0.000716 | 0.00422 | 0.0698 ± 0.014 | 0.0695 ± 0.014 |

- Clean data, largest dt = 0.2: Eulerian 0.000716, Lagrangian 0.00422 -> **Eulerian** is lower; growth from smallest to largest dt: Eulerian x1, Lagrangian x6.17.
- Noisy data, Eulerian: lowest mean error at dt = 0.2 (0.0698).
- Noisy data, Lagrangian: lowest mean error at dt = 0.2 (0.0695).

### Taylor-Green, convected (A=0.4)

| dt | Eulerian clean | Lagrangian clean | Eulerian noisy | Lagrangian noisy |
|---|---|---|---|---|
| 0.01 | 0.00176 | 0.00106 | 1.62 ± 0.57 | 1.7 ± 0.43 |
| 0.02 | 0.00163 | 0.00109 | 0.928 ± 0.39 | 1.05 ± 0.36 |
| 0.05 | 0.000951 | 0.00113 | 0.593 ± 0.27 | 0.596 ± 0.28 |
| 0.1 | 0.00345 | 0.00156 | 0.522 ± 0.23 | 0.506 ± 0.23 |
| 0.2 | 0.018 | 0.00577 | 0.497 ± 0.22 | 0.476 ± 0.22 |

- Clean data, largest dt = 0.2: Eulerian 0.018, Lagrangian 0.00577 -> **Lagrangian** is lower; growth from smallest to largest dt: Eulerian x10.2, Lagrangian x5.44.
- Noisy data, Eulerian: lowest mean error at dt = 0.2 (0.497).
- Noisy data, Lagrangian: lowest mean error at dt = 0.2 (0.476).

## 4. Boundary conditions (run_bc.py)

### Taylor-Green

| variant | eps_L2 clean | eps_L2 noisy | eps_max noisy | time (ms) |
|---|---|---|---|---|
| All-Neumann, pin mean | 0.000716 | 0.0823 ± 0.019 | 0.178 | 15.1 |
| All-Neumann, pin node | 0.00582 | 0.112 ± 0.053 | 0.375 | 8.0 |
| Bernoulli, 1 edge | 0.291 | 0.304 ± 0.013 | 0.719 | 7.8 |
| Exact p, 1 edge | 0.00258 | 0.093 ± 0.033 | 0.146 | 7.9 |
| Bernoulli, 2 edges | 0.537 | 0.539 ± 0.0083 | 1.1 | 8.2 |
| Exact p, 2 edges | 0.000795 | 0.0444 ± 0.0083 | 0.113 | 8.1 |

- Lowest noisy error: **Exact p, 2 edges** (0.0444).

### Kovasznay Re=20

| variant | eps_L2 clean | eps_L2 noisy | eps_max noisy | time (ms) |
|---|---|---|---|---|
| All-Neumann, pin mean | 0.0148 | 0.059 ± 0.017 | 0.742 | 15.5 |
| All-Neumann, pin node | 0.0171 | 0.0661 ± 0.027 | 0.765 | 8.2 |
| Bernoulli, 1 edge | 0.0297 | 0.0613 ± 0.019 | 0.699 | 8.4 |
| Exact p, 1 edge | 0.0258 | 0.059 ± 0.02 | 0.696 | 9.4 |
| Bernoulli, 2 edges | 0.814 | 0.813 ± 0.0041 | 6.34 | 13.4 |
| Exact p, 2 edges | 0.00201 | 0.0325 ± 0.002 | 0.378 | 7.8 |

- Lowest noisy error: **Exact p, 2 edges** (0.0325).

## 5. Integration methods (run_method.py)

### Taylor-Green

| variant | eps_L2 clean | eps_L2 noisy | eps_max noisy | time (ms) |
|---|---|---|---|---|
| Direct path 1 | 0.000449 | 0.402 ± 0.1 | 0.867 | 0.3 |
| Direct path 2 | 0.000443 | 0.729 ± 0.22 | 0.701 | 0.3 |
| Direct (avg.) | 0.00036 | 0.424 ± 0.11 | 0.653 | 0.3 |
| Poisson | 0.000716 | 0.0823 ± 0.019 | 0.178 | 15.0 |
| Poisson (velocity source) | 0.000708 | 0.418 ± 0.14 | 0.669 | 15.2 |
| Omni | 0.000931 | 0.0627 ± 0.0066 | 0.156 | 423.3 |

- Lowest noisy error: **Omni** (0.0627).

### Kovasznay Re=20

| variant | eps_L2 clean | eps_L2 noisy | eps_max noisy | time (ms) |
|---|---|---|---|---|
| Direct path 1 | 0.0334 | 0.545 ± 0.18 | 2.62 | 0.3 |
| Direct path 2 | 0.0192 | 0.57 ± 0.14 | 3.26 | 0.3 |
| Direct (avg.) | 0.0205 | 0.404 ± 0.098 | 2.35 | 0.3 |
| Poisson | 0.0148 | 0.059 ± 0.017 | 0.742 | 14.6 |
| Poisson (velocity source) | 0.00574 | 0.485 ± 0.15 | 2.98 | 15.1 |
| Omni | 0.00324 | 0.0666 ± 0.0071 | 0.672 | 422.3 |

- Lowest noisy error: **Poisson** (0.059).

- Taylor-Green (noisy): RMS difference between direct paths 1 and 2 = 0.156; RMS spread of single-direction omni estimates = 0.0463 (pressure units).
- Kovasznay Re=20 (noisy): RMS difference between direct paths 1 and 2 = 0.566; RMS spread of single-direction omni estimates = 0.199 (pressure units).

## 6. Interrogation window W (run_window.py)

Core ratio = reconstructed / reference vortex-core pressure (after mean removal); < 1 means the suction peak is under-predicted.

### Taylor-Green (core window)

| W / h | eps_L2 clean | core ratio clean | eps_L2 noisy | core ratio noisy |
|---|---|---|---|---|
| 0 | 8.63e-05 | 1 | 0.0812 ± 0.029 | 1.02 ± 0.01 |
| 1 | 0.000219 | 1 | 0.0812 ± 0.029 | 1.02 ± 0.01 |
| 2 | 0.000817 | 0.999 | 0.0812 ± 0.029 | 1.02 ± 0.01 |
| 4 | 0.00325 | 0.997 | 0.0812 ± 0.029 | 1.02 ± 0.01 |

- Clean core ratio goes from 1 (W = 0 h) to 0.9967 (W = 4 h); decreases monotonically with W.

### Kovasznay Re=20

| W / h | eps_L2 clean | core ratio clean | eps_L2 noisy | core ratio noisy |
|---|---|---|---|---|
| 0 | 0.0148 | - | 0.0586 ± 0.016 | - |
| 1 | 0.0179 | - | 0.059 ± 0.017 | - |
| 2 | 0.0265 | - | 0.0611 ± 0.019 | - |
| 4 | 0.0613 | - | 0.0793 ± 0.024 | - |

## 7. Denoising (run_denoise.py)

### Taylor-Green (core window)

| sigma (%) | none | Gaussian | Savitzky-Golay |
|---|---|---|---|
| 0 | 0.000219 | 0.00677 | 0.00033 |
| 0.5 | 0.0407 ± 0.014 | 0.0149 ± 0.0022 | 0.0154 ± 0.003 |
| 1 | 0.0812 ± 0.029 | 0.0274 ± 0.0045 | 0.0308 ± 0.006 |
| 2 | 0.162 ± 0.057 | 0.0534 ± 0.009 | 0.0616 ± 0.012 |
| 5 | 0.406 ± 0.14 | 0.132 ± 0.023 | 0.155 ± 0.03 |
| 10 | 0.828 ± 0.29 | 0.264 ± 0.047 | 0.312 ± 0.062 |

- Gaussian: lower error than no denoising at sigma = 0.5 %, 1 %, 2 %, 5 %, 10 %; higher at 0 %.
- Savitzky-Golay: lower error than no denoising at sigma = 0.5 %, 1 %, 2 %, 5 %, 10 %; higher at 0 %.

### Kovasznay Re=20

| sigma (%) | none | Gaussian | Savitzky-Golay |
|---|---|---|---|
| 0 | 0.0179 | 0.0576 | 0.0212 |
| 0.5 | 0.0323 ± 0.011 | 0.0578 ± 0.00084 | 0.0218 ± 0.0014 |
| 1 | 0.059 ± 0.017 | 0.0581 ± 0.0017 | 0.0234 ± 0.0026 |
| 2 | 0.116 ± 0.031 | 0.0596 ± 0.0033 | 0.0292 ± 0.0042 |
| 5 | 0.289 ± 0.073 | 0.069 ± 0.0076 | 0.055 ± 0.0065 |
| 10 | 0.581 ± 0.15 | 0.0956 ± 0.012 | 0.104 ± 0.0096 |

- Gaussian: lower error than no denoising at sigma = 1 %, 2 %, 5 %, 10 %; higher at 0 %, 0.5 %.
- Savitzky-Golay: lower error than no denoising at sigma = 0.5 %, 1 %, 2 %, 5 %, 10 %; higher at 0 %.

## 8. CFD cylinder wake

CFD grid-independence results not found (run cfd/run_cases.sh).

## 9. Iterative Poisson solvers and residual norms (run_residuals.py)

Taylor-Green 64x64, all-Neumann (zero mean), tol on the scaled L2_n residual.

| solver | omega | iterations | wall-clock (s) | converged | final scaled residual |
|---|---|---:|---:|---|---:|
| Jacobi (omega 0.95) | 0.95 | 25366 | 1.25 | True | 9.99e-09 |
| Gauss-Seidel | 1 | 12328 | 0.928 | True | 9.99e-09 |
| SOR | 1.932 | 314 | 0.0238 | True | 9.78e-09 |
| CG | nan | 180 | 0.0135 | True | 9.67e-09 |
| Jacobi, classic (omega 1) | 1 | 50732 | 2.6 | False | 4.15e-08 |
- Gauss-Seidel needed 2.06x fewer iterations than Jacobi (omega = 0.95).
- SOR needed 80.8x fewer iterations than Jacobi (omega = 0.95).
- CG needed 141x fewer iterations than Jacobi (omega = 0.95).
- Classic Jacobi (omega = 1) did not reach the tolerance: its residual levelled off at 4.15e-08.

| n | N | Jacobi (omega 0.95) it. (s) | Gauss-Seidel it. (s) | SOR it. (s) | CG it. (s) |
|---|---|---|---|---|---|
| 32 | 1024 | 6275 (0.16) | 3049 (0.11) | 152 (0.0057) | 93 (0.0035) |
| 64 | 4096 | 25366 (1.3) | 12328 (0.94) | 314 (0.025) | 180 (0.013) |
| 128 | 16384 | 100753 (47) | 48991 (76) | 644 (0.95) | 333 (1.3) |
| 256 | 65536 | - | - | 1314 (9.5) | 546 (2.5) |

- Jacobi (omega 0.95): iterations grow like N^1.00 (fitted over n = 32, 64, 128).
- Gauss-Seidel: iterations grow like N^1.00 (fitted over n = 32, 64, 128).
- SOR: iterations grow like N^0.52 (fitted over n = 32, 64, 128, 256).
- CG: iterations grow like N^0.43 (fitted over n = 32, 64, 128, 256).

- SOR: fewest iterations (314) at omega = 1.932 among the tested values.
- SOR with omega = 2: scaled residual went from 1 to 6.94 after 300 iterations (no convergence).
- SOR with omega = 2.1: scaled residual went from 1 to 1.09e+10 after 145 iterations (diverged).

- Norms of the same initial residual (interior rows) on grids 32, 64, 128, 256: growth exponents with N - L1 0.99, L2 0.50, L1N -0.01, L2N -0.00.
- Including the ghost-point boundary rows (flux terms 2G/h) the L2_n exponent is 0.24.

## 10. First- vs second-order scheme (run_order.py)

Slopes of error (%) vs dx; 'fine' = fit over the finest grids (asymptotic range).

| case | order | slope (fine) | slope (all grids) | dx where error first < 1 % | 1 % crossing (interpolated) | noisy: min error % at dx |
|---|---|---:|---:|---:|---:|---|
| Taylor-Green | 1 | 0.98 | 0.88 | 0.01222 | 0.01229 | 12.4 % at 0.05 |
| Taylor-Green | 2 | 1.96 | 1.82 | 0.1467 | 0.1526 | 8.23 % at 0.0349 |
| Kovasznay Re=20 | 1 | 1.00 | 0.84 | nan | 0.003338 (extrapolated) | 17.8 % at 0.0455 |
| Kovasznay Re=20 | 2 | 1.96 | 1.70 | 0.02222 | 0.02584 | 5.62 % at 0.0455 |
