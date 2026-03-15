```@meta
CurrentModule = SeaGap
```

# Software Architecture

This page provides an analysis of the SeaGap software architecture, covering module organization, data flow, key algorithms, and design principles.

## Overview

SeaGap is implemented as a Julia package (`module SeaGap`) consisting of more than 60 source files, each responsible for a specific aspect of GNSS-Acoustic positioning or post-processing.
The package follows a text-based, pipeline-oriented design: every major function reads its inputs from text files on disk and writes its results to text files, so individual processing steps can be run independently and inspected at each stage.

## Module Structure

The source files are organized into six logical groups, as reflected in `SeaGap.jl`:

### 1. Core Utilities

Foundational helpers used throughout the package.

| File | Responsibility |
|------|----------------|
| `dateprocessing.jl` | Convert between date strings, cumulative seconds, year fractions, and day-of-year |
| `ll2xy.jl` | Geographic ↔ projected coordinate conversion (Transverse Mercator via GMT) |
| `anttena2tr.jl` | Transform GNSS antenna position to transducer position using lever-arm offset and vessel attitude |
| `interpolate_gps.jl` | Interpolate GPS position/attitude time-series to arbitrary times |
| `unixsort.jl` | Unix-style multi-column sort for Julia arrays and matrices |
| `running_filter.jl` | Running-median and running-average filters |
| `LineFitting.jl` | Weighted least-squares line fitting with confidence intervals |
| `perturbation.jl` | Uniform random perturbation for MCMC step proposals |
| `select_distribution.jl` | Estimate the mode of a distribution for parameter reporting |

### 2. Input / Output

Functions that read raw instrument data and format it into the unified observation structure.

| File | Responsibility |
|------|----------------|
| `read_gnssa.jl` | Parse all input file types: site info, sound-speed profile, GNSS trajectory, travel-time records, transponder positions |
| `obsdata_format.jl` | Combine GNSS and travel-time data into a single `obsdata.inp` observation file |
| `make_initial_grad.jl` | Generate `initial.inp` with MCMC step widths for the sound-speed gradient inversion |
| `make_initial_gradv.jl` | Variant of the above for the full 3-D gradient variant |

### 3. Mathematical Core

The building blocks shared by all positioning methods.

| File | Responsibility |
|------|----------------|
| `traveltime.jl` | Exact ray-traced travel-time (`xyz2tt`), rapid Taylor-series approximation (`xyz2tt_rapid`), version with horizontal gradients (`xyz2ttg_rapid`), and correction-parameter computation (`ttcorrection`) |
| `ntdbasis.jl` | Construction and evaluation of temporally-uniform cubic B-spline bases for Nadir Total Delay (NTD) modelling |
| `simple_inversion.jl` | Unweighted least-squares solver: ``\hat{\bf a}=({\bf H}^T{\bf H})^{-1}{\bf H}^T{\bf d}`` |
| `inv_func.jl` | 1-D Laplacian smoothing matrix used as a regularization operator |
| `perturbation.jl` | Step-proposal functions for MCMC samplers |

### 4. Positioning Methods

The core scientific algorithms; all are iterative non-linear inverse problems solved by Gauss-Newton or MCMC.

| File | Method | Key outputs |
|------|--------|-------------|
| `denoise.jl` | Pre-processing outlier removal (static context) | Cleaned `obsdata.inp` |
| `denoise_kinematic.jl` | Pre-processing outlier removal (kinematic context) | Cleaned `obsdata.inp` |
| `ttres.jl` | Forward-only travel-time residual check | `residual.out` |
| `kinematic_array.jl` | Kinematic 2-D array positioning (shot-group by shot-group) | `kinematic.out` |
| `kinematic_array_3d.jl` | Kinematic 3-D array positioning | `kinematic3d.out` |
| `static_array.jl` | Static 3-D array positioning with B-spline NTD | `solve.out`, `position.out`, `residual.out`, `bspline.out` |
| `static_array_s.jl` | Static positioning with separate B-spline NTD for each transponder | same set of outputs |
| `static_array_AICBIC.jl` | Static positioning swept over a range of B-spline knot counts with AIC/BIC output | `AICBIC.out` |
| `static_array_s_ABIC.jl` | ABIC-regularised static positioning | `ABIC.out` |
| `static_array_grad.jl` | Static positioning with a horizontal sound-speed gradient (linear gradient estimation) | `solve.out` with gradient parameters |
| `static_array_TR.jl` | Static positioning including estimation of transducer–antenna offset | expanded `solve.out` |
| `static_individual.jl` | Individual transponder positioning from a single campaign | `position_individual.out` |
| `static_array_mcmcgrad.jl` | MCMC static positioning with horizontal sound-speed gradient | `mcmc.out`, `mcmcparam.out` |
| `static_array_mcmcgradc.jl` | MCMC variant with constrained gradient | same outputs |
| `static_array_mcmcgradv.jl` | MCMC variant with 3-D velocity gradient | same outputs |

### 5. Post-Processing and Visualisation

More than 20 plot functions generate figures from the output files produced by the positioning methods.
All plotting uses [Plots.jl](https://github.com/JuliaPlots/Plots.jl) and [GMT.jl](https://github.com/GenericMappingTools/GMT.jl).

| Group | Files |
|-------|-------|
| Residuals | `plot_ttres.jl` |
| NTD time-series | `plot_ntd.jl`, `plot_ntd_s.jl`, `plot_ntd_grad.jl`, `plot_ntd_gradv.jl`, `plot_ntdgrad.jl` |
| Positioning track | `plot_prof.jl`, `plot_track.jl`, `plot_multi-track.jl`, `plot_kinematic_array.jl`, `plot_array_each.jl` |
| MCMC diagnostics | `plot_mcmcres.jl`, `plot_mcmcres_grad.jl`, `plot_mcmcres_gradc.jl`, `plot_mcmcres_gradv.jl` |
| MCMC parameters | `plot_mcmcparam.jl`, `plot_mcmcparam_grad.jl`, `plot_mcmcparam_gradv.jl`, `plot_mcmcparam_each.jl` |
| Parameter maps | `plot_gradmap.jl`, `plot_gradmap_grad.jl`, `plot_gradmap_gradv.jl`, `plot_bspline_gradv.jl` |
| Histograms | `plot_histogram.jl`, `plot_histogram_grad.jl`, `plot_histogram_gradv.jl`, `plot_histogram2d.jl`, `plot_histogram2d_grad.jl`, `plot_histogram2d_gradv.jl`, `plot_histogram2d_each.jl` |
| Correlation maps | `plot_cormap.jl`, `plot_cormap_grad.jl`, `plot_cormap_gradv.jl` |
| Model selection | `plot_AICBIC.jl`, `plot_ABIC.jl` |
| Displacement | `convert_displacement.jl`, `plot_displacement.jl` |
| Kinematic | `position_kinematic.jl` |

### 6. Testing / Development Utilities

| File | Purpose |
|------|---------|
| `forward_test.jl` | Compute synthetic travel-times for a user-supplied set of sea-surface positions |
| `forward_gradv.jl` | Variant of forward calculation including horizontal sound-speed gradients |

---

## Data Flow

```
Raw Instrument Data
│
│  site_info.txt          – geographic position of the GNSS-A site
│  tr-ant.inp             – transducer–antenna lever-arm offsets
│  pxp-ini.inp            – initial transponder positions (lon, lat, depth)
│  ss_prof.inp            – sound-speed profile (depth vs. speed)
│  gps.jxyhhpr            – GNSS platform trajectory + vessel attitude
│  pxp-*.jttq             – two-way travel-time observations per transponder
│
▼  [obsdata_format]
│
Unified Observation File (obsdata.inp)
│  – shot epoch, transponder ID, observed two-way travel-time
│  – interpolated sea-surface position and attitude at transmit and receive times
│  – quality flags
│
▼  [denoise]  (optional outlier removal)
│
Cleaned Observations
│
├─▶ [ttres]               – travel-time residuals at fixed transponder positions
│       └─ residual.out
│
├─▶ [kinematic_array]     – shot-group-by-shot-group 2-D positioning
│       └─ kinematic.out
│
└─▶ [static_array]        – whole-campaign 3-D positioning
        ├─ solve.out          (estimated displacement + formal errors)
        ├─ position.out       (array position time-series)
        ├─ residual.out       (normalised travel-time residuals)
        ├─ bspline.out        (B-spline NTD coefficient time-series)
        └─ AICBIC.out         (model-selection statistics)

Post-Processing
│
├─ plot_* functions       – figures from the above output files
└─ convert_displacement   – convert position increments to geodetic displacements
```

---

## Key Algorithms

### Travel-Time Calculation (`traveltime.jl`)

Two complementary approaches are provided:

1. **Exact ray-tracing (`xyz2tt`)** – a shooting method that integrates the ray through a layered sound-speed medium defined by the profile in `ss_prof.inp`.
   It uses the local Earth radius computed from the WGS84 ellipsoid via `localradius(lat)`.

2. **Rapid Taylor approximation (`xyz2tt_rapid`)** – a first-order Taylor expansion about a pre-computed reference ray.
   The correction parameters (computed once by `ttcorrection`) encapsulate the vertical travel-time and the sensitivity to depth and horizontal offset.
   This approximation is used inside iterative solvers where many evaluations per iteration are needed.

3. **Gradient-aware rapid approximation (`xyz2ttg_rapid`)** – extends the rapid approximation to also return horizontal partial derivatives with respect to the array position, enabling gradient-based positioning.

### B-Spline NTD Modelling (`ntdbasis.jl`)

The Nadir Total Delay (NTD) represents the integrated effect of temporal sound-speed fluctuations.
It is modelled as a time-series of cubic B-splines with a uniform knot spacing ``\Delta s``:

```math
\text{NTD}(t) = \sum_{j=1}^{N_{\rm PB}} c_j \, \Phi_j(t)
```

where ``\Phi_j`` are the B-spline basis functions evaluated by `tbspline3` and ``c_j`` are the unknown coefficients estimated in every positioning run.
The number of basis functions ``N_{\rm PB}`` (knots) is a free parameter that controls temporal resolution; `static_array_AICBIC` sweeps over a range of ``N_{\rm PB}`` values and selects the optimal one using AIC and BIC.

### Static Array Positioning (`static_array.jl`)

The observation equation for a single travel-time observation is:

```math
\frac{T_n^{\rm obs}}{M(\xi_n)} = \frac{T_n^{\rm cal}(\Delta\mathbf{p})}{M(\xi_n)} + \sum_{j=1}^{N_{\rm PB}} c_j \Phi_j(t_n)
```

where:
- ``T_n^{\rm obs}`` is the observed two-way travel-time,
- ``M(\xi_n) = 1/\cos(\xi_n)`` is the mapping function (``\xi_n`` = nadir angle),
- ``\Delta\mathbf{p}`` is the unknown 3-D array displacement,
- ``t_n`` is the shot epoch.

The unknowns ``(\Delta\mathbf{p},\, c_1,\ldots,c_{N_{\rm PB}})`` are solved jointly by **Gauss-Newton iteration** with convergence threshold ``\varepsilon = 10^{-4}`` and a maximum of 50 iterations.
The Jacobian with respect to ``\Delta\mathbf{p}`` is computed by finite differences (perturbation ``10^{-4}`` m).

### Kinematic Array Positioning (`kinematic_array.jl`)

Each short time window (shot group) is processed independently.
For group ``k``, the unknowns are the horizontal array displacement ``(\Delta x_k, \Delta y_k)`` and one NTD value ``C_k``:

```math
T_{kn}^{\rm obs} = T_{kn}^{\rm cal}(\Delta\mathbf{p}_k,\, C_k) + \epsilon_{kn}
```

A Gauss-Newton solver is applied to each group separately, with a minimum of `NR = 3` observations required per group.

### MCMC Positioning (`static_array_mcmcgrad*.jl`)

When the horizontal sound-speed gradient is unknown, a Markov Chain Monte Carlo sampler is used.
At each MCMC step, all unknown parameters (array displacement, NTD coefficients, gradient components) are perturbed by a uniform random offset whose width is derived from the posterior errors of the gradient-free solution.
The posterior distribution is sampled over a user-specified number of steps, and the results are written to `mcmc.out` for post-processing.

---

## Design Principles

1. **Text-based pipeline** – every function reads from and writes to plain-text files.
   This makes it straightforward to inspect intermediate results, restart a failed run from any intermediate stage, and integrate SeaGap with shell scripts or other tools.

2. **Single-responsibility modules** – each `.jl` source file defines one primary exported function with the same name as the file (e.g. `static_array.jl` defines `static_array()`).
   Shared helpers are factored into dedicated utility files.

3. **Separation of concerns** – positioning (inverse problem), travel-time physics, NTD basis functions, visualisation, and I/O are implemented in separate files with minimal coupling.

4. **Progressive complexity** – simpler methods (kinematic, basic static) are implemented first; more complex variants (with gradients, MCMC, individual transponders, ABIC regularisation) are layered on top by re-using the same mathematical primitives.

5. **Model-selection support** – rather than prescribing the number of B-spline knots, `static_array_AICBIC` sweeps the parameter space and returns AIC/BIC statistics so the user can make an informed choice.

---

## External Dependencies

| Package | Role |
|---------|------|
| `LinearAlgebra` (stdlib) | Matrix operations, least-squares, LU decomposition |
| `Statistics` (stdlib) | Mean, variance |
| `Random` (stdlib) | Random number generation for MCMC |
| `Dates` (stdlib) | Date/time arithmetic |
| `Printf` (stdlib) | Formatted text output |
| `Distributed` (stdlib) | Parallel computation across workers |
| `Base.Threads` (stdlib) | Shared-memory multithreading |
| `Dierckx` | Spline fitting and interpolation |
| `Optim` | General-purpose numerical optimisation |
| `Plots` | 2-D visualisation |
| `GMT` | Map projections and map-based plots |
| `Distributions` | Statistical distributions (MCMC, mode estimation) |
| `DelimitedFiles` | Reading and writing delimiter-separated files |
| `PDFmerger` | Merging multi-page PDF output |
