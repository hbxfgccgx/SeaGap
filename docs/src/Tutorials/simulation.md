```@meta
CurrentModule = SeaGap
```

# Simulating GNSS-Acoustic Data

This tutorial describes how to generate a synthetic GNSS-Acoustic dataset
using `simulate_gnssa()` and then perform static array positioning on the
simulated observations.

## Scenario Overview

The simulated scenario represents a GNSS-Acoustic site in the **Western
Pacific Ocean** at **5000 m water depth**:

| Parameter | Value |
|-----------|-------|
| Site latitude | 20 °N |
| Water depth | 5000 m |
| Number of transponders | 4 |
| Transponder array | Square, half-width 700 m |
| Transducer depth | 2 m below sea surface |
| GNSS antenna height | 25 m above sea surface |

### Vessel Track

The sea-surface vessel track is composed of two parts:

1. **Circle** – radius equal to the water depth (5000 m = 1 × dep), centred
   directly over the seafloor transponder array.  The ship completes one
   full revolution at 3 m s⁻¹ in approximately 2.9 hours.

2. **Cross** (过顶十字) – an east–west leg and a north–south leg, each
   running from one side of the circle to the other (−5000 m → +5000 m)
   and passing **directly over the centre of the array**.  Each leg takes
   approximately 55 minutes.

The total track length is roughly 41 km; the complete observation session
lasts about 4.7 hours, yielding approximately 280 shots per transponder
(1120 observations total).

### Sound-Speed Profile

A realistic Western-Pacific Munk-like profile is embedded in
`simulate_gnssa()`:

| Depth (m) | Speed (m/s) |
|-----------|-------------|
| 0         | 1521        |
| 900       | 1476 (SOFAR minimum) |
| 5000      | 1529        |

### Observation Model

Each synthetic two-way travel time is generated as:

```math
T^{\rm obs}_{n,k} = 2\,T^{\rm exact}_{n,k} + \frac{\rm NTD}(t_n)}{M(\xi_{n,k})} + \varepsilon_n
```

where:
- ``T^{\rm exact}_{n,k}`` is the one-way exact travel time from `xyz2tt()`
- ``{\rm NTD}(t)`` is a sinusoidal Nadir Total Delay model:
  ``5 \times 10^{-3} \sin(2\pi t / 3\,{\rm h}) + 10^{-3} \sin(2\pi t / 30\,{\rm min} + 1.5)``
- ``M(\xi_{n,k}) = 1/\cos(\xi_{n,k})`` is the mapping function
- ``\varepsilon_n \sim \mathcal{N}(0,\,10^{-5}\ {\rm s})`` is Gaussian noise

## Running the Simulation

### From the Julia REPL

```julia
using SeaGap

# Change to an empty working directory first
cd("/path/to/work/directory")

# Step 1 – Generate all SeaGap input files and obsdata.inp
simulate_gnssa(20.0, 5000.0)

# Step 2 – Estimate the array displacement
# NPB = 50 gives a knot spacing of ~8 min, sufficient to resolve
# the 30-min short-term NTD component in the simulated observations.
static_array(20.0, [2.0], 50)
```

### From the Command Line

A ready-to-run script is provided at `simulation/run_simulation.jl`:

```sh
cd simulation
julia --project=.. run_simulation.jl
```

## Function Reference

```@docs
simulate_gnssa
```

## Generated Files

`simulate_gnssa()` writes the following files to the current directory:

| File | Description |
|------|-------------|
| `ss_prof.inp` | Western-Pacific sound-speed profile |
| `pxp-ini.inp` | Initial seafloor transponder positions (4 × [x, y, z]) |
| `tr-ant.inp` | GNSS antenna – transducer lever-arm offset |
| `site_info.txt` | Site information (name, lon, lat, depth, numk) |
| `gps.jxyhhpr` | GPS trajectory: time, x, y, z, heading, pitch, roll |
| `pxp-1.jttq` … `pxp-4.jttq` | Per-transponder travel-time time series |
| `obsdata.inp` | Formatted observation data (input to positioning functions) |

## Expected Results

After calling `static_array(20.0, [2.0], 50)`, the estimated array
displacement in `position.out` should be close to zero (the true
displacement is zero by construction).  Typical values:

```
$ cat position.out
5.00...e8  ~0.001  ~0.002  ~-0.009  ~0.004  ~0.004  ~0.004
```

The sub-centimetre residuals arise from the interaction between the
injected NTD and the B-spline basis.  Use `static_array_AICBIC()` to
search for the optimal number of bases:

```julia
static_array_AICBIC(20.0, [2.0]; NPB_min=20, NPB_max=80)
```

## Customising the Simulation

All parameters can be overridden via keyword arguments:

```julia
# Shallower site (3000 m), 3 transponders, faster ship
simulate_gnssa(25.0, 3000.0;
    numk       = 3,
    ship_speed = 4.0,
    dt_shot    = 30.0,
    noise_tt   = 5.0e-5,
)
```

See the `simulate_gnssa` docstring for the full list of keyword arguments.
