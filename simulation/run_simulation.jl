#!/usr/bin/env julia
# =============================================================================
# run_simulation.jl
#
# Demonstrate a complete GNSS-Acoustic workflow:
#   1. Simulate synthetic observation data for a Western-Pacific 5000 m site.
#   2. Run static array positioning with SeaGap.
#   3. Print the estimated array displacement.
#
# Usage (from the simulation/ directory):
#   julia --project=.. run_simulation.jl
#
# Or from a Julia REPL with the SeaGap environment active:
#   include("simulation/run_simulation.jl")
# =============================================================================

using SeaGap
import DelimitedFiles

# --------------------------------------------------------------------------- #
#  Simulation scenario                                                        #
# --------------------------------------------------------------------------- #
lat      = 20.0    # site latitude [°N], representative Western Pacific
dep      = 5000.0  # water depth [m]  →  also used as the circle radius
TR_DEPTH = 2.0     # transducer depth below sea surface [m]

# Number of B-spline bases for NTD modelling in static_array().
# The 6.3-hour observation window contains two NTD components:
#   - long-term (3-hour period) and short-term (30-minute period).
# NPB = 50 gives a knot spacing of ~8 minutes, adequate to resolve both.
NPB = 50

println("=" ^ 60)
println(" GNSS-A Simulation – Western Pacific, $(dep) m depth")
println("=" ^ 60)

# --------------------------------------------------------------------------- #
#  Step 1 – Generate synthetic observation data                              #
# --------------------------------------------------------------------------- #
println("\n[Step 1]  Generating synthetic observation data …")

simulate_gnssa(lat, dep;
    numk       = 4,         # 4 seafloor transponders
    TR_DEPTH   = TR_DEPTH,
    ant_height = 25.0,      # GNSS antenna 25 m above sea surface
    ship_speed = 3.0,       # vessel speed [m/s]
    dt_gps     = 1.0,       # GPS sampling rate [s]
    dt_shot    = 60.0,      # acoustic shot interval [s]
    noise_tt   = 1.0e-5,    # travel-time random noise [s]
    array_half = 700.0,     # transponder array half-width [m]
)

# --------------------------------------------------------------------------- #
#  Step 2 – Static array positioning                                         #
# --------------------------------------------------------------------------- #
println("\n[Step 2]  Static array positioning ($(NPB) B-spline bases) …")

static_array(lat, [TR_DEPTH], NPB)

# --------------------------------------------------------------------------- #
#  Step 3 – Report results                                                   #
# --------------------------------------------------------------------------- #
println("\n[Results]  Estimated array displacement (true value = 0 m):")

if isfile("position.out")
    d = DelimitedFiles.readdlm("position.out")
    @printf("  EW  = %+.4f ± %.4f m\n", d[1,2], d[1,5])
    @printf("  NS  = %+.4f ± %.4f m\n", d[1,3], d[1,6])
    @printf("  UD  = %+.4f ± %.4f m\n", d[1,4], d[1,7])
else
    println("  (position.out not found)")
end

println("\nOutput files: solve.out, position.out, residual.out, bspline.out, AICBIC.out")
