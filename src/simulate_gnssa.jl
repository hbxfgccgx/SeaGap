#=
  simulate_gnssa: Simulate GNSS-Acoustic observation data for SeaGap.

  Generates a synthetic dataset for a seafloor transponder array in a
  Western-Pacific scenario at a user-supplied water depth.  The sea-surface
  vessel track is composed of:
    • A circle of radius equal to the water depth (1× dep) centred over
      the transponder array.
    • A cross pattern – an EW leg and a NS leg – that passes directly over
      the centre of the array ("过顶十字").

  All SeaGap-compatible input files are written to the current working
  directory.  After calling this function, run obsdata_format (called
  automatically) and then any positioning function (e.g. static_array).

  written 2024
=#

export simulate_gnssa

"""
    simulate_gnssa(lat, dep; numk, TR_DEPTH, ant_height, ship_speed,
                   dt_gps, dt_shot, noise_tt, array_half, t0, sitename,
                   fno_ss, fno_pxp, fno_ant, fno_site, fno_obs,
                   fn_prefix, fn_suffix)

Simulate a synthetic GNSS-Acoustic dataset for a seafloor transponder array
and write all SeaGap-compatible input files to the current working directory.

The simulation environment is representative of the Western Pacific Ocean
at the given water depth `dep`.  The sea-surface vessel track consists of:

* A **circle** of radius `dep` metres (1 × water depth) centred over the
  transponder array.
* A **cross** pattern – an EW leg and a NS leg – that passes directly over
  the centre of the array.

Travel times are computed by exact ray-tracing (`xyz2tt`) through a
realistic Western-Pacific sound-speed profile.  A sinusoidal Nadir Total
Delay (NTD) and Gaussian random noise are added to each observation.

After the function returns, call `static_array` (or any other SeaGap
positioning function) in the same directory to estimate the array position.

# Arguments
* `lat` : Site latitude [°N] (default `20.0`)
* `dep` : Water depth [m]; also used as the circle radius (default `5000.0`)

# Keyword Arguments
* `numk`       : Number of seafloor transponders, 1–4 (default `4`)
* `TR_DEPTH`   : Transducer depth below sea surface [m] (default `2.0`)
* `ant_height` : GNSS antenna height above sea surface [m] (default `25.0`)
* `ship_speed` : Vessel speed [m s⁻¹] (default `3.0`)
* `dt_gps`     : GPS sampling interval [s] (default `1.0`)
* `dt_shot`    : Acoustic shot interval [s] (default `60.0`)
* `noise_tt`   : One-sigma random travel-time noise [s] (default `1.0e-5`)
* `array_half` : Half-width of the square transponder array [m] (default `700.0`)
* `t0`         : Reference start time [s from epoch] (default `5.0e8`)
* `sitename`   : Site identifier written to `site_info.txt` (default `"SIM1"`)
* `fno_ss`     : Output path for the sound-speed profile (default `"ss_prof.inp"`)
* `fno_pxp`    : Output path for transponder initial positions (default `"pxp-ini.inp"`)
* `fno_ant`    : Output path for antenna–transducer offset (default `"tr-ant.inp"`)
* `fno_site`   : Output path for site information (default `"site_info.txt"`)
* `fno_obs`    : Output path for the formatted observation data (default `"obsdata.inp"`)
* `fn_prefix`  : Filename prefix for per-transponder travel-time files (default `"pxp-"`)
* `fn_suffix`  : Filename suffix for per-transponder travel-time files (default `".jttq"`)

# Example
    # Run from an empty working directory:
    simulate_gnssa(20.0, 5000.0)
    static_array(20.0, [2.0], 20)
"""
function simulate_gnssa(lat=20.0, dep=5000.0;
        numk::Int64  = 4,
        TR_DEPTH     = 2.0,
        ant_height   = 25.0,
        ship_speed   = 3.0,
        dt_gps       = 1.0,
        dt_shot      = 60.0,
        noise_tt     = 1.0e-5,
        array_half   = 700.0,
        t0           = 5.0e8,
        sitename     = "SIM1",
        fno_ss       = "ss_prof.inp",
        fno_pxp      = "pxp-ini.inp",
        fno_ant      = "tr-ant.inp",
        fno_site     = "site_info.txt",
        fno_obs      = "obsdata.inp",
        fn_prefix    = "pxp-",
        fn_suffix    = ".jttq")

    println(stderr, " === GNSS-A simulation: simulate_gnssa ===")

    if numk < 1 || numk > 4
        error("simulate_gnssa: numk must be between 1 and 4")
    end
    if dep <= 0.0
        error("simulate_gnssa: dep must be positive")
    end

    Random.seed!(42)   # reproducible noise

    # ------------------------------------------------------------------ #
    #  Step 1 – Write the Western-Pacific sound-speed profile             #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Writing $fno_ss")

    # Typical Munk-like profile with SOFAR minimum near 900 m.
    # The last depth entry is dep + 300 m to guarantee the profile covers
    # the transponder depth for any dep ≤ 5000 m.
    ss_z = [0.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0,
            600.0, 700.0, 800.0, 900.0, 1000.0, 1200.0, 1500.0,
            2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0,
            5000.0, 5300.0]
    ss_v = [1521.0, 1520.0, 1519.0, 1515.0, 1509.0, 1499.0, 1491.0,
            1485.0, 1481.0, 1479.0, 1477.0, 1476.0, 1476.0, 1477.0,
            1480.0, 1486.0, 1494.0, 1502.0, 1509.0, 1515.0, 1520.0,
            1525.0, 1529.0, 1531.0]

    # For water depths > 5000 m, extend the profile linearly.
    if dep > 5000.0
        grad = (ss_v[end] - ss_v[end-1]) / (ss_z[end] - ss_z[end-1])
        old_end_z = ss_z[end]
        push!(ss_z, dep + 300.0)
        push!(ss_v, ss_v[end] + grad * ((dep + 300.0) - old_end_z))
    end

    open(fno_ss, "w") do f
        for (d, s) in zip(ss_z, ss_v)
            @printf(f, "%.1f %.4f\n", d, s)
        end
    end

    # ------------------------------------------------------------------ #
    #  Step 2 – Write the transponder initial positions (square array)   #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Writing $fno_pxp")

    # Square array centred at (0, 0, −dep).  Only the first numk are used.
    all_px = [ array_half, -array_half, -array_half,  array_half]
    all_py = [ array_half,  array_half, -array_half, -array_half]
    pxp_x  = all_px[1:numk]
    pxp_y  = all_py[1:numk]
    pxp_z  = fill(-dep, numk)

    open(fno_pxp, "w") do f
        for k in 1:numk
            @printf(f, "%.6f %.6f %.6f\n", pxp_x[k], pxp_y[k], pxp_z[k])
        end
    end

    # ------------------------------------------------------------------ #
    #  Step 3 – Write the antenna–transducer lever-arm offset            #
    # ------------------------------------------------------------------ #
    # The transducer is located directly below the GNSS antenna.
    # ant_height [m above sea surface] + TR_DEPTH [m below sea surface]
    # gives the total vertical offset.
    println(stderr, " --- Writing $fno_ant")
    ez = -(ant_height + TR_DEPTH)   # e.g. -(25 + 2) = -27 m
    open(fno_ant, "w") do f
        @printf(f, "0.0 0.0 %.4f\n", ez)
    end

    # ------------------------------------------------------------------ #
    #  Step 4 – Write the site information file                          #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Writing $fno_site")
    lon_site = 143.0   # representative Western Pacific longitude
    open(fno_site, "w") do f
        @printf(f, "%s %.6f %.6f %.1f %d\n", sitename, lon_site, lat, dep, numk)
    end

    # ------------------------------------------------------------------ #
    #  Step 5 – Load the sound-speed profile for travel-time computation #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Loading sound-speed profile for ray-tracing")
    z_prof, v_prof, nz_st, numz = read_prof(fno_ss, TR_DEPTH)
    Rg, _  = localradius(lat)
    # Transducer z-height [m above mean sea level]:
    # GNSS antenna at +ant_height, transducer offset ez below → zd_tr < 0
    zd_tr  = Float64(ant_height) + ez   # = -TR_DEPTH

    # ------------------------------------------------------------------ #
    #  Step 6 – Build the vessel track (circle + EW cross + NS cross)   #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Generating vessel track (circle + EW + NS legs)")
    radius = dep   # circle radius = 1 × water depth

    # 6a. Circle – counter-clockwise, 72 waypoints (5° spacing), closed.
    n_circ  = 72
    theta_c = [2π * i / n_circ for i in 0:n_circ]   # n_circ + 1 to close
    circ_x  = radius .* cos.(theta_c)
    circ_y  = radius .* sin.(theta_c)

    # 6b. EW cross leg: west edge → east edge (passes directly over centre)
    n_cross = 21   # 20 equal segments → 21 points
    ew_x    = collect(range(-radius, radius, length=n_cross))
    ew_y    = zeros(n_cross)

    # 6c. NS cross leg: south edge → north edge (passes directly over centre)
    ns_x    = zeros(n_cross)
    ns_y    = collect(range(-radius, radius, length=n_cross))

    # Concatenate all waypoints.  The gaps between the circle end and the
    # start of the EW leg, and between the EW end and the NS start, are
    # traversed at ship_speed as straight-line transits.
    all_x = vcat(circ_x, ew_x, ns_x)
    all_y = vcat(circ_y, ew_y, ns_y)

    # ------------------------------------------------------------------ #
    #  Step 7 – Interpolate the track to a uniform GPS time series       #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Interpolating track to GPS time series (dt=$(dt_gps) s)")

    # Cumulative arc length along all waypoints.
    n_wp  = length(all_x)
    arc   = zeros(n_wp)
    for i in 2:n_wp
        arc[i] = arc[i-1] + hypot(all_x[i] - all_x[i-1], all_y[i] - all_y[i-1])
    end
    total_length = arc[end]
    total_time   = total_length / ship_speed

    # Interpolate (x, y) at arc-length position s.
    function pos_at_s(s)
        s_cl = clamp(s, 0.0, arc[end])
        idx  = clamp(searchsortedlast(arc, s_cl), 1, n_wp - 1)
        denom = arc[idx+1] - arc[idx]
        frac  = denom > 0.0 ? (s_cl - arc[idx]) / denom : 0.0
        x = all_x[idx] + frac * (all_x[idx+1] - all_x[idx])
        y = all_y[idx] + frac * (all_y[idx+1] - all_y[idx])
        return x, y
    end

    # GPS sample times.  Add a 30-second buffer so that receive times
    # after the final shot are covered.
    n_gps  = round(Int, (total_time + 30.0) / dt_gps) + 1
    gps_t  = [t0 + i * dt_gps for i in 0:n_gps-1]
    gps_x  = zeros(n_gps)
    gps_y  = zeros(n_gps)
    gps_z  = fill(ant_height, n_gps)  # GNSS antenna height [m above sea level]
    gps_h  = zeros(n_gps)   # heading [deg, clockwise from North]
    gps_p  = zeros(n_gps)   # pitch   [deg]
    gps_r  = zeros(n_gps)   # roll    [deg]

    for i in 1:n_gps
        s = min((i - 1) * dt_gps * ship_speed, arc[end])
        gps_x[i], gps_y[i] = pos_at_s(s)
    end

    # Heading from successive positions (compass convention).
    for i in 1:n_gps-1
        dx = gps_x[i+1] - gps_x[i]
        dy = gps_y[i+1] - gps_y[i]
        if hypot(dx, dy) > 0.01
            gps_h[i] = mod(atan(dx, dy) * 180.0 / π, 360.0)
        elseif i > 1
            gps_h[i] = gps_h[i-1]
        end
    end
    gps_h[n_gps] = gps_h[n_gps-1]

    # ------------------------------------------------------------------ #
    #  Step 8 – Write the GPS trajectory file (gps.jxyhhpr)             #
    # ------------------------------------------------------------------ #
    fn_gps = "gps.jxyhhpr"
    println(stderr, " --- Writing $fn_gps ($n_gps points)")
    open(fn_gps, "w") do f
        for i in 1:n_gps
            @printf(f, "%.2f %.6f %.6f %.6f %.5f %.5f %.5f\n",
                gps_t[i], gps_x[i], gps_y[i], gps_z[i],
                gps_h[i], gps_p[i], gps_r[i])
        end
    end

    # ------------------------------------------------------------------ #
    #  Step 9 – Compute synthetic travel times                          #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Computing synthetic travel times")

    # Linear interpolation of GPS position at arbitrary time t.
    function gps_pos_at_t(t)
        idx  = clamp(searchsortedlast(gps_t, t), 1, n_gps - 1)
        frac = (t - gps_t[idx]) / (gps_t[idx+1] - gps_t[idx])
        x = gps_x[idx] + frac * (gps_x[idx+1] - gps_x[idx])
        y = gps_y[idx] + frac * (gps_y[idx+1] - gps_y[idx])
        return x, y
    end

    # NTD model: long-period + short-period sinusoidal variation.
    # NTD(t) is the nadir-direction travel-time anomaly [s].
    # The contribution to the slant-path two-way TT is NTD / vert_avg,
    # consistent with the normalisation used in static_array.
    A_long  = 5.0e-3;  T_long  = 3.0 * 3600.0   # 5 ms, 3-hour period
    A_short = 1.0e-3;  T_short = 30.0 * 60.0    # 1 ms, 30-min period
    ntd_func(t) = (A_long  * sin(2π * (t - t0) / T_long) +
                   A_short * sin(2π * (t - t0) / T_short + 1.5))

    # Shot times: every dt_shot seconds, within the GPS coverage window.
    t_shot_end = gps_t[end] - 20.0
    shot_times = collect(t0 + dt_shot : dt_shot : t_shot_end)
    n_shots    = length(shot_times)

    # Per-transponder storage: (shot_time, tt_with_delay, tt, quality)
    jttq = [Vector{NTuple{4,Float64}}() for _ in 1:numk]

    for ts in shot_times
        xs1, ys1 = gps_pos_at_t(ts)
        ntd     = ntd_func(ts)
        for k in 1:numk
            # Forward one-way exact travel time at the transmit position.
            tc_ow1, _, vert1 = xyz2tt(pxp_x[k], pxp_y[k], pxp_z[k],
                                       xs1, ys1, zd_tr,
                                       z_prof, v_prof, nz_st, numz,
                                       Rg, TR_DEPTH)
            # Receive time (approximate): t_transmit + 2 × one-way TT.
            t_receive  = ts + 2.0 * tc_ow1
            xs2, ys2   = gps_pos_at_t(t_receive)
            # Backward one-way exact travel time at the receive position.
            tc_ow2, _, vert2 = xyz2tt(pxp_x[k], pxp_y[k], pxp_z[k],
                                       xs2, ys2, zd_tr,
                                       z_prof, v_prof, nz_st, numz,
                                       Rg, TR_DEPTH)
            # Two-way travel time = forward leg + backward leg.
            tc_two_way = tc_ow1 + tc_ow2
            # Average mapping factor (consistent with static_array).
            vert_avg   = (vert1 + vert2) / 2.0
            # Observed TT = geometric two-way TT + NTD contribution + noise.
            # NTD / vert_avg converts nadir-direction anomaly to slant path.
            tt_obs = tc_two_way + ntd / vert_avg + randn() * noise_tt
            push!(jttq[k], (ts, tt_obs, tt_obs, 1.0))
        end
    end

    # Write one file per transponder.
    for k in 1:numk
        fn = "$(fn_prefix)$(k)$(fn_suffix)"
        println(stderr, " --- Writing $fn ($(length(jttq[k])) shots)")
        open(fn, "w") do f
            for (ts, tt0_v, tt_v, q) in jttq[k]
                @printf(f, "%.6f %.6f %.6f %d\n", ts, tt0_v, tt_v, Int(q))
            end
        end
    end

    # ------------------------------------------------------------------ #
    #  Step 10 – Format all data into obsdata.inp                        #
    # ------------------------------------------------------------------ #
    println(stderr, " --- Formatting observation data to $fno_obs")
    obsdata_format(numk, 1;
        fno  = fno_obs,
        fn1  = fn_gps,
        fn21 = fn_prefix,
        fn22 = fn_suffix)

    # ------------------------------------------------------------------ #
    #  Summary                                                           #
    # ------------------------------------------------------------------ #
    println(stderr, " === Simulation complete ===")
    @printf(stderr, "     GPS points     : %d  (%.2f hours)\n", n_gps, total_time / 3600)
    @printf(stderr, "     Shots per PXP  : %d\n", n_shots)
    @printf(stderr, "     Transponders   : %d\n", numk)
    @printf(stderr, "     Total obs      : %d\n", n_shots * numk)
    println(stderr, " Files written:")
    println(stderr, "   $fno_ss")
    println(stderr, "   $fno_pxp")
    println(stderr, "   $fno_ant")
    println(stderr, "   $fno_site")
    println(stderr, "   $fn_gps")
    for k in 1:numk
        println(stderr, "   $(fn_prefix)$(k)$(fn_suffix)")
    end
    println(stderr, "   $fno_obs")
    return (fno_ss=fno_ss, fno_pxp=fno_pxp, fno_ant=fno_ant,
            fno_site=fno_site, fn_gps=fn_gps, fno_obs=fno_obs)
end
