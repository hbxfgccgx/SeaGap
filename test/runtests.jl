using SeaGap
using Test

@testset "SeaGap.jl" begin
    # Write your tests here.
end

@testset "TEOS-10 Sound Speed" begin
    # Reference values from the GSW (Gibbs SeaWater) Oceanographic Toolbox
    # gsw_sound_speed(SA, CT, p) — Roquet et al. (2015)

    @testset "soundspeed_teos10" begin
        # Standard conditions: SA = 35 g/kg, CT = 10 °C, p = 0 dbar
        c = soundspeed_teos10(35.0, 10.0, 0.0)
        @test isapprox(c, 1489.643273, atol=1e-3)

        # Cold seawater at surface
        c = soundspeed_teos10(35.0, 0.0, 0.0)
        @test isapprox(c, 1448.782907, atol=1e-3)

        # Warm seawater at surface
        c = soundspeed_teos10(35.0, 20.0, 0.0)
        @test isapprox(c, 1521.214277, atol=1e-3)

        # Deep ocean (1000 dbar)
        c = soundspeed_teos10(35.0, 10.0, 1000.0)
        @test isapprox(c, 1506.402330, atol=1e-3)

        # Deep ocean (3000 dbar)
        c = soundspeed_teos10(35.0, 10.0, 3000.0)
        @test isapprox(c, 1540.509385, atol=1e-3)

        # Very deep ocean (6000 dbar)
        c = soundspeed_teos10(35.0, 10.0, 6000.0)
        @test isapprox(c, 1592.734062, atol=1e-3)

        # Fresh water
        c = soundspeed_teos10(0.0, 0.0, 0.0)
        @test isapprox(c, 1402.424920, atol=1e-3)

        # Typical deep-ocean conditions
        c = soundspeed_teos10(34.0, 2.0, 2000.0)
        @test isapprox(c, 1489.837390, atol=1e-3)
    end

    @testset "sp_to_sa" begin
        # Practical Salinity 35 PSU at surface → SA ≈ 35.165 g/kg
        SA = sp_to_sa(35.0, 0.0, 180.0, 0.0)
        @test isapprox(SA, 35.165284, atol=1e-3)

        # Fresh water remains fresh
        SA = sp_to_sa(0.0, 0.0, 0.0, 0.0)
        @test isapprox(SA, 0.0, atol=1e-6)
    end

    @testset "t_to_ct" begin
        # At the surface, CT ≈ t
        CT = t_to_ct(35.0, 10.0, 0.0)
        @test isapprox(CT, 9.992855, atol=1e-3)

        # CT and t diverge at depth (adiabatic heating)
        CT = t_to_ct(35.0, 10.0, 1000.0)
        @test CT < 10.0  # CT < t due to TEOS-10 potential enthalpy definition
    end

    @testset "soundspeed_teos10_from_ctd" begin
        # CTD data: SP=35 PSU, t=10 °C, p=0 dbar
        c = soundspeed_teos10_from_ctd(35.0, 10.0, 0.0)
        @test isapprox(c, 1489.587941, atol=0.5)  # tolerance allows for SA conversion

        # Sound speed increases with pressure
        c1 = soundspeed_teos10_from_ctd(35.0, 10.0, 0.0)
        c2 = soundspeed_teos10_from_ctd(35.0, 10.0, 3000.0)
        @test c2 > c1

        # Sound speed increases with temperature (at surface)
        c_cold = soundspeed_teos10_from_ctd(35.0, 5.0, 0.0)
        c_warm = soundspeed_teos10_from_ctd(35.0, 20.0, 0.0)
        @test c_warm > c_cold

        # Sound speed increases with salinity (at surface)
        c_fresh = soundspeed_teos10_from_ctd(0.0, 10.0, 0.0)
        c_salt  = soundspeed_teos10_from_ctd(35.0, 10.0, 0.0)
        @test c_salt > c_fresh
    end
end
