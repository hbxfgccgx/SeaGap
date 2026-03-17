#=
  Functions for calculating seawater sound speed using TEOS-10
    (International Thermodynamic Equation of Seawater - 2010)

  Theory
  ------
  The speed of sound in seawater is derived from the thermodynamic relation:

      c² = -(∂p/∂ρ)_η = v² (∂p/∂v)_η

  where:
    c   : speed of sound [m/s]
    p   : sea pressure [Pa]
    ρ   : in-situ density [kg/m³]
    v   : specific volume = 1/ρ [m³/kg]
    η   : specific entropy [J/(kg·K)]
    _η  : derivative at constant entropy (isentropic)

  Using the seawater Gibbs function g(S_A, T, p) from TEOS-10, the speed of
  sound is expressed as (IOC et al. 2010, eq. 2.17.1):

      c² = -g_p² / (g_pp - g_Tp² / g_TT)

  where subscripts denote partial derivatives:
    g_p  = (∂g/∂p)_{S_A,T}  = v                  (specific volume)
    g_pp = (∂²g/∂p²)_{S_A,T} = -v κ_T            (isothermal compressibility)
    g_Tp = (∂²g/∂T∂p)_{S_A} = v α_T              (thermal expansion)
    g_TT = (∂²g/∂T²)_{S_A,p} = -C_p / T          (isobaric heat capacity)

  This simplifies to (using the 75-term polynomial, Roquet et al. 2015):

      c² = -v² / (v_P + v_{CT}² / C_{p0})

  where:
    v_P   = (∂v/∂p)_{S_A,CT}  : isothermal pressure derivative of specific volume
    v_{CT} = (∂v/∂CT)_{S_A,p} : specific volume derivative w.r.t. Conservative
                                  Temperature (CT)
    C_{p0} = 3991.86795711963 J/(kg·K) : reference specific heat capacity
                                           (isobaric, at fixed SA and CT)

  The 75-term polynomial for specific volume (Roquet et al. 2015) is expressed
  as a function of normalised Absolute Salinity x_S, Conservative Temperature
  y_T, and pressure z:

      x_S = √(S_A / S_{SO})  where S_{SO} = 35.16504 g/kg
      y_T = CT / 40
      z   = p / 10000

  References
  ----------
  - IOC, SCOR and IAPSO (2010): The international thermodynamic equation of
    seawater - 2010: Calculation and use of thermodynamic properties.
    Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
    UNESCO (English), 196 pp.
  - McDougall, T.J. and Barker, P.M. (2011): Getting started with TEOS-10 and
    the Gibbs Seawater (GSW) Oceanographic Toolbox. 28 pp., SCOR/IAPSO WG127.
  - Roquet, F., G. Madec, T.J. McDougall, P.M. Barker (2015): Accurate
    polynomial expressions for the density and specific volume of seawater
    using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
  - McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz, P.M. Barker
    (2012): A global algorithm for estimating Absolute Salinity. Ocean Science,
    8(6), pp. 1123-1134.
  written by SeaGap contributors 2024
=#

import GibbsSeaWater as GSW

export soundspeed_teos10, sp_to_sa, t_to_ct, soundspeed_teos10_from_ctd

# === Convert Practical Salinity to Absolute Salinity
"""
    sp_to_sa(SP, p, lon, lat)

Convert Practical Salinity `SP` (PSU) to Absolute Salinity `SA` (g/kg) using
the TEOS-10 global salinity anomaly atlas (McDougall et al. 2012).

Absolute Salinity S_A accounts for the spatially varying composition of
seawater and is related to Practical Salinity S_P by:

    S_A = (35.16504 / 35) × S_P + δS_A(lon, lat, p)

where δS_A is the salinity anomaly from the reference composition. At the
standard ocean composition (open ocean), δS_A ≈ 0, giving:

    S_A ≈ (35.16504 / 35) × S_P  ≈ 1.004715 × S_P

Arguments:
* `SP`  : Practical Salinity [PSU]
* `p`   : Sea pressure (absolute pressure - 10.1325 dbar) [dbar]
* `lon` : Longitude [°E], used for the salinity anomaly atlas
* `lat` : Latitude [°N], used for the salinity anomaly atlas

Output:
* `SA`  : Absolute Salinity [g/kg]

# Example
    SA = sp_to_sa(35.0, 0.0, 180.0, 0.0)  # SA ≈ 35.165 g/kg
"""
function sp_to_sa(SP, p, lon, lat)
  return GSW.gsw_sa_from_sp(SP, p, lon, lat)
end

# === Convert in-situ temperature to Conservative Temperature
"""
    t_to_ct(SA, t, p)

Convert in-situ temperature `t` (°C, ITS-90) to Conservative Temperature `CT`
(°C) using TEOS-10 (McDougall & Barker 2011).

Conservative Temperature Θ is proportional to potential enthalpy h₀:

    Θ = h₀(S_A, θ) / C_{p0}

where:
  - θ   : potential temperature referenced to the surface [°C]
  - h₀  : potential enthalpy referenced to p = 0 [J/kg]
  - C_{p0} = 3991.86795711963 J/(kg·K) : reference specific heat capacity

At the sea surface (p = 0), CT ≈ t (Conservative Temperature ≈ in-situ
temperature). With increasing pressure the two diverge due to adiabatic
heating.

Arguments:
* `SA` : Absolute Salinity [g/kg]
* `t`  : In-situ temperature (ITS-90) [°C]
* `p`  : Sea pressure [dbar]

Output:
* `CT` : Conservative Temperature [°C]

# Example
    CT = t_to_ct(35.0, 10.0, 0.0)   # CT ≈ 9.993 °C at the surface
    CT = t_to_ct(35.0, 10.0, 1000.0) # CT ≈ 9.868 °C at 1000 dbar
"""
function t_to_ct(SA, t, p)
  return GSW.gsw_ct_from_t(SA, t, p)
end

# === TEOS-10 sound speed from Absolute Salinity and Conservative Temperature
"""
    soundspeed_teos10(SA, CT, p)

Calculate the speed of sound in seawater from Absolute Salinity `SA`,
Conservative Temperature `CT`, and sea pressure `p` using the TEOS-10
computationally efficient 75-term polynomial (Roquet et al. 2015).

## Theory

The speed of sound follows from the isentropic compressibility κ_S:

    c = 1 / √(ρ κ_S)

where ρ = 1/v is the in-situ density and κ_S (Pa⁻¹) is:

    κ_S = κ_T - T α_T² / (ρ C_p)

with:
  - κ_T = -(1/v)(∂v/∂p)_T     : isothermal compressibility [Pa⁻¹]
  - α_T = (1/v)(∂v/∂T)_p      : thermal expansion coefficient [K⁻¹]
  - C_p                         : isobaric heat capacity [J/(kg·K)]
  - T = t + 273.15              : absolute in-situ temperature [K]

In terms of the TEOS-10 Gibbs function g(S_A, T, p), this is equivalent to
(IOC et al. 2010, eq. 2.17.1):

    c² = -g_p² / (g_pp - g_Tp² / g_TT)

where g_p, g_pp, g_Tp, g_TT are first and second partial derivatives of the
Gibbs potential with respect to pressure p and in-situ temperature T.

Using the 75-term specific volume polynomial v(S_A, CT, p) of Roquet et al.
(2015), the sound speed can be computed efficiently as:

    c² = -v² / ( v_P + v_{CT}² / C_{p0} )

where:
  - v       = specific volume [m³/kg]
  - v_P     = (∂v/∂p)_{S_A,CT} [m³/(kg·Pa)]
  - v_{CT}  = (∂v/∂CT)_{S_A,p} [m³/(kg·K)]
  - C_{p0}  = 3991.86795711963 J/(kg·K) : reference isobaric heat capacity

The polynomial uses the normalised variables:

    x_S = √(S_A / S_{SO}),  y_T = CT / 40,  z = p / 10000

where S_{SO} = 35.16504 g/kg is the Standard Ocean Salinity.

Arguments:
* `SA` : Absolute Salinity [g/kg]  (range: 0 – 42 g/kg)
* `CT` : Conservative Temperature [°C]  (range: -2 – 40 °C)
* `p`  : Sea pressure [dbar]  (range: 0 – 10000 dbar)

Output:
* `c`  : Speed of sound in seawater [m/s]

# Example
    c = soundspeed_teos10(35.0, 10.0, 0.0)    # c ≈ 1489.6 m/s
    c = soundspeed_teos10(35.0, 10.0, 1000.0) # c ≈ 1506.4 m/s
    c = soundspeed_teos10(35.0,  2.0, 3000.0) # c ≈ 1520.6 m/s
"""
function soundspeed_teos10(SA, CT, p)
  return GSW.gsw_sound_speed(SA, CT, p)
end

# === Convenience: sound speed from CTD measurements
"""
    soundspeed_teos10_from_ctd(SP, t, p; lon=0.0, lat=0.0)

Calculate the speed of sound in seawater directly from CTD (Conductivity-
Temperature-Depth) measurements, converting to TEOS-10 variables internally.

## Conversion procedure

1. **Practical Salinity → Absolute Salinity** (McDougall et al. 2012):

       S_A = (35.16504 / 35) × S_P + δS_A(lon, lat, p)

   where δS_A is the spatially varying salinity anomaly from the TEOS-10 atlas.

2. **In-situ temperature → Conservative Temperature** (TEOS-10):

       Θ = h₀(S_A, θ) / C_{p0}

   where h₀ is potential enthalpy at the surface and θ is potential temperature.

3. **Sound speed from (S_A, Θ, p)** using the TEOS-10 75-term polynomial:

       c² = -v² / ( v_P + v_{Θ}² / C_{p0} )

Arguments:
* `SP`  : Practical Salinity (PSU)
* `t`   : In-situ temperature (ITS-90) [°C]
* `p`   : Sea pressure [dbar]

Keyword arguments:
* `lon` : Longitude [°E], default 0.0 (used for salinity anomaly correction)
* `lat` : Latitude [°N], default 0.0 (used for salinity anomaly correction)

Output:
* `c`   : Speed of sound in seawater [m/s]

# Example
    # Typical open-ocean surface conditions
    c = soundspeed_teos10_from_ctd(35.0, 10.0, 0.0)
    # With geographic correction for precise SA
    c = soundspeed_teos10_from_ctd(34.5, 2.0, 3000.0; lon=180.0, lat=30.0)
"""
function soundspeed_teos10_from_ctd(SP, t, p; lon=0.0, lat=0.0)
  SA = GSW.gsw_sa_from_sp(SP, p, lon, lat)
  CT = GSW.gsw_ct_from_t(SA, t, p)
  return GSW.gsw_sound_speed(SA, CT, p)
end
