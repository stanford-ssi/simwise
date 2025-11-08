using LinearAlgebra

using Simwise.Constants: seconds_per_day

##################################################
#       Classical Orbital Elements
##################################################

mutable struct OrbitalElements
    a::Float64 #  Semi-major axis [km]
    e::Float64 #  Eccentricity
    i::Float64 #  Inclination [rad]
    Ω::Float64 #  Right ascension of ascending node (RAAN) [rad]
    ω::Float64 #  Argument of periapsis [rad]
    ν::Float64 #  True anomaly [rad]
end


"""
    rv_to_orbital_elements(r_eci, v_eci, forces)

Convert ECI position and velocity vectors to classical Keplerian orbital elements.

# Arguments
- `r_eci::Vector{Float64}`: Position of satellite [km] (ECI)
- `v_eci::Vector{Float64}`: Velocity of satellite [km/s] (ECI)
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `OrbitalElements`: Orbital Elements struct
"""
function rv_to_orbital_elements(r_eci::Vector{Float64}, v_eci::Vector{Float64}, μ::Float64)
    # Specific relative angular momentum
    h_eci = cross(r_eci,v_eci) 
    h_mag = norm(h_eci)

    # Vector pointing to ascending node
    n_eci = cross([0.0,0.0,1.0], h_eci) 
    n_mag = norm(n_eci)
    equatorial = isapprox(n_mag, 0.0, atol=1e-8)
    
    # Norms to make this easier
    r_mag = norm(r_eci)
    v_mag = norm(v_eci)
    
    # Eccentricity
    e_eci = ((v_mag^2 - μ/r_mag) * r_eci - dot(r_eci, v_eci) * v_eci ) / μ
    e_mag = norm(e_eci)
    circular = isapprox(e_mag, 0.0, atol=1e-8)
    
    # Specific Orbital Energy [km^2/s^2]
    ϵ = (v_mag^2)/2 - μ/r_mag
    
    # Semi-major Axis [km]
    a = isapprox(e_mag, 1.0, atol=1e-5) ? Inf : -μ/(2ϵ) 

    # Inclination [rad]
    i = acos(h_eci[3]/h_mag)

    # Right ascension of the ascending node [rad]
    if !equatorial
        Ω = acos(clamp( n_eci[1]/n_mag, -1, 1))
        if n_eci[2] < 0
            Ω = 2*pi - Ω
        end
    else
        Ω = 0.0
    end

    # Argument of periapse
    if !equatorial && !circular
        ω = acos(clamp( dot(n_eci,e_eci)/e_mag/n_mag, -1, 1))
        if e_eci[3] < 0
            ω = 2*pi - ω
        end
    else
        ω = 0.0
    end


    # True anomaly
    if !circular
        ν = acos(clamp( dot(r_eci,e_eci)/e_mag/r_mag, -1, 1))
        if dot(r_eci, v_eci) < 0
            ν = 2*pi - ν
        end
    else
        ν = 0.0
    end

    return OrbitalElements(a, e_mag, i, Ω, ω, ν)

end

"""
    orbital_elements_to_eci(coes)

Convert classical Keplerian orbital elements to ECI position and velocity vectors.

# Arguments
- `coes::OrbitalElements`: Classical orbital elements struct containing
    - `a::Float64`: Semi-major axis [m]
    - `e::Float64`: Eccentricity
    - `i::Float64`: Inclination [rad]
    - `Ω::Float64`: Right ascension of ascending node (RAAN) [rad]
    - `ω::Float64`: Argument of periapsis [rad]
    - `ν::Float64`: True anomaly [rad]
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `(r_eci, v_eci)`: Tuple of position [km] and velocity [km/s] vectors in ECI frame

"""
function orbital_elements_to_rv(coes::OrbitalElements, μ::Float64)

    # TODO: There should be exceptions for Circular Equatorial, Circular Inclined, and Elliptical Equatorial
    # but those require the r vector so we could?? iterate? but our state should be tracked in r,v anyways
    # so we'd realistically only be using the rv_to_coes function

    # Compute position and velocity in perifocal frame
    p = coes.a * (1 - coes.e^2)  # semi-latus rectum

    # TODO: check conditions at beginning

    # Position in perifocal frame with components (P, Q, W) where P points to periapsis
    r_mag = p / (1 + coes.e * cos(coes.ν))
    r_pf = r_mag * [cos(coes.ν), sin(coes.ν), 0.0]

    # Velocity in perifocal frame
    v_pf = sqrt(μ / p) * [-sin(coes.ν), coes.e + cos(coes.ν), 0.0]

    # Rotation matrix from perifocal to ECI
    # R = R3(-Ω) * R1(-i) * R3(-ω)
    sin_Ω = sin(coes.Ω)
    cos_Ω = cos(coes.Ω)
    sin_i = sin(coes.i)
    cos_i = cos(coes.i)
    sin_ω = sin(coes.ω)
    cos_ω = cos(coes.ω)

    # Combined rotation matrix (perifocal to ECI)
    R = [
        cos_Ω*cos_ω - sin_Ω*sin_ω*cos_i    -cos_Ω*sin_ω - sin_Ω*cos_ω*cos_i     sin_Ω*sin_i;
        sin_Ω*cos_ω + cos_Ω*sin_ω*cos_i    -sin_Ω*sin_ω + cos_Ω*cos_ω*cos_i    -cos_Ω*sin_i;
        sin_ω*sin_i                          cos_ω*sin_i                          cos_i
    ]

    # Transform to ECI
    r_eci = R * r_pf
    v_eci = R * v_pf

    return (r_eci, v_eci)
end

##################################################
#       Two Line Elements (TLE)
##################################################

"""
    sma_to_orbit_period(a_km)

Convert semi-major axis to orbit period

# Arguments
- `a::Float64`: Semi-major axis of orbit [km]
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `Float64`: Orbit Period [s]
"""
function sma_to_orbit_period(a_km::Float64, μ::Float64)
    return 2*pi*sqrt(a_km^3 / μ)
end

"""
    orbit_period_to_sma(T_sec)

Convert orbit period to semi-major axis

# Arguments
- `T_sec::Float64`: Orbit period [s]
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `Float64`: Semi-major axis of orbit [km]
"""
function orbit_period_to_sma(T_sec::Float64, μ::Float64)
    return ( (T_sec/2/pi)^2 * μ) ^ (1/3)
end

"""
    mean_motion_to_sma(n_rev_per_day, μ)

Converts mean motion into semi-major axis
Useful function for two-line-elements (TLE)

# Arguments
- `n_rev_per_day::Float64`: Mean Motion [rev/day]
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `Float64`: Semi-major axis of orbit [km]
"""
function mean_motion_to_sma(n_rev_per_day::Float64, μ::Float64)
    n_rps = n_rev_per_day / (seconds_per_day) * 2*pi # mean motion in rad/seconds

    return (μ*1000/n_rps)^(1/3)
end