# Coordinate transformations and utilities

using LinearAlgebra
using ..Math: Quat, quat_apply

"""
    ecef_to_geocentric(r_ecef)

Convert ECEF position vector to geocentric spherical coordinates.

# Arguments
- `r_ecef::Vector{Float64}`: Position in ECEF frame [m]

# Returns
- `(r, λ, Ω)`: Tuple of (radius [m], geocentric latitude [rad], longitude [rad])
"""
function ecef_to_geocentric(r_ecef::Vector{Float64})
    x, y, z = r_ecef

    # Radius from Earth center
    r = norm(r_ecef)

    # Geocentric latitude: angle from equatorial plane
    λ = asin(z / r)

    # Geocentric longitude: angle in equatorial plane from X-axis
    Ω = atan(y, x)

    return (r, λ, Ω)
end

"""
    ned_to_eci(B_ned, λ, Ω, gmst)

Transform vector from geocentric NED frame to ECI frame.

# Arguments
- `B_ned::Vector{Float64}`: Vector in NED frame (North, East, Down)
- `λ::Float64`: Geocentric latitude [rad]
- `Ω::Float64`: Geocentric longitude [rad]
- `gmst::Float64`: Greenwich Mean Sidereal Time [rad]

# Returns
- `Vector{Float64}`: Vector in ECI frame
"""
function ned_to_eci(B_ned::Vector{Float64}, λ::Float64, Ω::Float64, gmst::Float64)
    N, E, D = B_ned

    # Rotation from NED to ECEF (geocentric)
    # NED: X=North, Y=East, Z=Down
    # ECEF: X through equator/prime meridian, Z through north pole
    sin_λ = sin(λ)
    cos_λ = cos(λ)
    sin_Ω = sin(Ω)
    cos_Ω = cos(Ω)

    # Transform NED to ECEF
    B_ecef = [
        -sin_λ * cos_Ω * N - sin_Ω * E - cos_λ * cos_Ω * D,
        -sin_λ * sin_Ω * N + cos_Ω * E - cos_λ * sin_Ω * D,
         cos_λ * N - sin_λ * D
    ]

    # Rotation from ECEF to ECI (about Z-axis by GMST angle)
    sin_gmst = sin(gmst)
    cos_gmst = cos(gmst)

    B_eci = [
        cos_gmst * B_ecef[1] - sin_gmst * B_ecef[2],
        sin_gmst * B_ecef[1] + cos_gmst * B_ecef[2],
        B_ecef[3]
    ]

    return B_eci
end

"""
    eci_to_ecef(r_eci, gmst)

Transform vector from ECI to ECEF frame.

# Arguments
- `r_eci::Vector{Float64}`: Vector in ECI frame
- `gmst::Float64`: Greenwich Mean Sidereal Time [rad]

# Returns
- `Vector{Float64}`: Vector in ECEF frame
"""
function eci_to_ecef(r_eci::Vector{Float64}, gmst::Float64)
    # Rotation from ECI to ECEF (about Z-axis by -GMST angle)
    # This is the inverse of ECEF to ECI rotation
    sin_gmst = sin(gmst)
    cos_gmst = cos(gmst)

    r_ecef = [
        cos_gmst * r_eci[1] + sin_gmst * r_eci[2],
        -sin_gmst * r_eci[1] + cos_gmst * r_eci[2],
        r_eci[3]
    ]

    return r_ecef
end

"""
    eci_to_body(r_eci, q)

Transform vector from ECI to body frame using quaternion multiplication.

# Arguments
- `v_eci::Vector{Float64}`: Vector in ECI frame
- `q_eci_to_body::Vector{Float64}`: Attitude quaternion from ECI to body (scalar-first: [q0, q1, q2, q3])

# Returns
- `Vector{Float64}`: Vector in body frame
"""
function eci_to_body(v_eci::Vector{Float64}, q_eci_to_body::Vector{Float64})
    # Convert vector array to Quat
    q = Quat(q_eci_to_body)

    # Apply active rotation: v_body = q ⊗ v ⊗ q*
    # This rotates the vector from ECI to body frame
    v_body = quat_apply(q, v_eci, false)

    return v_body
end

"""
    body_to_eci(v_body, q_eci_to_body)

Transform vector from body frame to ECI frame using quaternion multiplication.

# Arguments
- `v_eci::Vector{Float64}`: Vector in ECI frame
- `q_eci_to_body::Vector{Float64}`: Attitude quaternion from ECI to body (scalar-first: [q0, q1, q2, q3])

# Returns
- `Vector{Float64}`: Vector in body frame
"""
function body_to_eci(v_body::Vector{Float64}, q_eci_to_body::Vector{Float64})
    # Convert vector array to Quat
    q = Quat(q_eci_to_body)
    q = quat_conj(q)  # Invert quaternion to go from body to ECI

    # Apply active rotation: v_eci = q ⊗ v ⊗ q*
    # This rotates the vector from body to ECI frame
    v_eci = quat_apply(q, v_body, false)

    return v_eci
end

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

    # Compute position and velocity in perifocal frame
    p = coes.a * (1 - coes.e^2)  # semi-latus rectum

    # TODO: check conditions at beginning

    # Position in perifocal frame with components (P, Q, W) where P points to periapsis
    r_mag = p / (1 + coes.e * cos(coes.ν))
    r_pf = r_mag * [cos(coes.ν), sin(coes.ν), 0.0]

    # Velocity in perifocal frame
    v_pf = sqrt(μ_earth / p) * [-sin(coes.ν), e + cos(coes.ν), 0.0]

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