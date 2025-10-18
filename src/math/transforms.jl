# Coordinate transformations and utilities

using LinearAlgebra
using Simwise.Math: Quat, quat_apply

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