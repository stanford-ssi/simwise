# Orbital dynamics

using LinearAlgebra

"""
    propagate_keplerian(state, dt)

Analytically propagate orbital elements forward in time using two-body dynamics.
Assumes no perturbations (ideal Keplerian orbit).

# Arguments
- `state::State`: Current state with orbital elements [a, e, i, Ω, ω, ν]
- `dt::Float64`: Time step [s]

# Returns
- `new_oe::Vector{Float64}`: Updated orbital elements [a, e, i, Ω, ω, ν]

# Notes
- Only the true anomaly ν changes in two-body motion
- All other elements (a, e, i, Ω, ω) remain constant
- Semi-major axis a is in meters
"""
function propagate_keplerian(state::State, dt::Float64)
    # Extract orbital elements (a is in meters)
    a, e, i, Ω, ω, ν = state.orbital_elements

    # Mean motion [rad/s]
    n = sqrt(μ_earth / a^3)

    # Convert true anomaly to eccentric anomaly using two-argument atan
    E0 = atan(sqrt(1-e^2) * sin(ν), e + cos(ν))

    # Convert eccentric anomaly to mean anomaly
    M0 = E0 - e * sin(E0)

    # Propagate mean anomaly
    M = M0 + n * dt

    # Solve Kepler's equation: E - e*sin(E) = M
    # Using Newton-Raphson iteration
    E = M  # initial guess
    for _ in 1:10  # typically converges in 3-5 iterations
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E))
    end

    # Convert eccentric anomaly back to true anomaly
    # Use two-argument atan for correct quadrant
    ν_new = atan(sqrt(1-e^2) * sin(E), cos(E) - e)

    # Ensure ν is in [0, 2π]
    if ν_new < 0
        ν_new += 2π
    end

    # Return updated orbital elements (a in meters)
    return [a, e, i, Ω, ω, ν_new]
end

"""
    orbital_elements_to_eci(a, e, i, Ω, ω, ν)

Convert Keplerian orbital elements to ECI position and velocity vectors.

# Arguments
- `a::Float64`: Semi-major axis [m]
- `e::Float64`: Eccentricity
- `i::Float64`: Inclination [rad]
- `Ω::Float64`: Right ascension of ascending node (RAAN) [rad]
- `ω::Float64`: Argument of periapsis [rad]
- `ν::Float64`: True anomaly [rad]

# Returns
- `(r_eci, v_eci)`: Tuple of position [m] and velocity [m/s] vectors in ECI frame
"""
function orbital_elements_to_eci(a::Float64, e::Float64, i::Float64, Ω::Float64, ω::Float64, ν::Float64)
    # Compute position and velocity in perifocal frame
    p = a * (1 - e^2)  # semi-latus rectum
    r_mag = p / (1 + e * cos(ν))  # orbital radius

    # Position in perifocal frame (P, Q, W where P points to periapsis)
    r_pf = r_mag * [cos(ν), sin(ν), 0.0]

    # Velocity in perifocal frame
    v_pf = sqrt(μ_earth / p) * [-sin(ν), e + cos(ν), 0.0]

    # Rotation matrix from perifocal to ECI
    # R = R3(-Ω) * R1(-i) * R3(-ω)
    sin_Ω = sin(Ω)
    cos_Ω = cos(Ω)
    sin_i = sin(i)
    cos_i = cos(i)
    sin_ω = sin(ω)
    cos_ω = cos(ω)

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

"""
    orbital_dynamics(state, forces)

Compute orbital state derivatives (position and velocity rates).

# Arguments
- `state::State`: Current state
- `forces::Vector{Float64}`: Total external forces [N] (ECI frame)

# Returns
- `r_dot::Vector{Float64}`: Velocity [m/s]
- `v_dot::Vector{Float64}`: Acceleration [m/s^2]
"""
function orbital_dynamics(state::State, forces::Vector{Float64})
    # TODO: Implement orbital dynamics
    # r_dot = v
    # v_dot = -μ/r^3 * r + forces/m
end
