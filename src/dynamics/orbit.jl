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
- Semi-major axis a is stored in km but converted to m for calculations
"""
function propagate_keplerian(state::State, dt::Float64)
    # Extract orbital elements (a is in km, convert to m)
    a_km, e, i, Ω, ω, ν = state.orbital_elements
    a = a_km * 1e3  # Convert km to m

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

    # Return updated orbital elements (a in km)
    return [a_km, e, i, Ω, ω, ν_new]
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
