# Orbital dynamics

using LinearAlgebra

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
