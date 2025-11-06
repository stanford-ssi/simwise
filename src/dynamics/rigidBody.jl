# Rigid Body dynamics

using LinearAlgebra

using ..Satellite: Parameters, State
using ..Constants: μ_earth
using ..Math: Quat, quat_mult


"""
    rigid_body(state, force_eci, torque_body)

Compute attitude state derivatives (quaternion and angular velocity rates).

# Arguments
- `state::State`: Current state
- `dt_s::Float64`: Timestep
- `force_eci::Vector{Float64}`: Total external forces [N·m] (ECI)
- `torque_body::Vector{Float64}`:  Total external torques in [N·m] (body frame)

# Returns
- `next_state::State`: Next state after timestep
- `ω_dot::Vector{Float64}`: Angular acceleration [rad/s^2] (body frame)

# Equations
- Quaternion kinematics: q_dot = 0.5 * Ω(ω) * q
- Euler's equation: ω_dot = I^-1 * (τ - ω × (I * ω))
"""

function rigid_body_derivative(t::Float64, state::Vector{Float64}, parameters::Parameters)
    r = state[1:3]
    v = state[4:6]
    q = state[7:10]
    ω = state[11:13]
    # Position derivative is velocity 
    r_dot = v

    # Velocity derivative is acceleration (Schaub 2.15)
    v_dot = parameters.force_eci / parameters.mass

    # Quaternion derivative is based on hamilton product (Schaub 3.111)
    q_dot = 0.5 * quat_mult(Quat(q), Quat(vcat(0, ω)))

    # Angular derivative based on (Schaub 4.34-35)
    I = parameters.inertia_body
    I_inv = parameters.inertia_body_inv
    τ = parameters.torque_body
    ω_dot = I_inv * (τ - cross(ω, I * ω))

    return vcat(r_dot, v_dot, q_dot, ω_dot)
end