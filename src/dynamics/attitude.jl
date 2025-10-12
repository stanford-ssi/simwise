# Attitude dynamics

using LinearAlgebra

using ..Satellite: State, Parameters
using ..Math: Quat

"""
    attitude_dynamics(state, torques, inertia)

Compute attitude state derivatives (quaternion and angular velocity rates).

# Arguments
- `state::State`: Current state
- `torques::Vector{Float64}`: Total external torques [N·m] (body frame)
- `inertia::Matrix{Float64}`: Inertia tensor [kg·m²] (3x3, body frame)

# Returns
- `q_dot::Quat`: Quaternion derivative [q0_dot, q1_dot, q2_dot, q3_dot] (NOT normalized)
- `ω_dot::Vector{Float64}`: Angular acceleration [rad/s^2] (body frame)

# Equations
- Quaternion kinematics: q_dot = 0.5 * Ω(ω) * q
- Euler's equation: ω_dot = I^-1 * (τ - ω × (I * ω))
"""
function attitude_dynamics(state::State, torques::Vector{Float64}, inertia::Matrix{Float64})
    # Extract angular velocity and quaternion from state
    ω = state.ω  # Angular velocity in body frame [rad/s]
    q = state.q  # Quaternion (scalar-first: q0, q1, q2, q3)

    # Quaternion kinematics: q_dot = 0.5 * Ω(ω) * q
    # For scalar-first quaternion [q0, q1, q2, q3] and ω = [ωx, ωy, ωz]:
    # Ω(ω) = [  0   -ωx  -ωy  -ωz ]
    #        [ ωx    0    ωz  -ωy ]
    #        [ ωy  -ωz    0    ωx ]
    #        [ ωz   ωy  -ωx    0  ]

    ωx, ωy, ωz = ω

    # Return as Quat (does NOT auto-normalize - safe for derivatives!)
    q_dot = Quat(
        0.5 * (-ωx * q.q1 - ωy * q.q2 - ωz * q.q3),  # q0_dot
        0.5 * ( ωx * q.q0 + ωz * q.q2 - ωy * q.q3),  # q1_dot
        0.5 * ( ωy * q.q0 - ωz * q.q1 + ωx * q.q3),  # q2_dot
        0.5 * ( ωz * q.q0 + ωy * q.q1 - ωx * q.q2)   # q3_dot
    )

    # Euler's rotation equation: I * ω_dot = τ - ω × (I * ω)
    # Therefore: ω_dot = I^-1 * (τ - ω × (I * ω))

    I_ω = inertia * ω                    # I * ω
    gyroscopic_torque = cross(ω, I_ω)    # ω × (I * ω)
    net_torque = torques - gyroscopic_torque  # τ - ω × (I * ω)

    ω_dot = inertia \ net_torque         # I^-1 * (τ - ω × (I * ω))

    return q_dot, ω_dot
end