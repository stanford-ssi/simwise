# Gravity gradient torque

using LinearAlgebra


using ..Satellite: Parameters, State

using ..Constants: μ_earth
"""
    gravity_gradient_torque(state, params)

Compute gravity gradient torque on spacecraft.

# Arguments
- `state::State`: Current state
- `params::Parameters`: Spacecraft parameters

# Returns
- `Vector{Float64}`: Gravity gradient torque [N·m] (body frame)
"""
function gravity_gradient_torque(state::State, params::Parameters)
    # TODO: Implement gravity gradient torque
    # τ_gg = (3μ/r^5) * (r × I * r)

    # I * r in body frame
    r_b = quat_apply(state.q, state.r_eci)
    r_p = quat_apply(params.principle_axes, r_b)
    Ir_p =  params.principle_diag *. r_p
    Ir_b = quat_apply(quat_inv(params.principle_axes), Ir_p)
    
    # Body frame torque
    (constants.μ_earth/(norm(r_b)^5)) * cross(state.r_b, Ir_b)
end
