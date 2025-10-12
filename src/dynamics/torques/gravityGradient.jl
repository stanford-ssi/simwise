# Gravity gradient torque

using LinearAlgebra

using ..Satellite: State, Parameters

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
end
