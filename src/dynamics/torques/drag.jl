# Aerodynamic drag torque

using LinearAlgebra

using ..Satellite: Parameters, State

"""
    drag_torque(state, params)

Compute aerodynamic drag torque on spacecraft.

# Arguments
- `state::State`: Current state
- `params::Parameters`: Spacecraft parameters

# Returns
- `Vector{Float64}`: Drag torque [N·m] (body frame)
"""
function drag_torque(state::State, params::Parameters)
    # TODO: Implement drag torque calculation
    # τ_drag = r_cp × F_drag
end
