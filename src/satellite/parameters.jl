# Spacecraft parameters

"""
    Parameters

Spacecraft physical parameters.

# Fields
- `mass::Float64`: Spacecraft mass [kg]
- `inertia::Matrix{Float64}`: Inertia tensor [kg·m^2] (body frame)
- `area::Float64`: Cross-sectional area [m^2]
- `Cd::Float64`: Drag coefficient [-]
- `force_eci::Vector{Float64}`: Total external forces [N·m] (ECI)
- `torque_body::Vector{Float64}`:  Total external torques in [N·m] (body frame)
- `inertia_body_inv::Matrix{Float64}`: # Inertia tensor inverse (faster to precompute)
"""
mutable struct Parameters
    mass::Float64
    inertia_body::Matrix{Float64}  # 3x3 inertia tensor
    area::Float64
    Cd::Float64
    force_eci::Vector{Float64}
    torque_body::Vector{Float64}
    inertia_body_inv::Matrix{Float64}  # 3x3 inertia tensor, faster to precompute
end


function Parameters(mass::Float64, inertia_body::Matrix{Float64}, area::Float64, Cd::Float64)
    @assert length(inertia_body) == 9 "Inertia tensor must have 9 elements for a 3x3 matrix"
    return Parameters(mass, inertia_body, area, Cd, zeros(3), zeros(3), inv(inertia_body))
end

# Defined for just translation
function Parameters(mass::Float64, inertia_body::Matrix{Float64}, force_eci::Vector{Float64}, torque_body::Vector{Float64})
    return Parameters(mass, inertia_body, 0.0, 0.0, force_eci, torque_body, inv(inertia_body))
end
