# Spacecraft parameters

"""
    Parameters

Spacecraft physical parameters.

# Fields
- `mass::Float64`: Spacecraft mass [kg]
- `inertia::Matrix{Float64}`: Inertia tensor [kg·m^2] (body frame)
- `principle_diag::Vector{Float64}`: Inertia tensor diagonal elements 
- `principle_axes::Quat`: Transformation from body frame to principle_axes frame - constant
- `area::Float64`: Cross-sectional area [m^2]
- `Cd::Float64`: Drag coefficient [-]
"""
struct Parameters
    mass::Float64
    inertia::Matrix{Float64}  # 3x3 inertia tensor
    principle_diag::Vector{Float64}
    principle_axes::Quat
    area::Float64
    Cd::Float64
end
