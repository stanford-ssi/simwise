# Spacecraft parameters

"""
    Parameters

Spacecraft physical parameters.

# Fields
- `mass::Float64`: Spacecraft mass [kg]
- `inertia::Matrix{Float64}`: Inertia tensor [kg·m^2] (body frame)
- `area::Float64`: Cross-sectional area [m^2]
- `Cd::Float64`: Drag coefficient [-]
"""
struct Parameters
    mass::Float64
    inertia::Matrix{Float64}  # 3x3 inertia tensor
    area::Float64
    Cd::Float64
end
