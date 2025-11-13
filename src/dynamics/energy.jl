# Attitude dynamics

using LinearAlgebra

# using ..Satellite: Parameters, State
# using ..Math: Quat

"""
    calc_potential_energy(r, mass_kg, μ)

Compute potential energy (currently only gravity). (Schaub 9.76)

# Arguments
- `r::Vector{Float64}`: Position vector relative to center of central body [km] (typically the earth)
- `mass_kg::Float64`: Mass of satellite / rigid body [kg]
- `μ::Float64`: Gravitational parameter of central body [km3/s2] (typically the earth)

# Returns
- `PE::Float64`: Potential Energy [MJ]
"""
function calc_potential_energy(r::Vector{Float64}, mass_kg::Float64, μ::Float64)
    r_norm = norm(r)
    return -μ*mass_kg/r_norm
end

"""
    calc_kinetic_energy(r, mass_kg, μ)

Compute kinetic energy (currently only gravity). (Schaub 4.55)

# Arguments
- `v::Vector{Float64}`: Velocity vector in ECI [km/s] 
- `ω::Vector{Float64}`: Angular velocity vector in body frame [rad/s]
- `mass_kg::Float64`: Mass of satellite / rigid body [kg]
- `I::Matrix{Float64}`: Moment of inertia tensor in body frame [kg * m * m]

# Returns
- `KE::Float64`: Kinetic Energy [MJ]
"""
function calc_kinetic_energy(v::Vector{Float64}, ω::Vector{Float64}, mass_kg::Float64, I::Matrix{Float64})
    I_km2 = I * 1e6
    v_norm = norm(v)
    KE_trans = 0.5 * mass_kg * v_norm * v_norm 
    KE_rot = 0.5 * (transpose(ω) * I_km2 * ω)
    return KE_trans + KE_rot
end