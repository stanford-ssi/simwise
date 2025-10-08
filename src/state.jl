# State vector definition

"""
    State

Spacecraft state representation combining orbital and attitude state.

# Fields
- `r::Vector{Float64}`: Position vector [m] (ECI frame)
- `v::Vector{Float64}`: Velocity vector [m/s] (ECI frame)
- `q::Vector{Float64}`: Attitude quaternion (scalar-last convention)
- `ω::Vector{Float64}`: Angular velocity [rad/s] (body frame)
- `t::Float64`: Time [s]
"""
mutable struct State
    r::Vector{Float64}  # position [m]
    v::Vector{Float64}  # velocity [m/s]
    q::Vector{Float64}  # quaternion [qx, qy, qz, qw]
    ω::Vector{Float64}  # angular velocity [rad/s]
    t::Float64          # time [s]
end
