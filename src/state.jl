# State vector definition

"""
    State

Satellite state representation combining orbital and attitude state.
x = [q_eci→body, ω_body, t_utc, q_ess (orbital elements)]

# Fields
- `q::Vector{Float64}`: Attitude quaternion from ECI to body [qx, qy, qz, qw]
- `ω::Vector{Float64}`: Angular velocity [rad/s] (body frame)
- `t::Float64`: Time [MJD - Modified Julian Date]
- `orbital_elements::Vector{Float64}`: [a, e, i, Ω, ω, ν] - semi-major axis [m], eccentricity, inclination [rad], RAAN [rad], argument of periapsis [rad], true anomaly [rad]
"""
mutable struct State
    q::Vector{Float64}                  # quaternion from eci to body, scalar-first [qw, qx, qy, qz]
    ω::Vector{Float64}                  # angular velocity in eci [rad/s]
    t::Float64                          # time [MJD]
    orbital_elements::Vector{Float64}   # keplerian orbital elements, km and rad [a, e, i, Ω, ω, ν] ]
end

# Operator overloading for RK4 integration
import Base: +, *

function +(s1::State, s2::State)
    return State(s1.q + s2.q, s1.ω + s2.ω, s1.t + s2.t, s1.orbital_elements + s2.orbital_elements)
end

function *(a::Real, s::State)
    return State(a * s.q, a * s.ω, a * s.t, a * s.orbital_elements)
end
