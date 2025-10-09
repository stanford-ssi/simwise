# State vector definition

include("math/quaternions.jl")

"""
    State

Satellite state representation combining orbital and attitude state.
x = [q_eci→body, ω_body, t_utc, q_ess (orbital elements)]

# Fields
- `q::Quat`: Attitude quaternion from ECI to body (scalar-first: q0, q1, q2, q3)
- `ω::Vector{Float64}`: Angular velocity [rad/s] (body frame)
- `t::Float64`: Time [MJD - Modified Julian Date]
- `orbital_elements::Vector{Float64}`: [a, e, i, Ω, ω, ν] - semi-major axis [m], eccentricity, inclination [rad], RAAN [rad], argument of periapsis [rad], true anomaly [rad]
"""
mutable struct State
    q::Quat                             # quaternion from eci to body, scalar-first (q0, q1, q2, q3)
    ω::Vector{Float64}                  # angular velocity in body [rad/s]
    t::Float64                          # time [MJD]
    orbital_elements::Vector{Float64}   # keplerian orbital elements, m and rad [a, e, i, Ω, ω, ν] ]
end

# Convenience constructor: accepts Vector{Float64} and converts to Quat
# NOTE: Does NOT auto-normalize - caller must normalize if needed
function State(q::Vector{Float64}, ω::Vector{Float64}, t::Float64, orbital_elements::Vector{Float64})
    @assert length(q) == 4 "Quaternion must have 4 elements (scalar-first: q0, q1, q2, q3)"
    @assert length(ω) == 3 "Angular velocity must have 3 elements"
    @assert length(orbital_elements) == 6 "Orbital elements must have 6 elements [a, e, i, Ω, ω, ν]"
    q_quat = Quat(q[1], q[2], q[3], q[4])  # scalar-first, does NOT normalize
    return State(q_quat, ω, t, orbital_elements)
end


# Operator overloading for RK4 integration
import Base: +, *

function +(s1::State, s2::State)
    # Add quaternions component-wise (for RK4 derivatives)
    q_sum = s1.q + s2.q  # Uses Quat's + operator
    return State(q_sum, s1.ω + s2.ω, s1.t + s2.t, s1.orbital_elements + s2.orbital_elements)
end

function *(a::Real, s::State)
    # Scale quaternion component-wise (for RK4 derivatives)
    q_scaled = a * s.q  # Uses Quat's * operator
    return State(q_scaled, a * s.ω, a * s.t, a * s.orbital_elements)
end

"""
    normalize_quaternion!(state)

Normalize the quaternion in a state to unit norm (in-place).
This should be called after each integration step to prevent numerical drift.
"""
function normalize_quaternion!(state::State)
    normalize!(state.q)  # Uses Quat's normalize! function
    return state
end
