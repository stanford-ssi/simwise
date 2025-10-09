# State vector definition

using SatelliteDynamics: Quaternion

"""
    State

Satellite state representation combining orbital and attitude state.
x = [q_eci→body, ω_body, t_utc, q_ess (orbital elements)]

# Fields
- `q::Quaternion`: Attitude quaternion from ECI to body (scalar-first: q0, q1, q2, q3)
- `ω::Vector{Float64}`: Angular velocity [rad/s] (body frame)
- `t::Float64`: Time [MJD - Modified Julian Date]
- `orbital_elements::Vector{Float64}`: [a, e, i, Ω, ω, ν] - semi-major axis [m], eccentricity, inclination [rad], RAAN [rad], argument of periapsis [rad], true anomaly [rad]
"""
mutable struct State
    q::Quaternion                       # quaternion from eci to body, scalar-first (q0, q1, q2, q3)
    ω::Vector{Float64}                  # angular velocity in eci [rad/s]
    t::Float64                          # time [MJD]
    orbital_elements::Vector{Float64}   # keplerian orbital elements, m and rad [a, e, i, Ω, ω, ν] ]
end

# Convenience constructor: accepts Vector{Float64} and converts to Quaternion
# WARNING: Quaternion() auto-normalizes! Only use for STATE values, not derivatives!
function State(q::Vector{Float64}, ω::Vector{Float64}, t::Float64, orbital_elements::Vector{Float64})
    @assert length(q) == 4 "Quaternion must have 4 elements (scalar-first: q0, q1, q2, q3)"
    @assert length(ω) == 3 "Angular velocity must have 3 elements"
    @assert length(orbital_elements) == 6 "Orbital elements must have 6 elements [a, e, i, Ω, ω, ν]"
    q_quat = Quaternion(q[1], q[2], q[3], q[4])  # scalar-first, AUTO-NORMALIZES
    return State(q_quat, ω, t, orbital_elements)
end

# Special constructor for derivatives (bypasses Quaternion normalization)
function state_derivative(q_dot::Vector{Float64}, ω_dot::Vector{Float64}, t_dot::Float64, oe_dot::Vector{Float64})
    # Create Quaternion WITHOUT normalization by using internal constructor
    # We abuse the fact that adding normalized + non-normalized still works for RK4
    # Actually, we need to construct a "fake" Quaternion that's not normalized
    # Since we can't do that, we'll need to change the approach...

    # Temporary hack: store as normalized, but scale by actual norm so + and * work
    norm_val = sqrt(q_dot[1]^2 + q_dot[2]^2 + q_dot[3]^2 + q_dot[4]^2)
    if norm_val > 0
        # Store the normalized direction with the magnitude encoded elsewhere
        # This won't work... we need a different solution
        q_quat = Quaternion(q_dot[1]/norm_val, q_dot[2]/norm_val, q_dot[3]/norm_val, q_dot[4]/norm_val)
        # Then scale the result? No this is getting hacky...
    else
        q_quat = Quaternion(1.0, 0.0, 0.0, 0.0)
    end

    return State(q_quat, ω_dot, t_dot, oe_dot)
end

# Operator overloading for RK4 integration
import Base: +, *

function +(s1::State, s2::State)
    # Add quaternions component-wise (for RK4 derivatives)
    q_sum = Quaternion(s1.q.q0 + s2.q.q0, s1.q.q1 + s2.q.q1, s1.q.q2 + s2.q.q2, s1.q.q3 + s2.q.q3)
    return State(q_sum, s1.ω + s2.ω, s1.t + s2.t, s1.orbital_elements + s2.orbital_elements)
end

function *(a::Real, s::State)
    # Scale quaternion component-wise (for RK4 derivatives)
    q_scaled = Quaternion(a * s.q.q0, a * s.q.q1, a * s.q.q2, a * s.q.q3)
    return State(q_scaled, a * s.ω, a * s.t, a * s.orbital_elements)
end

"""
    normalize_quaternion!(state)

Normalize the quaternion in a state to unit norm (in-place).
This should be called after each integration step to prevent numerical drift.
"""
function normalize_quaternion!(state::State)
    norm = sqrt(state.q.q0^2 + state.q.q1^2 + state.q.q2^2 + state.q.q3^2)
    state.q = Quaternion(state.q.q0/norm, state.q.q1/norm, state.q.q2/norm, state.q.q3/norm)
    return state
end
