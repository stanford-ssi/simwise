# RK4 propagator

"""
    rk4_step(state, dt, params)

Single step of 4th order Runge-Kutta integration.

# Arguments
- `state::State`: Current state
- `dt::Float64`: Time step [s]
- `params::Parameters`: Spacecraft parameters

# Returns
- `State`: Updated state at t + dt
"""
function rk4_step(f::Function, state::State, dt::Float64, params::Parameters)
    k1 = f(state, params)
    k2 = f(state + 0.5 * dt * k1, params)
    k3 = f(state + 0.5 * dt * k2, params)
    k4 = f(state + dt * k3, params)

    new_state = state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # Normalize quaternion to prevent numerical drift
    normalize_quaternion!(new_state)

    return new_state
end
