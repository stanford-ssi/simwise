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
function rk4_step(state::State, dt::Float64, params::Parameters)
    # TODO: Implement RK4 integration
end
