# RK4 propagator
using Printf

"""
    rk4_step(dt, t, x, x_dot_func, parameters)

Single step of 4th order Runge-Kutta integration.

# Arguments
- `dt::Float64`: Timestep length
- `t::Float64`: Current time (may not be needed)
- `x::X`: Variable to be integrated over
- `x_dot_func::Function`: Derivative function for the variable
- `parameters::Params`: Generic type for additional parameters

# Returns
- `x`: Runge-Kuttad variable at t + dt
"""
function rk4_step(dt::Float64, t::Float64, x::X, x_dot_func::Function, parameters::Params) where {X, Params}

    # k1 = 0
    # k2 = 0
    # k3 = 0
    # k4 = 0
    # try
    k1 = x_dot_func(t, x, parameters)
    k2 = x_dot_func(t + 0.5*dt, x + 0.5 * dt * k1, parameters)
    k3 = x_dot_func(t + 0.5*dt, x + 0.5 * dt * k2, parameters)
    k4 = x_dot_func(t + dt, x + dt * k3, parameters)
    # catch
    #     @printf "before: %s\n" x
    #     @printf "K1: %s\n" k1
    #     @printf "K2: %s\n" k2
    #     @printf "K3: %s\n" k3
    #     @printf "K4: %s\n" k4
    #     println()
    # end

    


    return x + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
end
