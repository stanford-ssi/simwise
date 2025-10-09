# High-level propagator - glue code between dynamics and integration

using SatelliteDynamics: Quaternion

"""
    propagate(state, dt, n_steps, params)

Propagate spacecraft state forward in time using numerical integration.

This function acts as the "glue" between low-level dynamics functions
(attitude_dynamics, orbital_dynamics) and the RK4 integrator. It creates
a wrapper function that combines all dynamics and passes it to rk4_step.

# Arguments
- `state::State`: Initial state
- `dt::Float64`: Time step [s]
- `n_steps::Int`: Number of integration steps
- `params::Parameters`: Spacecraft parameters

# Returns
- `State`: Final state after n_steps

# Architecture
```
propagate
├── coupled_dynamics (wrapper function)
│   ├── compute_torques (from torque models)
│   ├── attitude_dynamics (q_dot, ω_dot)
│   ├── compute_forces (from force models)
│   └── orbital_dynamics (doe)
└── rk4_step (numerical integration)
```

# Example
```julia
# Initial state
q0 = [1.0, 0.0, 0.0, 0.0]  # scalar-first quaternion
ω0 = [0.0, 0.0, 0.1]
t0 = 60000.0  # MJD
oe0 = [7000e3, 0.001, deg2rad(45), 0.0, 0.0, 0.0]
state0 = State(q0, ω0, t0, oe0)  # Automatically converts to Quaternion

# Spacecraft parameters
params = Parameters(10.0, diagm([1.0, 1.0, 1.0]), 0.01, 2.2)

# Propagate for 1 orbit
dt = 1.0  # 1 second
n_steps = 6000
final_state = propagate(state0, dt, n_steps, params)
```
"""
function propagate(state::State, dt::Float64, n_steps::Int, params::Parameters)
    """
    coupled_dynamics(s, p)

    Wrapper function that combines all dynamics into a single state derivative.
    This is the function passed to rk4_step.

    # Arguments
    - `s::State`: Current state
    - `p::Parameters`: Spacecraft parameters

    # Returns
    - `State`: State derivative (dq, dω, dt, doe)
    """
    function coupled_dynamics(s::State, p::Parameters)
        # TODO: Compute torques from environment models
        # torques = gravity_gradient_torque(s, p) + drag_torque(s, p) + ...
        torques = zeros(3)  # Placeholder: no external torques yet

        # Attitude dynamics: returns (q_dot::Vector, ω_dot::Vector)
        q_dot, ω_dot = attitude_dynamics(s, torques, p.inertia)

        # TODO: Compute forces from environment models
        # forces = drag_force(s, p) + solar_radiation_pressure(s, p) + ...
        forces = zeros(3)  # Placeholder

        # Orbital dynamics: returns doe (derivative of orbital elements)
        # doe = orbital_dynamics(s, forces)
        doe = zeros(6)  # Placeholder until orbital_dynamics is implemented

        # Time derivative (convert dt in seconds to MJD)
        dt_mjd = 1.0 / 86400.0  # 1 second in MJD units

        # Return complete derivative state
        return State(q_dot, ω_dot, dt_mjd, doe)
    end

    # Time-stepping loop
    current_state = state
    for i in 1:n_steps
        current_state = rk4_step(coupled_dynamics, current_state, dt, params)
    end

    return current_state
end
