# Tumbling CubeSat example - focus on attitude dynamics

using Simwise

# CubeSat parameters (1U)
params = Parameters(
    mass = 1.33,           # kg
    inertia = diagm([0.001, 0.001, 0.0005]),  # kg·m^2
    area = 0.01,           # m^2
    Cd = 2.2
)

# Circular orbit at 400 km
r0 = [R_earth + 400e3, 0.0, 0.0]
v0 = [0.0, 7669.0, 0.0]

# Initial tumble: non-zero angular velocity
q0 = [0.0, 0.0, 0.0, 1.0]
ω0 = [0.1, 0.05, 0.15]  # rad/s - tumbling

state0 = State(r0, v0, q0, ω0, 0.0)

# Simulation parameters
dt = 0.1          # s (smaller timestep for attitude)
t_end = 3600.0    # s (1 orbit)

# Run propagation
trajectory = propagate(state0, params, dt, t_end)

# Visualize attitude evolution
plot_attitude(trajectory)
