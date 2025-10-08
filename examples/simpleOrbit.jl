# Simple circular orbit example

using Simwise

# Define spacecraft parameters
params = Parameters(
    mass = 10.0,           # kg
    inertia = diagm([0.1, 0.1, 0.05]),  # kg·m^2
    area = 0.01,           # m^2
    Cd = 2.2               # drag coefficient
)

# Initial orbital state (circular LEO at 400 km)
r0 = [R_earth + 400e3, 0.0, 0.0]  # m
v0 = [0.0, 7669.0, 0.0]           # m/s (circular velocity)

# Initial attitude (identity quaternion, zero angular velocity)
q0 = [0.0, 0.0, 0.0, 1.0]
ω0 = [0.0, 0.0, 0.0]

state0 = State(r0, v0, q0, ω0, 0.0)

# Simulation parameters
dt = 1.0          # s
t_end = 5400.0    # s (1.5 orbits)

# Run propagation
trajectory = propagate(state0, params, dt, t_end)

# Visualize results
plot_orbit(trajectory)
