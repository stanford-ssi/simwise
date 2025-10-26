# Simple example showing how easy it is to visualize with Simwise

using Simwise
using Simwise.Visualization
using LinearAlgebra

# Define spacecraft
params = Parameters(
    10.0,                           # mass (kg)
    diagm([0.1, 0.1, 0.05]),       # inertia (kg·m^2)
    0.01,                           # area (m^2)
    2.2                             # drag coefficient
)

# Define initial orbit (circular LEO at 400 km)
a = R_earth + 400e3
e = 0.0
i = 0.0
Ω = 0.0
ω = 0.0
ν = 0.0
oe0 = [a, e, i, Ω, ω, ν]

# Initial attitude (tumbling)
q0 = [1.0, 0.0, 0.0, 0.0]
ω0 = [0.1, 0.05, 0.15]  # rad/s

# Create initial state
using Simwise.Dynamics: state_from_oe
state0 = state_from_oe(q0, ω0, 0.0, oe0)

# Propagate
println("Propagating...")
trajectory = State[]
push!(trajectory, state0)

dt = 1.0
n_steps = 5400

for i in 1:n_steps
    new_state = propagate(trajectory[end], dt, 1, params)
    push!(trajectory, new_state)
end

println("Creating visualization...")

# This is all you need!
visualize(trajectory)

println("Done! Interact with the visualization.")
wait()
