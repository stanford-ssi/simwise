# Demo: Visualize a tumbling satellite in orbit
# This example shows how to use Simwise's visualization capabilities

using Simwise
using Simwise.Satellite: Parameters, State
using Simwise.Simulation: propagate
using Simwise.Visualization
using Simwise.Dynamics: state_from_oe
using Simwise: R_earth, μ_earth
using LinearAlgebra

println("=== Simwise Satellite Visualization Demo ===\n")

# Define a CubeSat
println("Setting up CubeSat (1U)...")
params = Parameters(
    1.33,                           # mass: 1.33 kg (typical 1U CubeSat)
    diagm([0.001, 0.001, 0.0005]), # inertia tensor (kg·m^2)
    0.01,                           # cross-sectional area (m^2)
    2.2                             # drag coefficient
)

# Define orbit: ISS-like orbit at ~400 km altitude
println("Setting up orbit...")
altitude = 400e3  # 400 km
a = R_earth + altitude              # semi-major axis
e = 0.001                           # nearly circular
i = deg2rad(51.6)                   # ISS inclination
Ω = deg2rad(45)                     # RAAN
ω = deg2rad(0)                      # argument of periapsis
ν = deg2rad(0)                      # true anomaly
oe0 = [a, e, i, Ω, ω, ν]

# Initial attitude: start with identity, but give it a spin!
println("Setting up tumbling motion...")
q0 = [1.0, 0.0, 0.0, 0.0]          # identity quaternion
ω0 = [0.15, 0.08, 0.12]            # tumbling at ~10 deg/s on each axis

# Create initial state
state0 = state_from_oe(q0, ω0, 0.0, oe0)

# Propagate for multiple orbits
dt = 1.0          # 1 second timestep
n_orbits = 2.0    # simulate 2 complete orbits

# Calculate orbital period
orbital_period = 2π * sqrt(a^3 / μ_earth)
n_steps = Int(floor(n_orbits * orbital_period / dt))

println("\nPropagating satellite...")
println("  Altitude: $(altitude/1e3) km")
println("  Orbital period: $(round(orbital_period/60, digits=1)) minutes")
println("  Simulation time: $(round(n_orbits * orbital_period/60, digits=1)) minutes")
println("  Number of steps: $n_steps")
println()

trajectory = State[]
push!(trajectory, state0)

# Propagate with progress indicator
print("Progress: ")
progress_interval = max(1, n_steps ÷ 10)
for i in 1:n_steps
    if i % progress_interval == 0
        print("$(round(Int, 100*i/n_steps))% ")
        flush(stdout)
    end
    new_state = propagate(trajectory[end], dt, 1, params)
    push!(trajectory, new_state)
end
println("✓")

println("\nLaunching interactive visualization...")
println("\nControls:")
println("  • Drag the time slider to scrub through the simulation")
println("  • Click 'Play' to animate")
println("  • Use mouse to rotate the 3D plots:")
println("    - Left-click + drag: Rotate view")
println("    - Scroll: Zoom in/out")
println("    - Right-click + drag: Pan")
println("\nClose the window when done.\n")

# Launch the visualization!
visualize(trajectory)

wait()
