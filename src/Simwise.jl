module Simwise

using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox

# Exports
export State, Parameters, Quat, rk4_step, propagate, propagate_keplerian, attitude_dynamics, normalize_quaternion!, μ_earth

# Constants
include("constants.jl")

# Math pure utilities (pure functions, no type dependencies)
include("math/quaternions.jl")
include("math/transforms.jl")

# Core types
include("satellite/state.jl")
include("satellite/parameters.jl")

# Integrators (need State, Parameters)
include("math/integrators.jl")

# Dynamics
include("dynamics/attitude.jl")
include("dynamics/orbit.jl")
include("dynamics/torques/drag.jl")
include("dynamics/torques/gravityGradient.jl")

# World models
include("world/atmosphere.jl")

# Visualization
include("visualization/orbitViz.jl")
include("visualization/attitudeViz.jl")

# Simulator
include("simulator.jl")

end # module Simwise
