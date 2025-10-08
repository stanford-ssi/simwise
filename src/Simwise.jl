module Simwise

using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox

# Core modules
include("constants.jl")
include("state.jl")
include("parameters.jl")

# Exports
export State, Parameters, rk4_step

# Integration
include("integration/rk4.jl")

# Dynamics
include("dynamics/attitude.jl")
include("dynamics/orbit.jl")
include("dynamics/torques/drag.jl")
include("dynamics/torques/gravityGradient.jl")

# Environment
include("environment/atmosphere.jl")

# Visualization
include("visualization/orbitViz.jl")
include("visualization/attitudeViz.jl")

# Utils
include("utils/transforms.jl")

end # module Simwise
