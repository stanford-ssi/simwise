module Simwise

using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox
using SatelliteToolboxTransformations
using SatelliteToolboxGeomagneticField

# Exports
export State, Parameters, Quat, rk4_step, propagate, propagate_keplerian, attitude_dynamics, normalize_quaternion!, μ_earth, state_from_oe, orbital_elements_to_eci, magnetic_field_eci, eci_to_geocentric, ned_to_eci, jd_to_gmst

# Constants
include("constants.jl")

# Math pure utilities (pure functions, no type dependencies)
include("math/quaternions.jl")
include("math/transforms.jl")
include("math/time.jl")

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
include("world/magneticField.jl")

# Visualization
include("visualization/orbitViz.jl")
include("visualization/attitudeViz.jl")

# Simulator
include("simulator.jl")

end # module Simwise
