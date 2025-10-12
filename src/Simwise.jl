module Simwise

using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox
using SatelliteToolboxTransformations
using SatelliteToolboxGeomagneticField

# Optional in case we don't want it I guess
# __precompile__(false)

# Constants
include("constants.jl")

include(joinpath(@__DIR__, "math", "module.jl"))
export Quat, quat_apply, quat_conj, quat_inv, quat_mult

include(joinpath(@__DIR__, "satellite", "module.jl"))
export Parameters, State

include(joinpath(@__DIR__, "dynamics", "module.jl"))

include(joinpath(@__DIR__, "simulation", "module.jl"))
export propagate

include(joinpath(@__DIR__, "visualization", "module.jl"))

include(joinpath(@__DIR__, "world", "module.jl"))
export magnetic_field_eci, sun_vector_eci

# Exports
# export μ_earth

end # module Simwise
