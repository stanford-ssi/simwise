module Simwise

using LinearAlgebra
using SatelliteDynamics
using SatelliteToolbox
using SatelliteToolboxTransformations
using SatelliteToolboxGeomagneticField

# TODO: figure out what happens to these export things


# Optional in case we don't want it I guess
# __precompile__(false)

# Constants
include("constants.jl")
export RAD_TO_DEG, DEG_TO_RAD, μ_earth

include(joinpath(@__DIR__, "math", "_module.jl"))
export Quat, quat_apply, quat_conj, quat_inv, quat_mult

include(joinpath(@__DIR__, "satellite", "_module.jl"))
export Parameters, State

include(joinpath(@__DIR__, "dynamics", "_module.jl"))

include(joinpath(@__DIR__, "simulation", "_module.jl"))
export propagate

include(joinpath(@__DIR__, "visualization", "_module.jl"))

include(joinpath(@__DIR__, "world", "_module.jl"))
export magnetic_field_eci, sun_vector_eci

# Exports
# export μ_earth

end # module Simwise
