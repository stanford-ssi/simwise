module Satellite

using LinearAlgebra

# Call every non-exported name with Simwise.Satellite.[function]

include("parameters.jl")
export Parameters

include("state.jl")
export State, normalize_quaternion!


end