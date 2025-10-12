module Simulation

using LinearAlgebra

# Call every non-exported name with Simwise.Simulation.[function]

include("simulator.jl")
export propagate



end