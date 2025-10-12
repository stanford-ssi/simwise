module Dynamics

using LinearAlgebra

# Call every non-exported name with Simwise.Dynamics.[function]

include("torques/drag.jl")
include("torques/gravityGradient.jl")

include("attitude.jl")
export attitude_dynamics

include("orbit.jl")
export propagate_keplerian, state_from_oe, orbital_elements_to_eci



end