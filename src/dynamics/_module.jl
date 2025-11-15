module Dynamics

using LinearAlgebra

include("../constants.jl")

include("energy.jl")
export calc_kinetic_energy, calc_potential_energy

# Call every non-exported name with Simwise.Dynamics.[function]

include("torques/drag.jl")
include("torques/gravityGradient.jl")

include("attitude.jl")
export attitude_dynamics

include("orbit.jl")
export propagate_keplerian, state_from_oe, orbital_elements_to_eci

include("rigidBody.jl")
export rigid_body_derivative


end