module Math

using LinearAlgebra

# Call every non-exported name with Simwise.Math.[function]

include("integrators.jl")
export rk4_step

include("quaternions.jl")
export Quat, quat_apply, quat_conj, quat_inv, quat_mult

include("time.jl")
include("transforms.jl")
export eci_to_ecef


end