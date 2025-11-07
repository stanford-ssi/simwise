module Math

using LinearAlgebra

# Call every non-exported name with Simwise.Math.[function]

include("integrators.jl")
export rk4_step

include("quaternions.jl")
export Quat, quat_apply, quat_conj, quat_inv, quat_mult, q_to_axis_angle, angle_axis_to_q, hamilton_product

include("time.jl")
export jd_to_gmst

include("transforms.jl")
export eci_to_ecef, ned_to_eci, ecef_to_geocentric, rv_to_orbital_elements


end