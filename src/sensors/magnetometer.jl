
#include("../math/transforms.jl")
#include("../satellite/state.jl")

using ..Satellite: State
using ..Math: eci_to_body, Quat


# Standard Deviation of Magnetometer [T]; from datasheet https://www.tri-m.com/products/pni/RM3100-User-Manual.pdf
const σ = 1.5e-5   # value between 15 and 30 μT depending on the count number; using the largest number to be safe

# Function transforms the velocity in ECI frame to body frame, then adds noise
function read_magnetometer(b_ref_eci::Vector{Float64}, quat_eci_to_body::Vector{Float64})
    b_ref_body = eci_to_body(b_ref_eci, quat_eci_to_body)
    return add_gauss(b_ref_body)
end

# Adds gaussian noise to gyro measurements
function add_gauss(vec::Vector{Float64})
    noise = randn(length(vec)) * σ  # mean = 0, stddev = your value
    return vec .+ noise

end