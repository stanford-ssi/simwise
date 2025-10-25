
#include("../math/transforms.jl")
#include("../satellite/state.jl")

using ..Satellite: State
using ..Math: eci_to_body, Quat


#Standard Deviation of Magnetometer [T]
const σ = 1.5e-8 #value between 15 and 30 nT depending on the count number, using the largest number to be safe


##Function transforms the velocity in ECI frame to body frame then adds noise
function simulate_magnetometer(eci_magno::Vector{Float64}, quat_eci_to_body::Vector{Float64})
    body_magno = eci_to_body(eci_magno, quat_eci_to_body)
    return add_gauss(body_magno)
end

##Adds gaussian noise to gyro measurements
    
function add_gauss(vec::Vector{Float64})
    noise = randn(length(vec)) * σ  # mean = 0, stddev = your value
    return vec .+ noise

end


#=


# Example usage and plotting of magnetometer noise distribution
# Need to have quat_state to test properly.

using Plots

true_gyro = [1.0, 1.0, 1.0]
# Run many simulations
N = 1000
# results = [simulate_magnetometer(true_gyro, true_gyro) for _ in 1:N]  ##IMPORTANT: to use this, need quat state and working eci__to_body
results = [add_gauss(true_gyro) for _ in 1:N]  # tests just the distrobution of the noise


# Convert to matrix for plotting
result_matrix = hcat(results...)  # 3 × N matrix: rows are x/y/z axes

# Plot histograms for each axis
plot(
    histogram(result_matrix[1, :], bins=50, label="X-axis", title="Gyro Noise Distribution", xlabel="Value", ylabel="Count"),
    histogram(result_matrix[2, :], bins=50, label="Y-axis"),
    histogram(result_matrix[3, :], bins=50, label="Z-axis"),
    layout = (3,1),
    size=(600, 800)
)
=#