#=
IMU gyroscope simulation that does two things:
Takes ECI vector, converts to body frame using quaternion, 
adds Gaussian noise and bias to simulate gyro measurements
=#
using ..Math: eci_to_body, Quat


const σ = sqrt(1.523209e-6) ##Standard deviation for gyro noise
const bias = randn(3) .* (0.3 * σ) ##Bias for gyro noise in body TODO: find proper value... currently arbitrary

##Function transfroms the angular velocity in ECI frame to body frame then adds noise
function simulate_imu(eci_gyro::Vector{Float64}, quat_eci_to_body::Vector{Float64})
    body_gyro = eci_to_body(eci_gyro, quat_eci_to_body)
    return add_gauss(body_gyro)
end

##Adds gaussian noise to gyro measurements
    
function add_gauss(vec::Vector{Float64})
    noise = randn(length(vec)) .* σ  # mean = 0, stddev = your value
    return (vec + bias) .+ noise

end



#=
true_gyro = [1.0, 1, 1.0]
# Run many simulations
N = 1000
results = [add_gauss(true_gyro) for _ in 1:N]  # array of vectors

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