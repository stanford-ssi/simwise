using Test
using LinearAlgebra
using Statistics
using SatelliteToolboxTransformations

using Simwise.Sensors: simulate_magnetometer

# constants to define standard deviation
const σ = 1.5e-5  # in μT
const ε = 0.3e-5  # in μT

@testset "Magnetometer Sensor Tests" begin
    # Test that noise is added correctly
    b_true_eci = [1.0, 0.0, 0.0]  # True magnetic field in ECI frame
    example_quat = [1.0, 0.0, 0.0, 0.0]  # Identity quaternion (no rotation)
    b_measured = simulate_magnetometer(b_true_eci, example_quat)

    @test length(b_measured) == 3
    N = 1000
    measurements = [simulate_magnetometer(b_true_eci, example_quat) for _ in 1:N]

    # Test that noise has expected standard deviation
    std_first = std([v[1] for v in measurements]) <= (σ + ε) && std([v[1] for v in measurements]) >= (σ - ε)
    std_second = std([v[2] for v in measurements]) <= (σ + ε) && std([v[2] for v in measurements]) >= (σ - ε)
    std_third = std([v[3] for v in measurements]) <= (σ + ε) && std([v[3] for v in measurements]) >= (σ - ε)

    @test std_first
    @test std_second
    @test std_third
end

@testset "Plotting measurements" begin
    # Example usage and plotting of magnetometer noise distribution
    # Need to have quat_state to test properly.

    using Plots

    true_gyro = [1.0, 1.0, 1.0]
    # Run many simulations
    N = 1000
    # results = [simulate_magnetometer(true_gyro, true_gyro) for _ in 1:N]  
    # IMPORTANT: to use this, need quat state and working eci_to_body
    results = [add_gauss(true_gyro) for _ in 1:N]  # tests just the distribution of the noise

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
end

nothing
