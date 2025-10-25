using Test
using LinearAlgebra
using Statistics
using SatelliteToolboxTransformations

using Simwise.Sensors: simulate_magnetometer

@testset "IMU" begin

    # Test that noise is added correctly
    true_magno = [1.0, 0.0, 0.0]  # True magnetic field in ECI frame
    example_quat = [1.0, 0.0, 0.0, 0.0]  # Identity quaternion (no rotation)
    measured_magno = simulate_magnetometer(true_magno, example_quat)

    @test length(measured_magno) == 3
    # Test that noise has expected standard deviation
    N = 1000
    measurements = [simulate_magnetometer(true_magno, example_quat) for _ in 1:N]
    std_first = std([v[1] for v in measurements]) <= 1.8e-8 && std([v[1] for v in measurements]) >= 1.2e-8
    std_second = std([v[2] for v in measurements]) <= 1.8e-8 && std([v[2] for v in measurements]) >= 1.2e-8
    std_third = std([v[3] for v in measurements]) <= 1.8e-8 && std([v[3] for v in measurements]) >= 1.2e-8

    @test std_first
    @test std_second
    @test std_third


    ##Testing Bias
    mean_first = mean([v[1] for v in measurements])
    mean_second = mean([v[2] for v in measurements])
    mean_third = mean([v[3] for v in measurements])
    @test abs(mean_first) + abs(mean_second) + abs(mean_third) != 0.0
    

end

nothing
