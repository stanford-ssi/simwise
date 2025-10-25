using Test

@testset "Sensor Tests" begin
    include("magnetometer.jl")
end
@testset "IMU Tests" begin
    include("imu.jl")
end
