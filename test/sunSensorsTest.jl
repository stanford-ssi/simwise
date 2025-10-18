using Test
using LinearAlgebra
using Random

include("../src/satellite/sensors/sunSensors.jl")

@testset "Sun Sensor Tests" begin
    
    @testset "Basic Functionality" begin
        # Test with sun vector along +X axis
        sun_vec = [1.0, 0.0, 0.0]
        intensities = sensor_relative_intensities(sun_vec; noise_std_intensity=0.0)
        
        @test length(intensities) == NUM_SUN_SENSORS
        @test all(0.0 .<= intensities .<= 1.0)
        
        # Sensors 0-3 (+X pyramid) should have non-zero intensity
        @test intensities[1] > 0.0  # SQRT_2_INV projection
        @test intensities[2] > 0.0
        @test intensities[3] > 0.0
        @test intensities[4] > 0.0
        
        # Sensors 4-7 (-X pyramid) should have zero intensity (facing away)
        @test intensities[5] == 0.0
        @test intensities[6] == 0.0
        @test intensities[7] == 0.0
        @test intensities[8] == 0.0
    end
    
    @testset "Sun Vector Normalization" begin
        # Test that different magnitude vectors give same result
        sun_vec1 = [1.0, 0.0, 0.0]
        sun_vec2 = [5.0, 0.0, 0.0]
        
        intensities1 = sensor_relative_intensities(sun_vec1; noise_std_intensity=0.0)
        intensities2 = sensor_relative_intensities(sun_vec2; noise_std_intensity=0.0)
        
        @test intensities1 ≈ intensities2 atol=1e-10
    end
    
    @testset "Cosine Projection Accuracy" begin
        # Test specific known projections
        sun_vec = [1.0, 0.0, 0.0]  # Along +X
        intensities = sensor_relative_intensities(sun_vec; noise_std_intensity=0.0)
        
        # Sensor 1: normal [SQRT_2_INV, 0, SQRT_2_INV] dot [1,0,0] = SQRT_2_INV
        @test intensities[1] ≈ 1/sqrt(2) atol=1e-10
        
        # Test with sun along +Z
        sun_vec_z = [0.0, 0.0, 1.0]
        intensities_z = sensor_relative_intensities(sun_vec_z; noise_std_intensity=0.0)
        
        # Z+ sensors (indices 13, 14) should have intensity = 1
        @test intensities_z[13] ≈ 1.0 atol=1e-10
        @test intensities_z[14] ≈ 1.0 atol=1e-10
        
        # Z- sensors (indices 15, 16) should have intensity = 0
        @test intensities_z[15] == 0.0
        @test intensities_z[16] == 0.0
    end
    
    @testset "Backside Illumination Clamping" begin
        # Sun vector pointing in -X direction
        sun_vec = [-1.0, 0.0, 0.0]
        intensities = sensor_relative_intensities(sun_vec; noise_std_intensity=0.0)
        
        # +X pyramid sensors should be zero (backside illumination)
        @test intensities[1] == 0.0
        @test intensities[2] == 0.0
        @test intensities[3] == 0.0
        @test intensities[4] == 0.0
        
        # -X pyramid sensors should be non-zero
        @test intensities[5] > 0.0
        @test intensities[6] > 0.0
        @test intensities[7] > 0.0
        @test intensities[8] > 0.0
    end
    
    @testset "Noise Addition" begin
        sun_vec = [1.0, 0.0, 0.0]
        rng = MersenneTwister(1234)
        
        # Test without noise
        intensities_clean = sensor_relative_intensities(sun_vec; noise_std_intensity=0.0, rng=rng)
        
        # Test with noise
        Random.seed!(rng, 1234)  # Reset RNG
        intensities_noisy = sensor_relative_intensities(sun_vec; noise_std_intensity=0.01, rng=rng)
        
        # Should be different due to noise
        @test !(intensities_clean ≈ intensities_noisy)
        
        # But still in valid range
        @test all(0.0 .<= intensities_noisy .<= 1.0)
    end
    
    @testset "Reproducible Noise" begin
        sun_vec = [1.0, 0.0, 0.0]
        
        # Two calls with same seed should give same result
        rng1 = MersenneTwister(1234)
        rng2 = MersenneTwister(1234)
        
        intensities1 = sensor_relative_intensities(sun_vec; noise_std_intensity=0.05, rng=rng1)
        intensities2 = sensor_relative_intensities(sun_vec; noise_std_intensity=0.05, rng=rng2)
        
        @test intensities1 ≈ intensities2 atol=1e-15
    end
    
    @testset "Edge Cases" begin
        # Zero vector (should not crash due to eps() handling)
        sun_vec_zero = [0.0, 0.0, 0.0]
        intensities = sensor_relative_intensities(sun_vec_zero; noise_std_intensity=0.0)
        @test length(intensities) == NUM_SUN_SENSORS
        
        # Very small vector
        sun_vec_small = [1e-15, 0.0, 0.0]
        intensities_small = sensor_relative_intensities(sun_vec_small; noise_std_intensity=0.0)
        @test all(isfinite.(intensities_small))
    end
    
    @testset "Input Validation" begin
        # Wrong vector length should throw assertion error
        @test_throws AssertionError sensor_relative_intensities([1.0, 0.0]; noise_std_intensity=0.0)
        @test_throws AssertionError sensor_relative_intensities([1.0, 0.0, 0.0, 0.0]; noise_std_intensity=0.0)
    end
    
end