using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: calc_potential_energy

μ_earth_km = 398600.0 # km3/s2

@testset "Kinetic Energy Tests" begin
    
    @testset "Values of 0" begin

        r = [7000.0 ,0.0 ,0.0]
        mass = 100.0

        # Radius of 0
        result = calc_potential_energy([0.0, 0.0, 0.0], mass, μ_earth_km)
        @test result == 0

        # Mass of 0
        result = calc_potential_energy(r, 0.0, μ_earth_km)
        @test result == 0

        # μ of 0
        result = calc_potential_energy(r, mass, 0.0)
        @test result == 0
    end

    @testset "Reasonable Values" begin
        r = [7000.0 ,0.0 ,0.0]
        mass = 100.0

        result = calc_potential_energy(r, mass, μ_earth_km)
        @test isapprox(result, -5694.285714, atol=0.00001)

        r = [8e9, 4e10, 3e8]
        result = calc_potential_energy(r, mass, μ_earth_km)
        @test isapprox(result, -9.77122219e-4, atol=1e-10)
        
        mass = 1e12
        result = calc_potential_energy(r, mass, μ_earth_km)
        @test isapprox(result, -9771222.19, atol=1)

        r = [0.1, 0.1, 0.1]
        result = calc_potential_energy(r, mass, μ_earth_km)
        @test isapprox(result, -2.301318172e18, atol=1e13)
    end
end

nothing