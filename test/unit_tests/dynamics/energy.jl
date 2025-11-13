using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: calc_potential_energy

μ_earth_km = μ_earth / 1e9 # since m3/s2

@testset "Kinetic Energy Tests" begin
    
    @testset "Values of 0" begin

        r = [7000,0,0]
        mass = 100

        # Radius of 0
        result = calc_potential_energy([0.0, 0.0, 0.0], mass, μ_earth_km)
        @test result == 0

        # Mass of 0
        result = calc_potential_energy(r, 0.0, μ_earth_km)
        @test result == 0

        # μ of 0
        result = calc_potential_energy(r, mass, 0)
        @test result == 0
    end

    @testset "Reasonable Values" begin
        r = [7000,0,0]
        mass = 100

        result = calc_potential_energy([0.0, 0.0, 0.0], mass, μ_earth_km)
    end
end

nothing