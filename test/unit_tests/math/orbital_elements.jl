using Test
using LinearAlgebra

using Simwise.Math: rv_to_orbital_elements, OrbitalElements
using Simwise.Constants: DEG2RAD
@testset "Orbital Element Transform Tests" begin

    @testset "rv to coes" begin
        
        @testset "Book Example" begin
            r = [6524.834, 6862.875, 6448.296] # km
            v = [4.901327, 5.533756, -1.976341] # km/s
            μ = 398600.4418 # km3/s2

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 36127.343, rtol=1e-5)
            @test isapprox(orbital_elements.e, 0.832853, atol=1e-5)
            @test isapprox(orbital_elements.i, 87.870 * DEG_TO_RAD, rtol=1e-5)
            @test isapprox(orbital_elements.Ω, 227.898 * DEG_TO_RAD, rtol=1e-5)
            @test isapprox(orbital_elements.ω, 53.38 * DEG_TO_RAD, rtol=1e-4)
            @test isapprox(orbital_elements.ν, 92.335 * DEG_TO_RAD, rtol=1e-5)
        end

        @testset "Random Test" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [100.0, 1300.0, -6000.0] # km
            v = [0.1, -1.0, 2.1] # km/s
            μ = 398600.4418 # km3/s2

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 3203.754267903962, rtol=1e-5)
            @test isapprox(orbital_elements.e, 0.9955257908191099, atol=1e-5)
            @test isapprox(orbital_elements.i, 93.90569432920545 * DEG_TO_RAD, rtol=1e-5)
            @test isapprox(orbital_elements.Ω, 283.91249467651966 * DEG_TO_RAD, rtol=1e-5)
            @test isapprox(orbital_elements.ω, 77.26657593515651 * DEG_TO_RAD, rtol=1e-4)
            @test isapprox(orbital_elements.ν, 181.10299292301673 * DEG_TO_RAD, rtol=1e-5)
        end


        @testset "0 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [6800.0, 0.0, 0.0] # km
            v = [0.0, 7.656278815328648, 0.0] # km/s (circular orbit)
            μ = 398600.4418 # km3/s2

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 6800.0, rtol=1e-4)
            @test isapprox(orbital_elements.e, 0.0, atol=1e-4)
            @test isapprox(orbital_elements.i, 0.0, rtol=1e-5)
            # Ω, ω, and ν are undefined for circular orbits (can be any value)
        end

        @testset "00 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [6800.0, 0.0, 0.0] # km
            v = [0.0, 0.0, 7.656278815328648] # km/s (circular orbit)
            μ = 398600.4418 # km3/s2

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 6800.0, rtol=1e-4)
            @test isapprox(orbital_elements.e, 0.0, atol=1e-4)
            @test isapprox(orbital_elements.i, 90.0 * DEG_TO_RAD, rtol=1e-5)
            # Ω, ω, and ν are undefined for circular orbits (can be any value)
        end

    end

end

nothing
