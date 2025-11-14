using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: calc_potential_energy, calc_kinetic_energy

@testset "Potential Energy Tests" begin
    
    @testset "Values of 0" begin

        r = [7000.0 ,0.0 ,0.0]
        mass = 100.0

        # Radius of 0
        result = calc_potential_energy([0.0, 0.0, 0.0], mass, μ_earth)
        @test result == 0

        # Mass of 0
        result = calc_potential_energy(r, 0.0, μ_earth)
        @test result == 0

        # μ of 0
        result = calc_potential_energy(r, mass, 0.0)
        @test result == 0
    end

    @testset "Reasonable Values" begin
        r = [7000.0 ,0.0 ,0.0] * 1e6
        mass = 100.0

        result = calc_potential_energy(r, mass, μ_earth)
        @test isapprox(result, -5694285.714, rtol=1e-4)

        r = [8e9, 4e10, 3e8] * 1e6
        result = calc_potential_energy(r, mass, μ_earth)
        @test isapprox(result, -.977122219, rtol=1e-4)
        
        mass = 1e12
        result = calc_potential_energy(r, mass, μ_earth)
        @test isapprox(result, -9771222190., rtol=1e-4)

        r = [0.1, 0.1, 0.1]
        result = calc_potential_energy(r, mass, μ_earth)
        @test isapprox(result, -2.301318172e27, rtol=1e-4)
    end
end

@testset "Kinetic Energy Tests" begin
    
    @testset "Values of 0" begin

        v = [7.0, 1.0, 0.0]
        ω = [0.1, 0.1, 0.1]
        mass = 100.0
        I = [
            1.0 2.0 0.0;
            0.0 2.0 0.0;
            0.0 0.9 1.0;
        ]

        # 0 v and ω
        result = calc_kinetic_energy([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], mass, I)
        @test result == 0.0

        # 0 ω only
        result = calc_kinetic_energy(v, [0.0, 0.0, 0.0], mass, I)
        @test result == 2500.0

        # 0 v only
        result = calc_kinetic_energy([0.0, 0.0, 0.0], ω, mass, I)
        @test isapprox(result, 3.45e-2, atol=1e-14)

        # 0 mass only
        result = calc_kinetic_energy(v, ω, 0.0, I)
        @test isapprox(result, 3.45e-2, atol=1e-14)

        # 0 MOI 
        result = calc_kinetic_energy(v, ω, mass, zeros(Float64, 3, 3))
        @test result == 2500.0
    end
    @testset "Reasonable Values" begin

        v = [7.0, 1.0, 0.0]
        ω = [1.0, 1.0, 1.0]
        mass = 100.0
        I = [
            1.0 2.0 0.0;
            0.0 2.0 0.0;
            0.0 0.9 1.0;
        ]

        result = calc_kinetic_energy(v, ω, mass, I)
        @test isapprox(result, 2503.45, atol=1e-2)

        v = [8e9, 4e10, 3e8]
        ω = [1e10, 1e10, 1e10]
        result = calc_kinetic_energy(v, ω, mass, I)
        @test isapprox(result, 8.35495e22, rtol=1e-4)
    end
end

nothing