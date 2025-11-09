using Test
using LinearAlgebra

using Simwise.Math: rv_to_orbital_elements, orbital_elements_to_rv, OrbitalElements, sma_to_orbit_period, orbit_period_to_sma, mean_motion_to_sma
using Simwise.Constants: DEG_TO_RAD, RAD_TO_DEG

μ = 398600.4418 # km3/s2 TODO: change when the constant in Constants is in correct units

@testset "Orbital Element Transform Tests" begin

    @testset "rv to coes" begin
        
        @testset "Vallado Example" begin
            # Vallado 4th Edition, pg 115
            r = [6524.834, 6862.875, 6448.296] # km
            v = [4.901327, 5.533756, -1.976341] # km/s

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 36127.343, rtol=1e-5)
            @test isapprox(orbital_elements.e, 0.832853, atol=1e-5)
            @test isapprox(orbital_elements.i * RAD_TO_DEG, 87.870, rtol=1e-5)
            @test isapprox(orbital_elements.Ω * RAD_TO_DEG, 227.898, rtol=1e-5)
            @test isapprox(orbital_elements.ω * RAD_TO_DEG, 53.38, rtol=1e-4)
            @test isapprox(orbital_elements.ν * RAD_TO_DEG, 92.335, rtol=1e-5)
        end

        @testset "Random Test" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [100.0, 1300.0, -6000.0] # km
            v = [0.1, -1.0, 2.1] # km/s

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 3203.754267903962, rtol=1e-5)
            @test isapprox(orbital_elements.e, 0.9955257908191099, atol=1e-5)
            @test isapprox(orbital_elements.i * RAD_TO_DEG, 93.90569432920545, rtol=1e-5)
            @test isapprox(orbital_elements.Ω * RAD_TO_DEG, 283.91249467651966, rtol=1e-5)
            @test isapprox(orbital_elements.ω * RAD_TO_DEG, 77.26657593515651, rtol=1e-4)
            @test isapprox(orbital_elements.ν * RAD_TO_DEG, 181.10299292301673, rtol=1e-5)
        end

        @testset "Circular" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [6000.0, 4713.675910949817, 1000.0] # km
            v = [4.446, -5.658, 0.0] # km/s (circular orbit)

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 7692.648, rtol=1e-3)
            @test isapprox(orbital_elements.e, 0.0, atol=1e-2)
            @test isapprox(orbital_elements.i * RAD_TO_DEG, 172.53, rtol=1e-4)
            @test isapprox(orbital_elements.Ω * RAD_TO_DEG, 128.1599, rtol=1e-5) 
            # ω, and ν are undefined for circular orbits (can be any value)
        end


        @testset "0 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [6800.0, 0.0, 0.0] # km
            v = [0.0, 8.656278815328648, 0.0] # km/s (equitorial orbit)

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 9421.974578532005, rtol=1e-4)
            @test isapprox(orbital_elements.e, 0.2782829179465395, atol=1e-4)
            @test isapprox(orbital_elements.i, 0.0, rtol=1e-5)
            # Ω, ω, and ν are undefined for equatorial orbits (can be any value)
        end

        @testset "90 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html

            r = [6686, 0.0, 968.35] # km
            v = [0.0, 0.0, -5.899127972862638] # km/s (circular orbit)

            orbital_elements = rv_to_orbital_elements(r, v, μ)

            @test isapprox(orbital_elements.a, 4790.6431, rtol=1e-4)
            @test isapprox(orbital_elements.e, 0.4305248975978197, atol=1e-4)
            @test isapprox(orbital_elements.i * RAD_TO_DEG, 90.0, rtol=1e-5)
            @test isapprox(orbital_elements.Ω * RAD_TO_DEG, 180.0, rtol=1e-5) # Since going down initially, ascends on other side
            @test isapprox(orbital_elements.ω * RAD_TO_DEG, 340.552, rtol=1e-5)
            @test isapprox(orbital_elements.ν * RAD_TO_DEG, 191.206, rtol=1e-5)
        end
        
        @testset "180 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html
            
            r = [7000, 2.0, 0.0] # km
            v = [1.2, -2.0, 0.0] # km/s (equitorial orbit, but reverse)
            
            orbital_elements = rv_to_orbital_elements(r, v, μ)
            
            @test isapprox(orbital_elements.a, 3675.56885, rtol=1e-4)
            @test isapprox(orbital_elements.e, 0.930685, atol=1e-4)
            @test isapprox(orbital_elements.i * RAD_TO_DEG, 180.0, rtol=1e-5)
            # Ω, ω, and ν are undefined for equatorial orbits (can be any value)
        end

    end

    @testset "coes to rv" begin
        @testset "Vallado Example" begin
            # Vallado 4th Edition, pg 119

            coes = OrbitalElements(
                36126.64,
                0.83285,
                87.87 * DEG_TO_RAD,
                227.89 * DEG_TO_RAD,
                53.38 * DEG_TO_RAD,
                92.335 * DEG_TO_RAD,
            )

            r,v = orbital_elements_to_rv(coes, μ)

            @test isapprox(r, [6525.344, 6861.535, 6449.125], atol=0.1)
            @test isapprox(v, [4.902276, 5.533124, -1.975709], atol=0.0001)
        end
        
        @testset "Random Example" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html
            coes = OrbitalElements(
                3203.754267903962,
                0.9955257908191099,
                93.90569432920545 * DEG_TO_RAD,
                283.91249467651966 * DEG_TO_RAD,
                77.26657593515651 * DEG_TO_RAD,
                181.10299292301673 * DEG_TO_RAD,
            )

            r,v = orbital_elements_to_rv(coes, μ)

            @test isapprox(r, [100.0, 1300.0, -6000.0], atol=0.00001)
            @test isapprox(v, [0.1, -1.0, 2.1], atol=0.0001)
        end

        @testset "90 Inclination" begin
            # http://orbitsimulator.com/formulas/OrbitalElements.html
            coes = OrbitalElements(
                4848.9197268518665,
                0.3938736352565609,
                90.0 * DEG_TO_RAD,
                180.0 * DEG_TO_RAD,
                353.8847334059262 * DEG_TO_RAD,
                177.87427980631205 * DEG_TO_RAD,
            )

            r,v = orbital_elements_to_rv(coes, μ)

            @test isapprox(r, [6686, 0.0, 968.35], atol=1e-6)
            @test isapprox(v, [1.0, 0.0, -5.899127972862638], atol=1e-4)
        end

        # Not testing Circular Equatorial, Circular Inclined, and Elliptical Equatorial cases
        # we'd realistically only be using the rv_to_coes function since the source of truth
        # for position are our r and v vectors

    end

end

@testset "Orbital Element Time Tests" begin

    @testset "Mean Motion to Semi-Major Axis" begin
        # ISS data
        # (not much of a validation, since there is no semi-major axis in the TLE, but a sanity-check for ballpark)
        @test isapprox(mean_motion_to_sma(15.4978258, μ), 7071.85952, atol=.001)
    end
    @testset "Orbital Period to Semi-Major Axis" begin
        # Vallado 4th edition pg. 31
        @test isapprox(orbit_period_to_sma(86164.090518, μ), 42164.1696, atol=.001)

        period = 92.9 * 60 # ISS orbital period [s]
        sma = mean_motion_to_sma(15.4978258, μ) # semi-major axis [km]
        @test isapprox(orbit_period_to_sma(period, μ), sma, atol=.001)
    end
    
    @testset "Semi-Major Axis to Orbital Period" begin
        # Vallado 4th edition pg. 31
        @test isapprox(sma_to_orbit_period(42164.1696, μ), 86164.090518, atol=.001)

        period = 92.9 * 60 # ISS orbital period [s]
        sma = mean_motion_to_sma(15.4978258, μ) # semi-major axis [km]
        @test isapprox(sma_to_orbit_period(sma, μ), period, atol=.001)
    end

    
end
    

nothing
