using Test
using Plots
using LinearAlgebra

using Simwise.World: low_solar_activity, moderate_solar_activity, high_solar_activity, atmosphere_characteristics

# println(typeof(low_solar_activity))
# println()
p = plot([alt for (alt, _) in high_solar_activity], [entry.temp for (_, entry) in high_solar_activity],
        title="Temperature vs Altitude", xlabel="Altitude (km)", ylabel="Temperature (K)")
savefig(p, "../src/world/atmosphere_graphs/temp.png")
p = plot([alt for (alt, _) in high_solar_activity], [log(entry.density) for (_, entry) in high_solar_activity],
        title="Density vs Altitude", xlabel="Altitude (km)", ylabel="log(density)")
savefig(p, "../src/world/atmosphere_graphs/density.png")
p = plot([alt for (alt, _) in high_solar_activity], [log(entry.pressure) for (_, entry) in high_solar_activity],
        title="Pressure vs Altitude", xlabel="Altitude (km)", ylabel="log(pressure)")
savefig(p, "../src/world/atmosphere_graphs/pressure.png")
p = plot([alt for (alt, _) in high_solar_activity], [entry.mol_wt for (_, entry) in high_solar_activity],
        title="Mole Weight vs Altitude", xlabel="Altitude (km)", ylabel="mol_wt (kg/kmol)")
savefig(p, "../src/world/atmosphere_graphs/mol_wt.png")


@testset "Atmosphere Lookup Tests" begin
    
    @testset "Simple Lookup" begin
        altitude_km = 600.0

        entry = atmosphere_characteristics(altitude_km, "low")
        
        @test isapprox(entry.temp, 699.1631)
        @test isapprox(entry.density, 1.03E-14)
        @test isapprox(entry.pressure, 1.06E-08)
        @test isapprox(entry.mol_wt, 5.5149)

        entry = atmosphere_characteristics(altitude_km, "moderate")
        
        @test isapprox(entry.temp, 1011.5325)
        @test isapprox(entry.density, 1.56E-13)
        @test isapprox(entry.pressure, 1.01E-07)
        @test isapprox(entry.mol_wt, 13.0389)


        entry = atmosphere_characteristics(altitude_km) # moderate default
        
        @test isapprox(entry.temp, 1011.5325)
        @test isapprox(entry.density, 1.56E-13)
        @test isapprox(entry.pressure, 1.01E-07)
        @test isapprox(entry.mol_wt, 13.0389)

        entry = atmosphere_characteristics(altitude_km, "high")
        
        @test isapprox(entry.temp, 1622.0421)
        @test isapprox(entry.density, 6.20E-12)
        @test isapprox(entry.pressure, 5.32E-06)
        @test isapprox(entry.mol_wt, 15.7321)

    end

    @testset "Interp Lookup" begin
        altitude_km = 610.0

        entry = atmosphere_characteristics(altitude_km, "low")
        
        @test isapprox(entry.temp, 699.1631, atol=1e-10)
        @test isapprox(entry.density, 9.15e-15, atol = 1e-16)
        @test isapprox(entry.pressure, 9.98e-09, atol = 1e-10)
        @test isapprox(entry.mol_wt, 5.19065, atol=1e-3)

        entry = atmosphere_characteristics(altitude_km, "moderate")
        
        @test isapprox(entry.temp, 1011.5335, atol=1e-10)
        @test isapprox(entry.density, 1.365e-13, atol=1e-18)
        @test isapprox(entry.pressure, 8.98e-08, atol=1e-10)
        @test isapprox(entry.mol_wt, 12.7523, atol=1e-3)

        entry = atmosphere_characteristics(altitude_km) # moderate default
        
        @test isapprox(entry.temp, 1011.5335, atol=1e-10)
        @test isapprox(entry.density, 1.365e-13, atol=1e-18)
        @test isapprox(entry.pressure, 8.98e-08, atol=1e-10)
        @test isapprox(entry.mol_wt, 12.7523, atol=1e-3)

        entry = atmosphere_characteristics(altitude_km, "high")
        
        @test isapprox(entry.temp, 1622.05045, atol=1e-10)
        @test isapprox(entry.density, 5.65e-12, atol=1e-17)
        @test isapprox(entry.pressure, 4.86e-06, atol=1e-11)
        @test isapprox(entry.mol_wt, 15.6889, atol=1e-3)

    end

    @testset "Lower Clamp Lookup" begin
        altitude_km = -10.0

        entry = atmosphere_characteristics(altitude_km, "low")
        
        @test isapprox(entry.temp, 300.2511)
        @test isapprox(entry.density, 1.17E+00)
        @test isapprox(entry.pressure, 1.01E+05)
        @test isapprox(entry.mol_wt, 28.9502)

        entry = atmosphere_characteristics(altitude_km, "moderate")
        
        @test isapprox(entry.temp, 300.2511)
        @test isapprox(entry.density, 1.17E+00)
        @test isapprox(entry.pressure, 1.01E+05)
        @test isapprox(entry.mol_wt, 28.9502)


        entry = atmosphere_characteristics(altitude_km) # moderate default
        
        @test isapprox(entry.temp, 300.2511)
        @test isapprox(entry.density, 1.17E+00)
        @test isapprox(entry.pressure, 1.01E+05)
        @test isapprox(entry.mol_wt, 28.9502)

        entry = atmosphere_characteristics(altitude_km, "high")

        @test isapprox(entry.temp, 300.2511)
        @test isapprox(entry.density, 1.16E+00)
        @test isapprox(entry.pressure, 9.98E+04)
        @test isapprox(entry.mol_wt, 28.9502)

    end

    @testset "Upper Clamp Lookup" begin
        altitude_km = 100000.0

        entry = atmosphere_characteristics(altitude_km, "low")
        
        @test isapprox(entry.temp, 699.1631)
        @test isapprox(entry.density, 1.18E-15)
        @test isapprox(entry.pressure, 2.81E-09)
        @test isapprox(entry.mol_wt, 2.4470)

        entry = atmosphere_characteristics(altitude_km, "moderate")

        @test isapprox(entry.temp, 1011.5379)
        @test isapprox(entry.density, 5.46E-15)
        @test isapprox(entry.pressure, 9.47E-09)
        @test isapprox(entry.mol_wt, 4.8460)


        entry = atmosphere_characteristics(altitude_km) # moderate default
        
        @test isapprox(entry.temp, 1011.5379)
        @test isapprox(entry.density, 5.46E-15)
        @test isapprox(entry.pressure, 9.47E-09)
        @test isapprox(entry.mol_wt, 4.8460)

        entry = atmosphere_characteristics(altitude_km, "high")

        @test isapprox(entry.temp, 1622.0940)
        @test isapprox(entry.density, 4.03E-13)
        @test isapprox(entry.pressure, 3.97E-07)
        @test isapprox(entry.mol_wt, 13.7015)

    end
end

nothing