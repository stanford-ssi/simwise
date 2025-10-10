using Plots
using LinearAlgebra

@testset "Magnetic Field Tests" begin
    @testset "Magnetic Field ECI Heatmap" begin
        # Test magnetic field computation and visualization
        # Sample the magnetic field on a grid at a fixed altitude

        altitude = 400e3  # 400 km altitude [m]
        R_earth = 6378.137e3  # Earth radius [m]
        r_orbit = R_earth + altitude

        # Julian date for evaluation (Jan 1, 2024)
        jd = 2460310.5

        # Create a grid in latitude and longitude
        n_lat = 50
        n_lon = 100
        lats = range(-π/2, π/2, length=n_lat)  # -90° to +90°
        lons = range(-π, π, length=n_lon)      # -180° to +180°

        # Initialize arrays for magnetic field components
        B_x = zeros(n_lat, n_lon)
        B_y = zeros(n_lat, n_lon)
        B_z = zeros(n_lat, n_lon)
        B_mag = zeros(n_lat, n_lon)

        println("  Computing magnetic field on grid...")
        for (i, lat) in enumerate(lats)
            for (j, lon) in enumerate(lons)
                # Convert spherical to ECI coordinates
                # Note: This is a simplified conversion assuming lon is already in ECI frame
                # For a more accurate representation, we'd need to account for GMST
                x = r_orbit * cos(lat) * cos(lon)
                y = r_orbit * cos(lat) * sin(lon)
                z = r_orbit * sin(lat)

                r_eci = [x, y, z]

                # Compute magnetic field
                B_eci = magnetic_field_eci(r_eci, jd)

                # Store components (convert to nT for better visualization)
                B_x[i, j] = B_eci[1] * 1e9  # T to nT
                B_y[i, j] = B_eci[2] * 1e9
                B_z[i, j] = B_eci[3] * 1e9
                B_mag[i, j] = norm(B_eci) * 1e9
            end
        end

        println("  Magnetic field range:")
        println("    B_x: $(minimum(B_x)) to $(maximum(B_x)) nT")
        println("    B_y: $(minimum(B_y)) to $(maximum(B_y)) nT")
        println("    B_z: $(minimum(B_z)) to $(maximum(B_z)) nT")
        println("    |B|: $(minimum(B_mag)) to $(maximum(B_mag)) nT")

        # Test that magnetic field values are reasonable
        @test all(B_mag .> 0)  # Magnitude should be positive
        @test all(B_mag .< 100000)  # Should be less than 100 µT
        @test minimum(B_mag) > 1000  # Should be more than 1 µT at LEO altitude

        # Create heatmap plots
        try
            # Convert to degrees for better readability
            lats_deg = rad2deg.(lats)
            lons_deg = rad2deg.(lons)

            p1 = heatmap(lons_deg, lats_deg, B_x,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_x (ECI) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_x)), maximum(abs.(B_x))),
                        colorbar_title="nT")

            p2 = heatmap(lons_deg, lats_deg, B_y,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_y (ECI) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_y)), maximum(abs.(B_y))),
                        colorbar_title="nT")

            p3 = heatmap(lons_deg, lats_deg, B_z,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_z (ECI) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_z)), maximum(abs.(B_z))),
                        colorbar_title="nT")

            p4 = heatmap(lons_deg, lats_deg, B_mag,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field Magnitude at $(altitude/1e3) km",
                        c=:viridis,
                        colorbar_title="nT")

            p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1400, 1000))
            savefig(p, "plots/magnetic_field_eci_heatmap.png")
            println("✓ Magnetic field heatmap saved to test/plots/magnetic_field_eci_heatmap.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end
end
