using Plots
using LinearAlgebra

@testset "Magnetic Field Tests" begin
    @testset "Magnetic Field ECEF Heatmap" begin
        # Test magnetic field computation and visualization
        # Sample the magnetic field on a grid at a fixed altitude in ECEF

        altitude = 400e3  # 400 km altitude [m]
        R_earth = 6378.137e3  # Earth radius [m]
        r_orbit = R_earth + altitude

        # Julian date for evaluation (Jan 1, 2024)
        jd = 2460310.5

        # Get GMST for this time
        gmst = jd_to_gmst(jd)

        # Create a grid in latitude and longitude
        n_lat = 1000
        n_lon = 2000
        lats = range(-π/2, π/2, length=n_lat)  # -90° to +90°
        lons = range(-π, π, length=n_lon)      # -180° to +180°

        # Initialize arrays for magnetic field components (ECEF)
        B_x = zeros(n_lat, n_lon)
        B_y = zeros(n_lat, n_lon)
        B_z = zeros(n_lat, n_lon)
        B_mag = zeros(n_lat, n_lon)

        # Initialize arrays for magnetic field components (ENU)
        B_east = zeros(n_lat, n_lon)
        B_north = zeros(n_lat, n_lon)
        B_up = zeros(n_lat, n_lon)

        println("  Computing magnetic field on grid...")
        for (i, lat) in enumerate(lats)
            for (j, lon) in enumerate(lons)
                # Convert spherical ECEF to Cartesian ECEF coordinates
                x = r_orbit * cos(lat) * cos(lon)
                y = r_orbit * cos(lat) * sin(lon)
                z = r_orbit * sin(lat)

                r_ecef = [x, y, z]

                # Convert ECEF to ECI for magnetic_field_eci function
                # (reverse rotation by GMST)
                sin_gmst = sin(gmst)
                cos_gmst = cos(gmst)
                r_eci = [
                    cos_gmst * r_ecef[1] - sin_gmst * r_ecef[2],
                    sin_gmst * r_ecef[1] + cos_gmst * r_ecef[2],
                    r_ecef[3]
                ]

                # Compute magnetic field in ECI
                B_eci = magnetic_field_eci(r_eci, jd)

                # Convert back to ECEF
                B_ecef = eci_to_ecef(B_eci, gmst)

                # Store ECEF components (convert to nT for better visualization)
                B_x[i, j] = B_ecef[1] * 1e9  # T to nT
                B_y[i, j] = B_ecef[2] * 1e9
                B_z[i, j] = B_ecef[3] * 1e9
                B_mag[i, j] = norm(B_ecef) * 1e9

                # Convert ECEF to ENU (East-North-Up) frame
                # ENU transformation matrix at this lat/lon
                sin_lat = sin(lat)
                cos_lat = cos(lat)
                sin_lon = sin(lon)
                cos_lon = cos(lon)

                # Rotation matrix from ECEF to ENU
                R_enu_ecef = [
                    -sin_lon           cos_lon           0;
                    -sin_lat*cos_lon  -sin_lat*sin_lon   cos_lat;
                     cos_lat*cos_lon   cos_lat*sin_lon   sin_lat
                ]

                B_enu = R_enu_ecef * B_ecef

                # Store ENU components (convert to nT)
                B_east[i, j] = B_enu[1] * 1e9
                B_north[i, j] = B_enu[2] * 1e9
                B_up[i, j] = B_enu[3] * 1e9
            end
        end

        println("  Magnetic field range (ECEF):")
        println("    B_x: $(minimum(B_x)) to $(maximum(B_x)) nT")
        println("    B_y: $(minimum(B_y)) to $(maximum(B_y)) nT")
        println("    B_z: $(minimum(B_z)) to $(maximum(B_z)) nT")
        println("    |B|: $(minimum(B_mag)) to $(maximum(B_mag)) nT")

        println("  Magnetic field range (ENU):")
        println("    B_east: $(minimum(B_east)) to $(maximum(B_east)) nT")
        println("    B_north: $(minimum(B_north)) to $(maximum(B_north)) nT")
        println("    B_up: $(minimum(B_up)) to $(maximum(B_up)) nT")

        # Test that magnetic field values are reasonable
        @test all(B_mag .> 0)  # Magnitude should be positive
        @test all(B_mag .< 100000)  # Should be less than 100 µT
        @test minimum(B_mag) > 1000  # Should be more than 1 µT at LEO altitude

        # Create heatmap plots with contour lines
        try
            # Convert to degrees for better readability
            lats_deg = rad2deg.(lats)
            lons_deg = rad2deg.(lons)

            p1 = heatmap(lons_deg, lats_deg, B_x,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_x (ECEF) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_x)), maximum(abs.(B_x))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p1, lons_deg, lats_deg, B_x, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p2 = heatmap(lons_deg, lats_deg, B_y,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_y (ECEF) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_y)), maximum(abs.(B_y))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p2, lons_deg, lats_deg, B_y, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p3 = heatmap(lons_deg, lats_deg, B_z,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_z (ECEF) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_z)), maximum(abs.(B_z))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p3, lons_deg, lats_deg, B_z, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p4 = heatmap(lons_deg, lats_deg, B_mag,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field Magnitude at $(altitude/1e3) km",
                        c=:viridis,
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p4, lons_deg, lats_deg, B_mag, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1400, 1000))
            savefig(p, "plots/magnetic_field_ecef_heatmap.png")
            println("✓ Magnetic field ECEF heatmap saved to test/plots/magnetic_field_ecef_heatmap.png")

            # Create ENU plots
            p5 = heatmap(lons_deg, lats_deg, B_east,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_east (ENU) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_east)), maximum(abs.(B_east))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p5, lons_deg, lats_deg, B_east, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p6 = heatmap(lons_deg, lats_deg, B_north,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_north (ENU) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_north)), maximum(abs.(B_north))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p6, lons_deg, lats_deg, B_north, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p7 = heatmap(lons_deg, lats_deg, B_up,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field B_up (ENU) at $(altitude/1e3) km",
                        c=:RdBu, clims=(-maximum(abs.(B_up)), maximum(abs.(B_up))),
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p7, lons_deg, lats_deg, B_up, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p8 = heatmap(lons_deg, lats_deg, B_mag,
                        xlabel="Longitude [°]", ylabel="Latitude [°]",
                        title="Magnetic Field Magnitude at $(altitude/1e3) km",
                        c=:viridis,
                        colorbar_title="nT", aspect_ratio=1)
            contour!(p8, lons_deg, lats_deg, B_mag, levels=10, color=:black,
                    linewidth=0.5, linealpha=0.4, label=false)

            p_enu = plot(p5, p6, p7, p8, layout=(2, 2), size=(1400, 1000))
            savefig(p_enu, "plots/magnetic_field_enu_heatmap.png")
            println("✓ Magnetic field ENU heatmap saved to test/plots/magnetic_field_enu_heatmap.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end
end
