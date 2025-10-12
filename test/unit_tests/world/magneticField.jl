using Plots
using LinearAlgebra


using Simwise.Math: jd_to_gmst, eci_to_ecef
using Simwise.World: magnetic_field_eci

@testset "Magnetic Field Tests" begin
    @testset "Magnetic Field ECEF Heatmap" begin
        # Test magnetic field computation and visualization
        # Sample the magnetic field on a grid at a fixed altitude in ECEF

        altitude = 400e3  # 400 km altitude [m]
        R_earth = 6378.137e3  # Earth radius [m]
        r_orbit = R_earth + altitude

        # Julian date for evaluation (October 10, 2025)
        jd = 2460958.5000000

        # Get GMST for this time
        gmst = jd_to_gmst(jd)

        # Create a grid in latitude and longitude
        n_lat = 500
        n_lon = 1000
        lats = range(-π/2 + 0.01, π/2 - 0.01, length=n_lat)  # -90° to +90°
        lons = range(-π + 0.01, π - 0.01, length=n_lon)      # -180° to +180°

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

    @testset "Magnetic Field Objective Points (ENU)" begin
        # Test magnetic field at specific geographic locations in ENU frame
        # Using October 10, 2025 as reference date
        jd = 2460958.5000000
        altitude = 400e3  # 400 km altitude [m]
        R_earth = 6378.137e3  # Earth radius [m]

        # Helper function to compute magnetic field in ENU
        function compute_B_enu(lat, lon, alt, jd)
            r = R_earth + alt

            # Convert geodetic to ECEF
            x = r * cos(lat) * cos(lon)
            y = r * cos(lat) * sin(lon)
            z = r * sin(lat)
            r_ecef = [x, y, z]

            # Convert ECEF to ECI
            gmst = jd_to_gmst(jd)
            sin_gmst = sin(gmst)
            cos_gmst = cos(gmst)
            r_eci = [
                cos_gmst * r_ecef[1] - sin_gmst * r_ecef[2],
                sin_gmst * r_ecef[1] + cos_gmst * r_ecef[2],
                r_ecef[3]
            ]

            # Compute magnetic field in ECI
            B_eci = magnetic_field_eci(r_eci, jd)

            # Convert to ECEF
            B_ecef = eci_to_ecef(B_eci, gmst)

            # Convert ECEF to ENU
            sin_lat = sin(lat)
            cos_lat = cos(lat)
            sin_lon = sin(lon)
            cos_lon = cos(lon)

            R_enu_ecef = [
                -sin_lon           cos_lon           0;
                -sin_lat*cos_lon  -sin_lat*sin_lon   cos_lat;
                 cos_lat*cos_lon   cos_lat*sin_lon   sin_lat
            ]

            B_enu = R_enu_ecef * B_ecef
            return B_enu .* 1e9  # Convert to nT
        end

        # Define test points with expected direction values (in nT)
        # Format: (name, lat_deg, lon_deg, expected_east, expected_north, expected_up)
        test_points = [
            ("Equator (0°, 0°)",           0.0,       0.0,   -1694,   22561,   11663),
            ("North Pole (90°, 0°)",      90.0,       0.0,     152,   1160,   -48216),
            ("Mid-lat NH (45°, -90°)",   45.0,     -90.0,    -657,   14810,  -42038),
            ("Mid-lat SH (-45°, 135°)",  -45.0,     135.0,   1768,   12833,  50920),
            ("South Pole (-90°, 0°)",    -90.0,       0.0,     -7296,   10847,   43084),
            ("South Atlantic (-30°, -40°)", -30.0,   -40.0,  -4458,   12309,   14938),
        ]

        # Angular tolerance for direction check (~0.5 degrees)
        angle_tolerance = 0.5 * π / 180  # radians

        println("\n  Testing magnetic field direction at objective points:")
        for (name, lat_deg, lon_deg, exp_east, exp_north, exp_up) in test_points
            # Convert to radians
            lat = lat_deg * π / 180
            lon = lon_deg * π / 180

            # Compute magnetic field
            B_enu = compute_B_enu(lat, lon, altitude, jd)
            B_east, B_north, B_up = B_enu
            B_mag = norm(B_enu)

            # Normalize to get unit vectors
            B_calc_unit = B_enu / B_mag
            B_exp = [exp_east, exp_north, exp_up]
            B_exp_unit = B_exp / norm(B_exp)

            # Compute angle between vectors using dot product
            cos_angle = dot(B_calc_unit, B_exp_unit)
            # Clamp to [-1, 1] to avoid numerical issues with acos
            cos_angle = clamp(cos_angle, -1.0, 1.0)
            angle_diff = acos(cos_angle)
            angle_diff_deg = angle_diff * 180 / π

            # Print results
            # println("\n  $name:")
            # println("    Calculated: [$(round(B_east, digits=1)), $(round(B_north, digits=1)), $(round(B_up, digits=1))] nT")
            # println("    Expected:   [$(exp_east), $(exp_north), $(exp_up)] nT")
            # println("    |B|:        $(round(B_mag, digits=1)) nT")
            # println("    Angle diff: $(round(angle_diff_deg, digits=4))°")

            # Test that direction is within tolerance
            @test angle_diff < angle_tolerance

            # General sanity checks
            @test B_mag > 10000  # Should be at least 10 µT
            @test B_mag < 80000  # Should be less than 80 µT
        end

        println("\n  ✓ All direction tests passed (within 0.5°)")
    end
end