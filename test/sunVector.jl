using LinearAlgebra
using SatelliteDynamics

@testset "Sun Vector Tests" begin
    @testset "Comparison with SatelliteDynamics" begin
        # Test against SatelliteDynamics.sun_position at various dates
        println("\n  Comparing sun vector with SatelliteDynamics.sun_position:")

        # Test dates spanning different seasons and years (name, JD)
        test_dates = [
            ("Vernal Equinox 2025", 2460752.5),  # March 20, 2025
            ("Summer Solstice 2025", 2460845.5), # June 21, 2025
            ("Autumnal Equinox 2025", 2460935.5),# September 23, 2025
            ("Winter Solstice 2025", 2461025.5), # December 22, 2025
            ("Random Date 2024", 2460400.5),     # April 24, 2024
            ("Random Date 2026", 2461500.0),     # July 15, 2026
            ("J2000 Epoch", 2451545.0),          # January 1, 2000, 12:00
        ]

        # Angular tolerance (degrees) compared to SatelliteDynamics (TODO: this is a terrible test since we're using SatelliteDynamics lol)
        angle_tolerance_deg = 0.001
        angle_tolerance_rad = angle_tolerance_deg * π / 180

        for (name, jd) in test_dates
            # Compute using our implementation
            sun_vec_our = sun_vector_eci(jd)

            # Compute using SatelliteDynamics (returns position in meters)
            # Use utc date string to create Epoch
            epc = Epoch(SatelliteDynamics.jd_to_caldate(jd)...)
            sun_pos_sd = sun_position(epc)
            sun_vec_sd = sun_pos_sd ./ norm(sun_pos_sd)  # Normalize to unit vector

            # Compute angle between vectors
            cos_angle = dot(sun_vec_our, sun_vec_sd)
            cos_angle = clamp(cos_angle, -1.0, 1.0)  # Clamp for numerical stability
            angle_diff_rad = acos(cos_angle)
            angle_diff_deg = angle_diff_rad * 180 / π

            println("  $name (JD $jd):")
            println("    Our vector:    [$(round(sun_vec_our[1], digits=6)), $(round(sun_vec_our[2], digits=6)), $(round(sun_vec_our[3], digits=6))]")
            println("    SatDyn vector: [$(round(sun_vec_sd[1], digits=6)), $(round(sun_vec_sd[2], digits=6)), $(round(sun_vec_sd[3], digits=6))]")
            println("    Angular error: $(round(angle_diff_deg, digits=4))°")

            # Test that vectors are normalized
            @test isapprox(norm(sun_vec_our), 1.0, atol=1e-10)

            # Test angular error is within tolerance
            @test angle_diff_rad < angle_tolerance_rad
        end

        println("\n  ✓ All sun vector tests passed (within $(angle_tolerance_deg)°)")
    end

    @testset "Seasonal Position Verification" begin
        # Verify that sun vector points in expected general direction during solstices/equinoxes
        println("\n  Verifying seasonal sun positions:")

        # Vernal Equinox 2025 - Sun should be near equatorial plane
        jd_vernal = 2460752.5
        sun_vernal = sun_vector_eci(jd_vernal)
        z_component_vernal = sun_vernal[3]
        println("  Vernal Equinox: z-component = $(round(z_component_vernal, digits=6))")
        @test abs(z_component_vernal) < 0.1  # Should be small near equator

        # Summer Solstice 2025 - Sun should be north of equator
        jd_summer = 2460845.5
        sun_summer = sun_vector_eci(jd_summer)
        z_component_summer = sun_summer[3]
        println("  Summer Solstice: z-component = $(round(z_component_summer, digits=6))")
        @test z_component_summer > 0.3  # Should be significantly positive

        # Winter Solstice 2025 - Sun should be south of equator
        jd_winter = 2461025.5
        sun_winter = sun_vector_eci(jd_winter)
        z_component_winter = sun_winter[3]
        println("  Winter Solstice: z-component = $(round(z_component_winter, digits=6))")
        @test z_component_winter < -0.3  # Should be significantly negative

        println("  ✓ Seasonal positions verified")
    end

    @testset "Unit Vector Property" begin
        # Test that output is always a unit vector for various dates
        jd_range = range(2451545.0, 2461500.0, length=100)  # 2000-2026

        for jd in jd_range
            sun_vec = sun_vector_eci(jd)
            @test isapprox(norm(sun_vec), 1.0, atol=1e-10)
        end

        println("  ✓ Unit vector property verified for 100 dates")
    end

    @testset "Continuity Test" begin
        # Test that sun vector changes smoothly over time
        jd_start = 2460400.0
        dt = 0.1  # 0.1 day steps (2.4 hours)
        n_steps = 100

        max_angle_change = 0.0

        for i in 1:(n_steps-1)
            jd1 = jd_start + (i-1) * dt
            jd2 = jd_start + i * dt

            sun1 = sun_vector_eci(jd1)
            sun2 = sun_vector_eci(jd2)

            cos_angle = dot(sun1, sun2)
            cos_angle = clamp(cos_angle, -1.0, 1.0)
            angle_change = acos(cos_angle)

            max_angle_change = max(max_angle_change, angle_change)
        end

        max_angle_change_deg = max_angle_change * 180 / π
        println("  Maximum angular change over $(dt) days: $(round(max_angle_change_deg, digits=4))°")

        # Sun moves ~1° per day, so over 0.1 days it should move ~0.1°
        # Allow some margin
        @test max_angle_change_deg < 0.2

        println("  ✓ Continuity verified")
    end
end
