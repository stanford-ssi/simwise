using Test
using LinearAlgebra
using SatelliteToolboxTransformations

using Simwise.Math: ecef_to_geocentric, ned_to_eci, eci_to_ecef, eci_to_body, body_to_eci

@testset "Transform Tests" begin

    @testset "ECEF to Geocentric" begin

        @testset "On equator at prime meridian" begin
            r_ecef = [6371e3, 0.0, 0.0]  # Earth radius at equator
            r, λ, Ω = ecef_to_geocentric(r_ecef)

            @test isapprox(r, 6371e3)
            @test isapprox(λ, 0.0, atol=1e-10)
            @test isapprox(Ω, 0.0, atol=1e-10)
        end

        @testset "On equator at 90° longitude" begin
            r_ecef = [0.0, 6371e3, 0.0]
            r, λ, Ω = ecef_to_geocentric(r_ecef)

            @test isapprox(r, 6371e3)
            @test isapprox(λ, 0.0, atol=1e-10)
            @test isapprox(Ω, π/2, atol=1e-10)
        end

        @testset "At north pole" begin
            r_ecef = [0.0, 0.0, 6371e3]
            r, λ, Ω = ecef_to_geocentric(r_ecef)

            @test isapprox(r, 6371e3)
            @test isapprox(λ, π/2, atol=1e-10)
            # Longitude undefined at pole, but atan should give 0
            @test isapprox(Ω, 0.0, atol=1e-10)
        end

        @testset "At south pole" begin
            r_ecef = [0.0, 0.0, -6371e3]
            r, λ, Ω = ecef_to_geocentric(r_ecef)

            @test isapprox(r, 6371e3)
            @test isapprox(λ, -π/2, atol=1e-10)
            @test isapprox(Ω, 0.0, atol=1e-10)
        end

        @testset "45° latitude, 45° longitude" begin
            # Position at 45° lat, 45° lon
            r_mag = 6371e3
            λ_expected = π/4
            Ω_expected = π/4

            r_ecef = r_mag * [
                cos(λ_expected) * cos(Ω_expected),
                cos(λ_expected) * sin(Ω_expected),
                sin(λ_expected)
            ]

            r, λ, Ω = ecef_to_geocentric(r_ecef)

            @test isapprox(r, r_mag)
            @test isapprox(λ, λ_expected)
            @test isapprox(Ω, Ω_expected)
        end

    end

    @testset "ECI to ECEF" begin

        @testset "Zero GMST (aligned frames)" begin
            r_eci = [7000e3, 0.0, 0.0]
            gmst = 0.0
            r_ecef = eci_to_ecef(r_eci, gmst)

            @test isapprox(r_ecef, r_eci)
        end

        @testset "90° GMST rotation" begin
            r_eci = [7000e3, 0.0, 0.0]
            gmst = π/2
            r_ecef = eci_to_ecef(r_eci, gmst)

            @test isapprox(r_ecef, [0.0, -7000e3, 0.0], atol=1e-6)
        end

        @testset "180° GMST rotation" begin
            r_eci = [7000e3, 0.0, 0.0]
            gmst = Float64(π)
            r_ecef = eci_to_ecef(r_eci, gmst)

            @test isapprox(r_ecef, [-7000e3, 0.0, 0.0], atol=1e-6)
        end

        @testset "Z-axis unchanged" begin
            r_eci = [0.0, 0.0, 7000e3]
            gmst = π/3  # arbitrary angle
            r_ecef = eci_to_ecef(r_eci, gmst)

            @test isapprox(r_ecef, r_eci)
        end

        @testset "Magnitude preservation" begin
            r_eci = [5000e3, 3000e3, 4000e3]
            gmst = 1.234  # arbitrary angle
            r_ecef = eci_to_ecef(r_eci, gmst)

            @test isapprox(norm(r_ecef), norm(r_eci))
        end

    end

    @testset "NED to ECI" begin

        @testset "North vector at equator, prime meridian" begin
            B_ned = [1.0, 0.0, 0.0]  # North
            λ = 0.0  # Equator
            Ω = 0.0  # Prime meridian
            gmst = 0.0  # Aligned frames

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # North at equator/prime meridian points to +Z in ECEF/ECI
            @test isapprox(B_eci, [0.0, 0.0, 1.0], atol=1e-10)
        end

        @testset "East vector at equator, prime meridian" begin
            B_ned = [0.0, 1.0, 0.0]  # East
            λ = 0.0
            Ω = 0.0
            gmst = 0.0

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # East at prime meridian points to +Y in ECEF/ECI
            @test isapprox(B_eci, [0.0, 1.0, 0.0], atol=1e-10)
        end

        @testset "Down vector at equator, prime meridian" begin
            B_ned = [0.0, 0.0, 1.0]  # Down
            λ = 0.0
            Ω = 0.0
            gmst = 0.0

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # Down at equator/prime meridian points to -X in ECEF/ECI
            @test isapprox(B_eci, [-1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "North vector at north pole" begin
            B_ned = [1.0, 0.0, 0.0]  # North
            λ = π/2  # North pole
            Ω = 0.0
            gmst = 0.0

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # At north pole, NED frame is flipped
            # North points to -X in ECEF
            @test isapprox(B_eci, [-1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "Down vector at north pole" begin
            B_ned = [0.0, 0.0, 1.0]  # Down
            λ = π/2  # North pole
            Ω = 0.0
            gmst = 0.0

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # Down at north pole points to -Z in ECEF
            @test isapprox(B_eci, [0.0, 0.0, -1.0], atol=1e-10)
        end

        @testset "GMST rotation effect" begin
            B_ned = [0.0, 1.0, 0.0]  # East
            λ = 0.0
            Ω = 0.0
            gmst = π/2  # 90° rotation

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            # East at prime meridian with 90° GMST should point to -X in ECI
            @test isapprox(B_eci, [-1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "Magnitude preservation" begin
            B_ned = [1.0, 2.0, 3.0]
            λ = π/6
            Ω = π/4
            gmst = 0.5

            B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

            @test isapprox(norm(B_eci), norm(B_ned))
        end

    end

    @testset "ECI to Body" begin

        @testset "Identity quaternion" begin
            v_eci = [1.0, 2.0, 3.0]
            q = [1.0, 0.0, 0.0, 0.0]  # Identity

            v_body = eci_to_body(v_eci, q)

            @test isapprox(v_body, v_eci)
        end

        @testset "90° rotation about X-axis" begin
            v_eci = [0.0, 0.0, 1.0]  # +Z in ECI
            q = [sqrt(2)/2, sqrt(2)/2, 0.0, 0.0]  # 90° about X

            v_body = eci_to_body(v_eci, q)

            # Active rotation: +Z rotates to -Y when rotating 90° about X
            @test isapprox(v_body, [0.0, -1.0, 0.0], atol=1e-10)
        end

        @testset "90° rotation about Y-axis" begin
            v_eci = [1.0, 0.0, 0.0]  # +X in ECI
            q = [sqrt(2)/2, 0.0, sqrt(2)/2, 0.0]  # 90° about Y

            v_body = eci_to_body(v_eci, q)

            # Active rotation: +X rotates to -Z when rotating 90° about Y
            @test isapprox(v_body, [0.0, 0.0, -1.0], atol=1e-10)
        end

        @testset "90° rotation about Z-axis" begin
            v_eci = [1.0, 0.0, 0.0]  # +X in ECI
            q = [sqrt(2)/2, 0.0, 0.0, sqrt(2)/2]  # 90° about Z

            v_body = eci_to_body(v_eci, q)

            # Active rotation: +X rotates to +Y when rotating 90° about Z
            @test isapprox(v_body, [0.0, 1.0, 0.0], atol=1e-10)
        end

        @testset "180° rotation about Z-axis" begin
            v_eci = [1.0, 0.0, 0.0]
            q = [0.0, 0.0, 0.0, 1.0]  # 180° about Z

            v_body = eci_to_body(v_eci, q)

            @test isapprox(v_body, [-1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "Magnitude preservation" begin
            v_eci = [3.0, 4.0, 5.0]
            q = [0.7071068, 0.0204722, 0.6960552, 0.1228333]  # Arbitrary rotation

            v_body = eci_to_body(v_eci, q)

            @test isapprox(norm(v_body), norm(v_eci))
        end

    end

    @testset "Body to ECI" begin

        @testset "Identity quaternion" begin
            v_body = [1.0, 2.0, 3.0]
            q = [1.0, 0.0, 0.0, 0.0]  # Identity

            v_eci = body_to_eci(v_body, q)

            @test isapprox(v_eci, v_body)
        end

        @testset "90° rotation about X-axis" begin
            v_body = [0.0, -1.0, 0.0]  # -Y in body
            q = [sqrt(2)/2, sqrt(2)/2, 0.0, 0.0]  # 90° about X

            v_eci = body_to_eci(v_body, q)

            # Inverse of eci_to_body: -Y in body becomes +Z in ECI
            @test isapprox(v_eci, [0.0, 0.0, 1.0], atol=1e-10)
        end

        @testset "90° rotation about Y-axis" begin
            v_body = [0.0, 0.0, -1.0]  # -Z in body
            q = [sqrt(2)/2, 0.0, sqrt(2)/2, 0.0]  # 90° about Y

            v_eci = body_to_eci(v_body, q)

            # Inverse of eci_to_body: -Z in body becomes +X in ECI
            @test isapprox(v_eci, [1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "90° rotation about Z-axis" begin
            v_body = [0.0, 1.0, 0.0]  # +Y in body
            q = [sqrt(2)/2, 0.0, 0.0, sqrt(2)/2]  # 90° about Z

            v_eci = body_to_eci(v_body, q)

            # Inverse of eci_to_body: +Y in body becomes +X in ECI
            @test isapprox(v_eci, [1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "180° rotation about Z-axis" begin
            v_body = [-1.0, 0.0, 0.0]
            q = [0.0, 0.0, 0.0, 1.0]  # 180° about Z

            v_eci = body_to_eci(v_body, q)

            @test isapprox(v_eci, [1.0, 0.0, 0.0], atol=1e-10)
        end

        @testset "Magnitude preservation" begin
            v_body = [3.0, 4.0, 5.0]
            q = [0.7071068, 0.0204722, 0.6960552, 0.1228333]  # Arbitrary rotation

            v_eci = body_to_eci(v_body, q)

            @test isapprox(norm(v_eci), norm(v_body))
        end

    end

    @testset "Round-trip Consistency" begin

        @testset "ECI -> ECEF -> Geocentric consistency" begin
            # Start with a position in ECI
            r_eci = [7000e3, 2000e3, 3000e3]
            gmst = 0.5

            # Convert to ECEF
            r_ecef = eci_to_ecef(r_eci, gmst)

            # Convert to geocentric
            r_geo, λ, Ω = ecef_to_geocentric(r_ecef)

            # Verify magnitude
            @test isapprox(r_geo, norm(r_eci))

            # Reconstruct ECEF from geocentric
            r_ecef_reconstructed = r_geo * [
                cos(λ) * cos(Ω),
                cos(λ) * sin(Ω),
                sin(λ)
            ]

            @test isapprox(r_ecef_reconstructed, r_ecef)
        end

        @testset "ECI <-> Body round-trip" begin
            # Test that eci_to_body and body_to_eci are inverses
            v_eci_original = [1.5, 2.5, 3.5]
            q = [0.9799247, 0.0057721, 0.196252, 0.0346327]  # Arbitrary rotation

            # Convert ECI -> Body -> ECI
            v_body = eci_to_body(v_eci_original, q)
            v_eci_reconstructed = body_to_eci(v_body, q)

            @test isapprox(v_eci_reconstructed, v_eci_original, atol=1e-10)
        end

        @testset "Body <-> ECI round-trip" begin
            # Test the reverse direction: Body -> ECI -> Body
            v_body_original = [2.0, 3.0, 4.0]
            q = [0.7071068, 0.0204722, 0.6960552, 0.1228333]  # Arbitrary rotation

            # Convert Body -> ECI -> Body
            v_eci = body_to_eci(v_body_original, q)
            v_body_reconstructed = eci_to_body(v_eci, q)

            @test isapprox(v_body_reconstructed, v_body_original, atol=1e-10)
        end

    end

end

nothing
