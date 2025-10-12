# IGRF 14 magnetic field model

using SatelliteToolboxGeomagneticField
using SatelliteToolboxTransformations
using Dates

using Simwise.Math: eci_to_ecef, ned_to_eci, ecef_to_geocentric, jd_to_gmst
"""
    magnetic_field_eci(r_eci, jd)

Compute the geomagnetic field vector in ECI coordinates.

# Arguments
- `r_eci::Vector{Float64}`: Position in ECI frame [m]
- `jd::Float64`: Julian date

# Returns
- `Vector{Float64}`: Magnetic field in ECI frame [T]
"""
function magnetic_field_eci(r_eci::Vector{Float64}, jd::Float64)
    # Convert ECI to ECEF first
    gmst = jd_to_gmst(jd)
    r_ecef = eci_to_ecef(r_eci, gmst)

    # Convert ECEF position to geocentric coordinates
    r, λ, Ω = ecef_to_geocentric(r_ecef)
    λ = clamp(λ, -pi/2, pi/2)
    Ω = clamp(Ω, -pi, pi)

    # Convert Julian date to Year A.D. for IGRF
    # JD 2451545.0 = January 1, 2000, 12:00 TT
    year_ad = 2000.0 + (jd - 2451545.0) / 365.25

    # Compute magnetic field in NED frame using IGRF
    # igrf returns SVector{3} in nT
    B_ned_nT = igrf(year_ad, r, λ, Ω, Val(:geocentric))

    # Convert to regular Vector and to Tesla
    B_ned = [B_ned_nT[1], B_ned_nT[2], B_ned_nT[3]] ./ 1e9

    # Transform from NED to ECI
    B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

    return B_eci
end