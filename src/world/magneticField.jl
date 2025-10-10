# IGRF 14 magnetic field model

using SatelliteToolboxGeomagneticField
using SatelliteToolboxTransformations
using Dates

include("../math/transforms.jl")
include("../math/time.jl")

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
    # Convert ECI position to geocentric coordinates
    r, λ, Ω = eci_to_geocentric(r_eci)

    # Convert Julian date to Year A.D. for IGRF
    # JD 2451545.0 = January 1, 2000, 12:00 TT
    year_ad = 2000.0 + (jd - 2451545.0) / 365.25

    # Compute magnetic field in NED frame using IGRF
    # igrf returns SVector{3} in nT
    B_ned_nT = igrf(year_ad, r, λ, Ω, Val(:geocentric))

    # Convert to regular Vector and to Tesla
    B_ned = [B_ned_nT[1], B_ned_nT[2], B_ned_nT[3]] ./ 1e9

    # Get GMST for ECI/ECEF rotation
    gmst = jd_to_gmst(jd)

    # Transform from NED to ECI
    B_eci = ned_to_eci(B_ned, λ, Ω, gmst)

    return B_eci
end