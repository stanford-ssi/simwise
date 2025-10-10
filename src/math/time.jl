# Time conversion utilities

"""
    jd_to_gmst(jd)

Convert Julian date to Greenwich Mean Sidereal Time.

# Arguments
- `jd::Float64`: Julian date

# Returns
- `Float64`: GMST in radians
"""
function jd_to_gmst(jd::Float64)
    # Days from J2000.0
    d = jd - 2451545.0

    # GMST in seconds (simplified formula)
    gmst_sec = 67310.54841 + d * (876600.0 * 3600.0 + 8640184.812866) +
               0.093104 * (d / 36525.0)^2 - 6.2e-6 * (d / 36525.0)^3

    # Convert to radians (modulo 2π)
    gmst_rad = mod(gmst_sec * (2π / 86400.0), 2π)

    return gmst_rad
end
