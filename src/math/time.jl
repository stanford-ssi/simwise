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

function mjd_to_jd(mjd::Float64)
    return mjd + 2400000.5
end

function jd_to_mjd(jd::Float64)
    return jd - 2400000.5
end

function mjd_to_gmst(mjd::Float64)
    jd = mjd_to_jd(mjd)
    return jd_to_gmst(jd)
end

function jd_to_caldate(jd::Float64)
    return SatelliteDynamics.jd_to_caldate(jd)
end

function caldate_to_jd(date::String)
    return SatelliteDynamics.caldate_to_jd(date)
end

function mjd_to_caldate(mjd::Float64)
    return SatelliteDynamics.mjd_to_caldate(mjd)
end

function caldate_to_mjd(date::String)
    jd = SatelliteDynamics.caldate_to_jd(date)
    return jd_to_mjd(jd)
end