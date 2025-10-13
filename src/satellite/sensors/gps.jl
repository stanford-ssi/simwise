# GPS SENSOR OUTPUT DATA
# ECEF COORDINATES
# gps_lat::Float64
# gps_lon::Float64
# gps_time::Float64
# gps_data::String
# gps_alive::Boolean
struct GPS_data
    gps_lat::Float64
    gps_lon::Float64
    gps_alt::Float64
    gps_t::Float64

end

function simulate_gps(state::State)
    r_eci = state.r_eci # Satellite Position inn ECI
    t = state.t # Time in Julian Date

    gmst = julian_to_gmst(t)

    r_ecef = eci_to_ecef(r_eci, gmst)
end

function julian_to_gmst(t::Float64)
    # From https://www.astrogreg.com/snippets/greenwichMeanSiderealTime1982.html
    T = (t - 2451545.0)/36525
    th = 100.46061837 + 36000.770053608*T + (0.000387933*T*T) - T*T*T/38710000 # GMST in degree
    th = th % 360
    if th < 0
        th += 360
    end
    return th*pi/180
end



