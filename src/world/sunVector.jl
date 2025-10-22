# Sun vector model

using SatelliteDynamics

"""
    sun_vector_eci(jd)
Compute the Sun vector in ECI coordinates.
# Arguments
- `jd::Float64`: Julian date
# Returns
- `Vector{Float64}`: Sun vector in ECI frame (unit vector)
"""
function sun_vector_eci(jd::Float64)
    date = jd_to_caldate(jd)
    epoch = SatelliteDynamics.Epoch(date...)
    sun_pos = SatelliteDynamics.sun_position(epoch)
    # Normalize to unit vector and convert to Vector{Float64}
    return Vector{Float64}(sun_pos ./ norm(sun_pos))
end