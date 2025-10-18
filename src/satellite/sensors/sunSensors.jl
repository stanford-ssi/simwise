# Sun sensors model for simwise
using LinearAlgebra
using Printf
using Random

include("../parameters.jl")
include("../../constants.jl")

"""
    sensor_relative_intensities(sun_vec::Vector{Float64})

Compute the relative intensity for each of the 16 sun sensors given a sun vector.
Input `sun_vec` is expected in the ECI frame (3-vector) and converts to body frame. The function returns an
array of length 16 with values in [0, 1] representing the cosine projection (clamped at 0).
"""
function sensor_relative_intensities(sun_vec::Vector{Float64}; rng::AbstractRNG=Random.GLOBAL_RNG, noise_std_intensity::Real=0.0)
    @assert length(sun_vec) == 3 "sun_vec must be a 3-element vector"
    s = Vector{Float64}(sun_vec) ./ (norm(sun_vec) == 0 ? eps() : norm(sun_vec))

    intensities = zeros(Float64, NUM_SUN_SENSORS)
    for i in 1:NUM_SUN_SENSORS
        n = SUN_SENSOR_NORMALS[i]
        # Ensure normal is unit (they already are, but normalize for safety)
        n_unit = n ./ norm(n)
        # Cosine projection; clamp negative values to 0 (no illumination from backside)
        intensities[i] = max(0.0, dot(n_unit, s))
        # Add Gaussian noise to intensity measurement
        intensities[i] += randn(rng) * noise_std_intensity
        # ensure still in [0, 1] 
        intensities[i] = clamp(intensities[i], 0.0, 1.0)
    end
    return intensities
end

"""
    sensor_readings_from_sunvec(sun_vec::AbstractVector{<:Real}; I=SUN_SENSOR_CURRENT_A)

Convenience function that takes a sun vector (body frame), computes relative
intensities and returns a tuple (intensities, voltages).
"""
function sensor_readings_from_sunvec(sun_vec::AbstractVector{<:Real}; I::Real=SUN_SENSOR_CURRENT_A, intensity_noise_std::Real=0.0, voltage_noise_std::Real=0.0, rng::AbstractRNG=Random.GLOBAL_RNG)
    intensities = sensor_relative_intensities(sun_vec; noise_std=intensity_noise_std, rng=rng)
    voltages = sensor_voltages_from_intensity(intensities; I=I, voltage_noise_std=voltage_noise_std, rng=rng)
    return intensities, voltages
end

export NUM_SUN_SENSORS, SUN_SENSOR_NORMALS, RESISTANCE_OHMS,
       sensor_relative_intensities, sensor_voltages_from_intensity, sensor_readings_from_sunvec