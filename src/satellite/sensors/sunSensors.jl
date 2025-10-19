# Sun sensors model for simwise
using LinearAlgebra
using Printf
using Random

include("parameters.jl")
include("constants.jl")

"""
    sensor_relative_intensities(sun_vec::Vector{Float64})

Compute the relative intensity for each of the 16 sun sensors given a sun vector.
Input `sun_vec` is expected in the ECI frame (3-vector) and converts to body frame. The function returns 2 
array's of length 16 for the intensities (clamped at 0) and voltages.
"""
function sensor_relative_intensities(sun_vec_eci::Vector{Float64}, rng::AbstractRNG=Random.GLOBAL_RNG)
    sun_vec_body = eci_to_body(sun_vec_eci, q_eci_to_body)
    s = Vector{Float64}(sun_vec_body) ./ (norm(sun_vec_body) == 0 ? eps() : norm(sun_vec_body))

    intensities = zeros(Float64, NUM_SUN_SENSORS)
    voltages = zeros(Float64, NUM_SUN_SENSORS)

    for i in 1:NUM_SUN_SENSORS
        # 0-7 are pyramids, 8-15 are flat sensors
        n = SUN_SENSOR_NORMALS[i]
        # Ensure normal is unit (they already are, but normalize for safety)
        n_unit = n ./ norm(n)
        # Cosine projection; clamp negative values to 0 (no illumination from backside)
        intensities[i] = max(0.0, dot(n_unit, s))
        # Add Gaussian noise to intensity measurement
        intensities[i] += randn(rng) * noise_std_intensity
        # ensure still in [0, 1] 
        intensities[i] = clamp(intensities[i], 0.0, 1.0) 
        if intensities[i] < 0.0 
            @printf("Warning: Intensity for sensor %d is negative: %f\n", i, intensities[i])
    
        # calculate current and check against reference voltage
        irradiance = intensities[i] * solar_constant  # W/m^2
        current = irradiance * sensor_area * diode_responsitivity  # A
        voltage = current * resistance # V
        voltages[i] = voltage
        if voltage > (i < 8 ? pyramid_ref_voltage : flat_ref_voltage)
            @printf("Warning: Oversaturation! Voltage for sensor %d exceeds reference: %f V\n", i, voltage)
        end
    end
    return intensities, voltages
end




