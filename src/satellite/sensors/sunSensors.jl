# Sun sensors model for simwise

using LinearAlgebra
using Printf
using Random

const NUM_SUN_SENSORS = 16

# Precomputed constant: 1/sqrt(2)
const SQRT_2_INV = 1 / sqrt(2.0)

# 16 sensor normals in body frame (each row is a normal vector)
const SUN_SENSOR_NORMALS = [
    # Pyramid group 1 on +X face (0-3)
    [ SQRT_2_INV,  0.0,         SQRT_2_INV],  # sun_pyramid_1_1
    [ SQRT_2_INV,  SQRT_2_INV,  0.0        ],  # sun_pyramid_1_2
    [ SQRT_2_INV,  0.0,        -SQRT_2_INV],  # sun_pyramid_1_3
    [ SQRT_2_INV, -SQRT_2_INV,  0.0        ],  # sun_pyramid_1_4
    # Pyramid group 2 on -X face (4-7)
    [-SQRT_2_INV,  0.0,         SQRT_2_INV],  # sun_pyramid_2_1
    [-SQRT_2_INV, -SQRT_2_INV,  0.0        ],  # sun_pyramid_2_2
    [-SQRT_2_INV,  0.0,        -SQRT_2_INV],  # sun_pyramid_2_3
    [-SQRT_2_INV,  SQRT_2_INV,  0.0        ],  # sun_pyramid_2_4
    # Y+ sensors (8-9) (note: normals point in -y per your listing)
    [ 0.0, -1.0,  0.0 ],  # y+ sensor 1
    [ 0.0, -1.0,  0.0 ],  # y+ sensor 2
    # Y- sensors (10-11)
    [ 0.0,  1.0,  0.0 ],  # y- sensor 1
    [ 0.0,  1.0,  0.0 ],  # y- sensor 2
    # Z+ face sensors (12-13)
    [ 0.0,  0.0,  1.0 ],  # z+ sensor 1
    [ 0.0,  0.0,  1.0 ],  # z+ sensor 2
    # Z- face sensors (14-15)
    [ 0.0,  0.0, -1.0 ],  # z- sensor 1
    [ 0.0,  0.0, -1.0 ]   # z- sensor 2
]

# Global current (A) per sensor (given)
const SUN_SENSOR_CURRENT_A = 0.0016153846

# Global resistance (Ohms) - make mutable so user can change
global RESISTANCE_OHMS = 3700.0

"""
    sensor_relative_intensities(sun_vec::AbstractVector{<:Real})

Compute the relative intensity for each of the 16 sun sensors given a sun vector.
Input `sun_vec` is expected in the body frame (3-vector). The function returns an
array of length 16 with values in [0, 1] representing the cosine projection (clamped at 0).
"""
function sensor_relative_intensities(sun_vec::AbstractVector{<:Real}; noise_std::Real=0.0, rng::AbstractRNG=Random.GLOBAL_RNG)
    @assert length(sun_vec) == 3 "sun_vec must be a 3-element vector"
    s = Vector{Float64}(sun_vec) ./ (norm(sun_vec) == 0 ? eps() : norm(sun_vec))

    intensities = zeros(Float64, NUM_SUN_SENSORS)
    for i in 1:NUM_SUN_SENSORS
        n = SUN_SENSOR_NORMALS[i]
        # Ensure normal is unit (they already are, but normalize for safety)
        n_unit = n ./ norm(n)
        # Cosine projection; clamp negative values to 0 (no illumination from backside)
        intensities[i] = max(0.0, dot(n_unit, s))
        # Add Gaussian noise to intensity measurement if requested
        if noise_std > 0
            intensities[i] += randn(rng) * noise_std
            # ensure still in [0, 1]
            intensities[i] = clamp(intensities[i], 0.0, 1.0)
        end
    end
    return intensities
end

"""
    sensor_voltages_from_intensity(intensities::AbstractVector{<:Real}; I=SUN_SENSOR_CURRENT_A)

Given an array of relative intensities (from `sensor_relative_intensities`), compute
the voltage across each sensor using V = I * R, where I is the provided current scaled
by intensity (i.e., I_sensor = I * intensity) and R is the global `RESISTANCE_OHMS`.
Returns an array of voltages (length 16).
"""
function sensor_voltages_from_intensity(intensities::AbstractVector{<:Real}; I::Real=SUN_SENSOR_CURRENT_A, voltage_noise_std::Real=0.0, rng::AbstractRNG=Random.GLOBAL_RNG)
    @assert length(intensities) == NUM_SUN_SENSORS "intensities must have length $NUM_SUN_SENSORS"
    voltages = zeros(Float64, NUM_SUN_SENSORS)
    for i in 1:NUM_SUN_SENSORS
        I_sensor = I * max(0.0, intensities[i])
        voltages[i] = I_sensor * RESISTANCE_OHMS
        # Add additive Gaussian noise in volts if requested
        if voltage_noise_std > 0
            voltages[i] += randn(rng) * voltage_noise_std
        end
    end
    return voltages
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

# Demo runner: when the file is executed directly, print a sample set of readings
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running sunSensors demo...")
    # Example sun vector in body frame: along +X
    sun_vec = [1.0, 0.0, 0.0]
    println("Using RESISTANCE_OHMS = ", RESISTANCE_OHMS, " Ohms")
    # Show clean readings
    intensities_clean, voltages_clean = sensor_readings_from_sunvec(sun_vec)
    println("Sensor | Intensity (rel) | Voltage (V) (clean)")
    for i in 1:NUM_SUN_SENSORS
        @printf("%6d | %14.6f | %9.6f\n", i, intensities_clean[i], voltages_clean[i])
    end

    # Show noisy readings (example noise): intensity noise std = 0.02, voltage noise std = 0.001 V
    intensity_noise_std = 0.02
    voltage_noise_std = 0.001
    println("\nNow with noise: intensity_noise_std=$(intensity_noise_std), voltage_noise_std=$(voltage_noise_std)")
    rng = MersenneTwister(1234)
    intensities_noisy, voltages_noisy = sensor_readings_from_sunvec(sun_vec; intensity_noise_std=intensity_noise_std, voltage_noise_std=voltage_noise_std, rng=rng)
    println("Sensor | Intensity (rel) | Voltage (V) (noisy)")
    for i in 1:NUM_SUN_SENSORS
        @printf("%6d | %14.6f | %9.6f\n", i, intensities_noisy[i], voltages_noisy[i])
    end
    println("Demo complete")
end
