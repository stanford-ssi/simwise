# Spacecraft parameters

"""
    Parameters

Spacecraft physical parameters.

# Fields
- `mass::Float64`: Spacecraft mass [kg]
- `inertia::Matrix{Float64}`: Inertia tensor [kg·m^2] (body frame)
- `area::Float64`: Cross-sectional area [m^2]
- `Cd::Float64`: Drag coefficient [-]
"""
struct Parameters
    mass::Float64
    inertia::Matrix{Float64}  # 3x3 inertia tensor
    area::Float64
    Cd::Float64
end

const NUM_SUN_SENSORS = 16  # Total number of sun sensors
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