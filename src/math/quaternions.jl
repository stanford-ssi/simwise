# Custom quaternion type for attitude representation
# Scalar-first convention: [q0, q1, q2, q3] where q0 is the scalar part

using LinearAlgebra

"""
    Quat

Quaternion representation for 3D rotations (scalar-first convention).

A quaternion q = [q0, q1, q2, q3] represents a rotation, where:
- q0 is the scalar part
- q1, q2, q3 are the vector part
- Satisfies: q0² + q1² + q2² + q3² = 1 (unit quaternion)

# Fields
- `q0::Float64`: Scalar component
- `q1::Float64`: First vector component
- `q2::Float64`: Second vector component
- `q3::Float64`: Third vector component

# Note
This type does NOT auto-normalize on construction, allowing it to be used
for quaternion derivatives in RK4 integration.
"""
mutable struct Quat
    q0::Float64
    q1::Float64
    q2::Float64
    q3::Float64
end

# Convenience constructor from vector (scalar-first)
function Quat(v::Vector{Float64})
    @assert length(v) == 4 "Quaternion must have 4 elements [q0, q1, q2, q3]"
    return Quat(v[1], v[2], v[3], v[4])
end

# Identity quaternion
Quat() = Quat(1.0, 0.0, 0.0, 0.0)

"""
    normalize(q::Quat)

Return a normalized copy of the quaternion (unit norm).
"""
function normalize(q::Quat)
    norm = sqrt(q.q0^2 + q.q1^2 + q.q2^2 + q.q3^2)
    return Quat(q.q0/norm, q.q1/norm, q.q2/norm, q.q3/norm)
end

"""
    normalize!(q::Quat)

Normalize the quaternion in-place to unit norm.
"""
function normalize!(q::Quat)
    norm = sqrt(q.q0^2 + q.q1^2 + q.q2^2 + q.q3^2)
    q.q0 /= norm
    q.q1 /= norm
    q.q2 /= norm
    q.q3 /= norm
    return q
end

"""
    quat_mult(q1::Quat, q2::Quat)

Quaternion multiplication (Hamilton product): q1 ⊗ q2

For scalar-first quaternions [s1, v1] and [s2, v2]:
q1 ⊗ q2 = [s1*s2 - v1·v2, s1*v2 + s2*v1 + v1 × v2]
"""
function quat_mult(q1::Quat, q2::Quat)
    # Extract scalar and vector parts
    s1 = q1.q0
    v1 = [q1.q1, q1.q2, q1.q3]
    s2 = q2.q0
    v2 = [q2.q1, q2.q2, q2.q3]

    # Hamilton product
    s_result = s1 * s2 - dot(v1, v2)
    v_result = s1 * v2 + s2 * v1 + cross(v1, v2)

    return Quat(s_result, v_result[1], v_result[2], v_result[3])
end

# Operator overloading for quaternion multiplication
import Base: *
*(q1::Quat, q2::Quat) = quat_mult(q1, q2)

"""
    quat_conj(q::Quat)

Quaternion conjugate: [q0, q1, q2, q3] → [q0, -q1, -q2, -q3]
"""
function quat_conj(q::Quat)
    return Quat(q.q0, -q.q1, -q.q2, -q.q3)
end

"""
    quat_inv(q::Quat)

Quaternion inverse (for unit quaternions, this is the conjugate).
"""
function quat_inv(q::Quat)
    norm_sq = q.q0^2 + q.q1^2 + q.q2^2 + q.q3^2
    return Quat(q.q0/norm_sq, -q.q1/norm_sq, -q.q2/norm_sq, -q.q3/norm_sq)
end

"""
    to_vector(q::Quat)

Convert quaternion to vector [q0, q1, q2, q3] (scalar-first).
"""
function to_vector(q::Quat)
    return [q.q0, q.q1, q.q2, q.q3]
end

# Operator overloading for RK4 integration
import Base: +

"""
    +(q1::Quat, q2::Quat)

Component-wise quaternion addition (for RK4 derivatives).
Note: This is NOT a standard quaternion operation, but needed for RK4.
"""
function +(q1::Quat, q2::Quat)
    return Quat(q1.q0 + q2.q0, q1.q1 + q2.q1, q1.q2 + q2.q2, q1.q3 + q2.q3)
end

"""
    *(a::Real, q::Quat)

Scalar multiplication of quaternion (for RK4 derivatives).
"""
function *(a::Real, q::Quat)
    return Quat(a * q.q0, a * q.q1, a * q.q2, a * q.q3)
end

# Also support q * a
*(q::Quat, a::Real) = a * q
