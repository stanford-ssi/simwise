# Custom quaternion type for attitude representation
# Scalar-first convention: [q0, q1, q2, q3] where q0 is the scalar part

using LinearAlgebra

using Simwise: RAD_TO_DEG, DEG_TO_RAD

"""
    Quat

Quaternion representation for 3D rotations (scalar-first convention).

A quaternion q = [w, x, y, z] represents a rotation, where:
- w is the scalar part
- x, y, z are the vector part
- Satisfies: w² + x² + y² + z² = 1 (unit quaternion)

# Fields
- `w::Float64`: Scalar component
- `x::Float64`: First vector component
- `y::Float64`: Second vector component
- `z::Float64`: Third vector component

# Note
This type does NOT auto-normalize on construction, allowing it to be used
for quaternion derivatives in RK4 integration.
"""

###########################################################################
#           Quaternion Struct
###########################################################################
mutable struct Quat
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

# Convenience constructor from vector (scalar-first)
function Quat(v::Vector{Float64})
    @assert length(v) == 4 "Quaternion must have 4 elements [w, x, y, z]"
    return Quat(v[1], v[2], v[3], v[4])
end

# Identity quaternion
Quat() = Quat(1.0, 0.0, 0.0, 0.0)

"""
    to_vector(q::Quat)

Convert quaternion to vector [w, x, y, z] (scalar-first).
"""
function to_vector(q::Quat)
    return [q.w, q.x, q.y, q.z]
end

function quat_unit(q::Quat)
    return [q.w, q.x, q.y, q.z]
end

###########################################################################
#           Helpers
###########################################################################


"""
    quat_conj(q::Quat)

Quaternion conjugate: [w, x, y, z] → [w, -x, -y, -z]
"""
function quat_conj(q::Quat)
    return Quat(q.w, -q.x, -q.y, -q.z)
end

"""
    quat_inv(q::Quat)

Quaternion inverse (for unit quaternions, this is the conjugate).
"""
function quat_inv(q::Quat)
    norm_sq = q.w^2 + q.x^2 + q.y^2 + q.z^2
    return Quat(q.w/norm_sq, -q.x/norm_sq, -q.y/norm_sq, -q.z/norm_sq)
end

###########################################################################
#           Normalization
###########################################################################

"""
    normalize(q::Quat)

Return a normalized copy of the quaternion (unit norm).
"""
function normalize(q::Quat)
    norm = sqrt(q.w^2 + q.x^2 + q.y^2 + q.z^2)
    return Quat(q.w/norm, q.x/norm, q.y/norm, q.z/norm)
end

"""
    normalize!(q::Quat)

Normalize the quaternion in-place to unit norm.
"""
function normalize!(q::Quat)
    norm = sqrt(q.w^2 + q.x^2 + q.y^2 + q.z^2)
    q.w /= norm
    q.x /= norm
    q.y /= norm
    q.z /= norm
    return q
end

###########################################################################
#           Quaternion Multiplication/Applying
###########################################################################

"""
    quat_mult(q1::Quat, q2::Quat)

Quaternion multiplication (Hamilton product): q1 ⊗ q2

For scalar-first quaternions [s1, v1] and [s2, v2]:
q1 ⊗ q2 = [s1*s2 - v1·v2, s1*v2 + s2*v1 + v1 × v2]
"""
function quat_mult(q1::Quat, q2::Quat)
    # Extract scalar and vector parts
    s1 = q1.w
    v1 = [q1.x, q1.y, q1.z]
    s2 = q2.w
    v2 = [q2.x, q2.y, q2.z]

    # Hamilton product
    s_result = s1 * s2 - dot(v1, v2)
    v_result = s1 * v2 + s2 * v1 + cross(v1, v2)

    return Quat(s_result, v_result[1], v_result[2], v_result[3])
end

"""Equivalent to expressing the SAME vector in a DIFFERENT coordinate frame

Types of rotations specified by the boolean argument `passive`:
- "passive" (true) - Rotates the coordinate frame, and produces the vector representation in that new frame
- "active" (false) - Rotates a vector V around a FIXED axis to product V'
Defaults to "passive".

For active rotation, returns the new vector. For passive rotation, the vector expressed in the new coordinate system.
"""
function quat_apply(q::Quat, v::Vector{Float64}, passive::Bool = true)

    normalize!(q)

    # If passive rotation, the order of q and q* are reversed
    if passive
        # Passive: v' = q * v * q*
        result = quat_mult(quat_mult(quat_conj(q), Quat(0, v[1], v[2], v[3])), q)
    else
        # Active: v' = q* * v * q
        result = quat_mult(quat_mult(q, Quat(0, v[1], v[2], v[3])), quat_conj(q))
    end

    return Float64[result.x, result.y, result.z]
end

###########################################################################
#           Arithmetic Overloading
###########################################################################

function angle_axis_to_q(angle::Float64, axis::Vector, degrees::Bool = false)
    """Turn angle-axis rotation into quaternion.
    The `degrees` argument specified whether the angle is in degrees.
    Defaults to false
    """

    if degrees
        angle = angle * DEG_TO_RAD
    end

    unit_axis = axis / norm(axis)

    w = cos(angle/2)
    x,y,z = sin(angle / 2) * unit_axis

    return Quat(w,x,y,z)
end

function q_to_axis_angle(q::Quat, degrees::Bool = false)
    """Turn quaternion rotation into angle-axis.
    The `degrees` argument specified whether the angle is in degrees.
    Defaults to false.

    If angle is zero, returns 0 angle and [0,0,0]
    """
    angle = 2*acos(q.w)

    # Shouldn't realistically encounter this
    if angle < 1e-10
        return 0, [0,0,0]
    end

    i = q.x/sin(angle/2)
    j = q.y/sin(angle/2)
    k = q.z/sin(angle/2)

    if degrees
        angle = angle * RAD_TO_DEG
    end

    return angle, [i,j,k]
end

###########################################################################
#           Direction Cosine Matrix
###########################################################################

"""
    https://motoq.github.io/doc/tnotes/dcm
"""
function dcm_to_q(dcm::Matrix{Float64})
    trace = dcm

    c11, c12, c13 = dcm[1,:]
    c21, c22, c23 = dcm[2,:]
    c31, c32, c33 = dcm[3,:]

    tr = c11 + c22 + c33

    if tr > c11 && tr > c22 && tr > c33
        w = sqrt((1 + c11 + c22 + c33) / 4)
        x = (c23 - c32) / 4 / w
        y = (c31 - c13) / 4 / w
        z = (c12 - c21) / 4 / w
    elseif c11 > c22 && c11 > c33
        x = sqrt((1 + c11 - c22 - c33) / 4)
        w = (c23 - c32) / 4 / x
        y = (c12 + c21) / 4 / x
        z = (c31 + c13) / 4 / x
    elseif c22 > c33
        y = sqrt((1 - c11 + c22 - c33) / 4)
        w = (c31 - c13) / 4 / y
        x = (c12 + c21) / 4 / y
        z = (c23 + c32) / 4 / y
    else
        z = sqrt((1 - c11 - c22 + c33) / 4)
        w = (c12 - c21) / 4 / z
        x = (c31 + c13) / 4 / z
        y = (c23 + c32) / 4 / z
    end

    return Quat(w,x,y,z)
end


###########################################################################
#           Arithmetic Overloading
###########################################################################

# Operator overloading for quaternion multiplication
import Base: *
*(q1::Quat, q2::Quat) = quat_mult(q1, q2)

# Operator overloading for RK4 integration
import Base: +

"""
    +(q1::Quat, q2::Quat)

Component-wise quaternion addition (for RK4 derivatives).
Note: This is NOT a standard quaternion operation, but needed for RK4.
"""
function +(q1::Quat, q2::Quat)
    return Quat(q1.w + q2.w, q1.x + q2.x, q1.y + q2.y, q1.z + q2.z)
end

Base.:(==)(q1::Quat, q2::Quat) = isapprox(q1.w,q2.w, atol=1e-12) && isapprox(q1.x,q2.x, atol=1e-12) && isapprox(q1.y,q2.y, atol=1e-12) && isapprox(q1.z,q2.z, atol=1e-12)


"""
    *(a::Real, q::Quat)

Scalar multiplication of quaternion (for RK4 derivatives).
"""
function *(a::Real, q::Quat)
    return Quat(a * q.w, a * q.x, a * q.y, a * q.z)
end

# Also support q * a
*(q::Quat, a::Real) = a * q
