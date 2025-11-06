# Rigid Body dynamics

using LinearAlgebra

using ..Satellite: Parameters, State
using ..Constants: μ_earth


"""
    rigid_body(state, force_eci, torque_body)

Compute attitude state derivatives (quaternion and angular velocity rates).

# Arguments
- `state::State`: Current state
- `dt_s::Float64`: Timestep
- `force_eci::Vector{Float64}`: Total external forces [N·m] (ECI)
- `torque_body::Vector{Float64}`:  Total external torques in [N·m] (body frame)

# Returns
- `next_state::State`: Next state after timestep
- `ω_dot::Vector{Float64}`: Angular acceleration [rad/s^2] (body frame)

# Equations
- Quaternion kinematics: q_dot = 0.5 * Ω(ω) * q
- Euler's equation: ω_dot = I^-1 * (τ - ω × (I * ω))
"""
# function rigid_body_derivative(state::State, dt_s::Float64, force_eci::Vector{Float64}, torque_body::Vector{Float64})

# end


function rigid_body_derivative(t::Float64, state::Vector{Float64}, dt_s::Float64, force_eci::Vector{Float64}, torque_body::Vector{Float64})

end

"""

    Vector3d r = state.segment<3>(0);
    Vector3d v = state.segment<3>(3);
    Quaterniond q(state(6), state(7), state(8), state(9));
    Vector3d w = state.segment<3>(10);

    VectorXd dot(13);

    StateDerivative state_dot;

    // Position derivative is velocity 
    dot.segment<3>(0) = v;

    // Velocity derivative is acceleration (Schaub 2.15)
    dot.segment<3>(3) = this->force / this->mass_kg;

    // Quaternion derivative is based on hamilton prodict (Schaub 3.109)
    Quaterniond product = Cesium::Math::hamilton_product(this->state.q, this->state.w);
    dot(6) = product.w();
    dot(7) = product.x();
    dot(8) = product.y();
    dot(9) = product.z();

    // Angular velocity derivative based on (Schaub 4.34)
    // dwdt = I.inv * (τ - w x (Iw) )
    dot.segment<3>(10) = this->InertiaTensor.inverse() * (this->torque - this->state.w.cross(this->InertiaTensor * this->state.w));

    return dot;


"""