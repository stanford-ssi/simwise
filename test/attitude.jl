using Plots
using LinearAlgebra

@testset "Attitude Dynamics Tests" begin
    @testset "Torque-Free Rigid Body" begin
        # Test torque-free rigid body motion (conservation of angular momentum)
        # For a symmetric rigid body with no external torques, ω should remain constant

        # Initial state: spinning about z-axis
        q0 = [1.0, 0.0, 0.0, 0.0]  # Identity quaternion
        ω0 = [0.0, 0.0, 0.1]  # 0.1 rad/s about z-axis
        t0 = 60000.0  # MJD
        oe0 = [7000e3, 0.0, 0.0, 0.0, 0.0, 0.0]
        state0 = State(q0, ω0, t0, oe0)

        # Spherical inertia (symmetric)
        I = diagm([1.0, 1.0, 1.0])  # kg·m²

        # No external torques
        torques = zeros(3)

        # Compute dynamics
        q_dot, ω_dot = attitude_dynamics(state0, torques, I)

        # For symmetric inertia with no torques, ω_dot should be zero
        @test isapprox(ω_dot, zeros(3), atol=1e-12)

        # Quaternion derivative should follow q_dot = 0.5 * Ω(ω) * q
        # For q = [1, 0, 0, 0] and ω = [0, 0, 0.1]:
        # q_dot = 0.5 * [0, 0, 0, 0.1]
        expected_q_dot = [0.0, 0.0, 0.0, 0.05]
        @test isapprox([q_dot.q0, q_dot.q1, q_dot.q2, q_dot.q3], expected_q_dot, atol=1e-12)
    end

    @testset "Asymmetric Body (Gyroscopic Effect)" begin
        # Test asymmetric body with gyroscopic torque
        # Angular momentum should be conserved when no external torques

        # Initial state: spinning about principal axis
        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.5, 0.0, 0.0]  # Spin about x-axis
        t0 = 60000.0
        oe0 = [7000e3, 0.0, 0.0, 0.0, 0.0, 0.0]
        state0 = State(q0, ω0, t0, oe0)

        # Asymmetric inertia (like a satellite with extended solar panels)
        I = diagm([0.5, 2.0, 2.0])  # kg·m²

        # No external torques
        torques = zeros(3)

        # Compute dynamics
        q_dot, ω_dot = attitude_dynamics(state0, torques, I)

        # Angular momentum H = I * ω
        H = I * ω0

        # For no external torques, H_dot = 0
        # H_dot = I * ω_dot + ω × (I * ω) = 0
        H_dot = I * ω_dot + cross(ω0, I * ω0)
        @test isapprox(H_dot, zeros(3), atol=1e-12)
    end

    @testset "Propagation and Quaternion Evolution" begin
        # Propagate a tumbling satellite and plot quaternion components

        # Initial state: tumbling with angular velocity on all axes
        q0 = [1.0, 0.0, 0.0, 0.0]  # Start at identity
        ω0 = [0.05, 0.02, 0.01]  # rad/s - gentle tumbling motion
        t0 = 60000.0
        oe0 = [7000e3, 0.0, 0.0, 0.0, 0.0, 0.0]
        state0 = State(q0, ω0, t0, oe0)

        # CubeSat-like inertia (slightly asymmetric)
        I = diagm([0.01, 0.012, 0.015])  # kg·m²

        # Propagate with RK4 - one full rotation at ω=0.05 rad/s takes ~126 seconds
        dt = 0.1  # 0.1 second time step
        n_steps = 1260  # 126 seconds ≈ 1 rotation

        # Storage for plotting
        times = Float64[]
        q0_hist = Float64[]
        q1_hist = Float64[]
        q2_hist = Float64[]
        q3_hist = Float64[]
        ωx_hist = Float64[]
        ωy_hist = Float64[]
        ωz_hist = Float64[]

        # Dynamics function for RK4
        function attitude_only_dynamics(state::State, params)
            torques = zeros(3)  # No external torques
            q_dot, ω_dot = attitude_dynamics(state, torques, I)

            # For this test, we don't care about orbital dynamics
            doe = zeros(6)
            # Time derivative: d(MJD)/d(second) = 1/86400
            dt_mjd = 1.0 / 86400.0

            return State(q_dot, ω_dot, dt_mjd, doe)
        end

        # Dummy parameters (not used)
        params = Parameters(1.0, I, 0.01, 2.2)

        # Propagate
        state = state0
        for i in 1:n_steps
            t = (i-1) * dt
            push!(times, t)
            push!(q0_hist, state.q.q0)
            push!(q1_hist, state.q.q1)
            push!(q2_hist, state.q.q2)
            push!(q3_hist, state.q.q3)
            push!(ωx_hist, state.ω[1])
            push!(ωy_hist, state.ω[2])
            push!(ωz_hist, state.ω[3])

            state = rk4_step(attitude_only_dynamics, state, dt, params)
        end

        # Check quaternion norm stays at 1 (should be preserved with normalization)
        for i in 1:length(times)
            q_norm = sqrt(q0_hist[i]^2 + q1_hist[i]^2 + q2_hist[i]^2 + q3_hist[i]^2)
            @test isapprox(q_norm, 1.0, atol=1e-10)
        end

        # Check angular momentum magnitude conservation (no external torques)
        # Direction changes in body frame (torque-free precession), but magnitude is conserved
        H0 = I * ω0
        Hf = I * [ωx_hist[end], ωy_hist[end], ωz_hist[end]]
        H0_mag = norm(H0)
        Hf_mag = norm(Hf)
        @test isapprox(H0_mag, Hf_mag, rtol=1e-6)

        # Compute angular momentum magnitude over time
        H_mag_hist = Float64[]
        for i in 1:length(times)
            H = I * [ωx_hist[i], ωy_hist[i], ωz_hist[i]]
            push!(H_mag_hist, norm(H))
        end

        # Convert quaternions to Euler angles (roll, pitch, yaw) in degrees
        roll_hist = Float64[]
        pitch_hist = Float64[]
        yaw_hist = Float64[]
        for i in 1:length(times)
            q0, q1, q2, q3 = q0_hist[i], q1_hist[i], q2_hist[i], q3_hist[i]

            # Euler angles from quaternion (ZYX convention: yaw-pitch-roll)
            roll = atan(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2))
            pitch = asin(2*(q0*q2 - q3*q1))
            yaw = atan(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2))

            push!(roll_hist, rad2deg(roll))
            push!(pitch_hist, rad2deg(pitch))
            push!(yaw_hist, rad2deg(yaw))
        end

        # Plot combined figure with 4 subplots
        try
            p1 = plot(times, q0_hist, label="q0 (scalar)", lw=2, xlabel="", ylabel="Quaternion", title="Quaternion Evolution", legend=:right)
            plot!(p1, times, q1_hist, label="q1", lw=2)
            plot!(p1, times, q2_hist, label="q2", lw=2)
            plot!(p1, times, q3_hist, label="q3", lw=2)

            p2 = plot(times, roll_hist, label="Roll", lw=2, xlabel="", ylabel="Euler Angles [deg]", title="Euler Angles (ZYX)", legend=:right)
            plot!(p2, times, pitch_hist, label="Pitch", lw=2)
            plot!(p2, times, yaw_hist, label="Yaw", lw=2)

            p3 = plot(times, ωx_hist, label="ωx", lw=2, xlabel="", ylabel="ω [rad/s]", title="Angular Velocity", legend=:right)
            plot!(p3, times, ωy_hist, label="ωy", lw=2)
            plot!(p3, times, ωz_hist, label="ωz", lw=2)

            p4 = plot(times, H_mag_hist, label="", lw=2, ylim=(0, 1.1 * maximum(H_mag_hist)), xlabel="Time [s]", ylabel="|H| [kg·m²/s]", title="Angular Momentum Magnitude", color=:black)

            combined = plot(p1, p2, p3, p4, layout=(4,1), size=(800, 1000))
            savefig(combined, "plots/attitude_dynamics.png")
            println("✓ Attitude dynamics plot saved to test/plots/attitude_dynamics.png")
        catch e
            println("⚠ Could not generate plots: $e")
        end
    end
end
