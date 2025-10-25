using Plots

using Simwise.Satellite: State, Parameters
using Simwise.Dynamics: state_from_oe
using Simwise.Math: rk4_step, Quat

@testset "Propagator Tests" begin
    @testset "RK4 Integration" begin
        # Test with simple exponential: dx/dt = -k*x (solution: x(t) = x0*exp(-k*t))
        # Using a simple State with just one element varying

        # Simple derivative function: exponential decay in first orbital element
        function simple_dynamics(t::Float64, state::State, params)
            k = 0.1  # decay constant
            dq = Quat(0.0, 0.0, 0.0, 0.0)
            dω = zeros(3)
            dt_mjd = 1.0 / 86400.0  # 1 second in MJD units
            doe = [-k * state.orbital_elements[1], 0.0, 0.0, 0.0, 0.0, 0.0]
            dr_eci = zeros(3)
            dv_eci = zeros(3)
            return State(dq, dω, dt_mjd, doe, dr_eci, dv_eci)
        end

        # Initial state
        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0  # MJD
        oe0 = [7000e3, 0.0, 0.0, 0.0, 0.0, 0.0]  # only first element varies
        state0 = state_from_oe(q0, ω0, t0, oe0)

        # Dummy parameters
        params = Parameters(10.0, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], 0.01, 2.2)

        # Propagate over time
        dt = 1.0  # 1 second
        n_steps = 100
        k = 0.1

        times = Float64[]
        numerical = Float64[]
        analytical = Float64[]

        state = state0
        for i in 1:n_steps
            t = (i-1) * dt
            push!(times, t)
            push!(numerical, state.orbital_elements[1])
            push!(analytical, oe0[1] * exp(-k * t))

            state = rk4_step(dt, t, state, simple_dynamics, params)
        end

        # Test accuracy at final time (relaxed tolerance for multi-step)
        final_numerical = numerical[end]
        final_analytical = analytical[end]
        @test isapprox(final_numerical, final_analytical, rtol=1e-5)

        # Test single step accuracy
        state1 = rk4_step(dt, 0.0, state0, simple_dynamics, params)
        a_expected = oe0[1] * exp(-k * dt)
        @test isapprox(state1.orbital_elements[1], a_expected, rtol=1e-6)
        @test isapprox(state1.t, t0 + dt/86400.0, rtol=1e-10)

        # Plot comparison for visual inspection
        try
            p = plot(times, analytical, label="Analytical", lw=2, xlabel="Time [s]", ylabel="Semi-major axis [m]", title="Exponential Decay")
            plot!(p, times, numerical, label="RK4", lw=2, ls=:dash)
            savefig(p, "plots/rk4_validation_exponential.png")
            println("✓ RK4 exponential validation plot saved to test/plots/rk4_validation_exponential.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end

    @testset "RK4 Oscillatory Test" begin
        # Test with oscillatory function: dx/dt = A*omega*cos(omega*t)
        # Solution: x(t) = x0 + A*sin(omega*t)

        function oscillatory_dynamics(t::Float64, state::State, params)
            A = 1000e3  # amplitude [m]
            omega = 0.1  # angular frequency [rad/s]

            # Extract current time from state
            t = (state.t - 60000.0) * 86400.0  # convert MJD to seconds from t0

            dq = Quat(0.0, 0.0, 0.0, 0.0)
            dω = zeros(3)
            dt_mjd = 1.0 / 86400.0  # 1 second in MJD units
            doe = [A * omega * cos(omega * t), 0.0, 0.0, 0.0, 0.0, 0.0]
            dr_eci = zeros(3)
            dv_eci = zeros(3)
            return State(dq, dω, dt_mjd, doe, dr_eci, dv_eci)
        end

        # Initial state
        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0  # MJD
        a0 = 7000e3  # initial semi-major axis [m]
        oe0 = [a0, 0.0, 0.0, 0.0, 0.0, 0.0]
        state0 = state_from_oe(q0, ω0, t0, oe0)

        # Dummy parameters
        params = Parameters(10.0, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], 0.01, 2.2)

        # Propagate over time
        dt = 1.0  # 1 second
        n_steps = 100
        A = 1000e3
        omega = 0.1

        times = Float64[]
        numerical = Float64[]
        analytical = Float64[]

        state = state0
        for i in 1:n_steps
            t = (i-1) * dt
            push!(times, t)
            push!(numerical, state.orbital_elements[1])
            push!(analytical, a0 + A * sin(omega * t))

            state = rk4_step(dt, t, state, oscillatory_dynamics, params)
        end

        # Test accuracy at final time
        final_numerical = numerical[end]
        final_analytical = analytical[end]
        @test isapprox(final_numerical, final_analytical, rtol=1e-5)

        # Plot comparison for visual inspection
        try
            p = plot(times, analytical, label="Analytical", lw=2, xlabel="Time [s]", ylabel="Semi-major axis [m]", title="Sinusoidal Oscillation")
            plot!(p, times, numerical, label="RK4", lw=2, ls=:dash)
            savefig(p, "plots/rk4_validation_oscillatory.png")
            println("✓ RK4 oscillatory validation plot saved to test/plots/rk4_validation_oscillatory.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end
end

nothing