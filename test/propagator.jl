using Plots

@testset "Propagator Tests" begin
    @testset "RK4 Integration" begin
        # Test with simple exponential: dx/dt = -k*x (solution: x(t) = x0*exp(-k*t))
        # Using a simple State with just one element varying

        # Simple derivative function: exponential decay in first orbital element
        function simple_dynamics(state::State, params)
            k = 0.1  # decay constant
            dq = zeros(4)
            dω = zeros(3)
            dt_mjd = 1.0 / 86400.0  # 1 second in MJD units
            doe = [-k * state.orbital_elements[1], 0.0, 0.0, 0.0, 0.0, 0.0]
            return State(dq, dω, dt_mjd, doe)
        end

        # Initial state
        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0  # MJD
        oe0 = [7000e3, 0.0, 0.0, 0.0, 0.0, 0.0]  # only first element varies
        state0 = State(q0, ω0, t0, oe0)

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

            state = rk4_step(simple_dynamics, state, dt, params)
        end

        # Test accuracy at final time (relaxed tolerance for multi-step)
        final_numerical = numerical[end]
        final_analytical = analytical[end]
        @test isapprox(final_numerical, final_analytical, rtol=1e-5)

        # Test single step accuracy
        state1 = rk4_step(simple_dynamics, state0, dt, params)
        a_expected = oe0[1] * exp(-k * dt)
        @test isapprox(state1.orbital_elements[1], a_expected, rtol=1e-6)
        @test isapprox(state1.t, t0 + dt/86400.0, rtol=1e-10)

        # Plot comparison for visual inspection
        try
            p = plot(times, analytical, label="Analytical", lw=2, xlabel="Time [s]", ylabel="Semi-major axis [m]")
            plot!(p, times, numerical, label="RK4", lw=2, ls=:dash)
            savefig(p, "rk4_validation.png")
            println("✓ RK4 validation plot saved to test/rk4_validation.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end
end
