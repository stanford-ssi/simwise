using Plots

using Simwise.Satellite: State, Parameters
using Simwise.Dynamics: state_from_oe
using Simwise.Math: rk4_step, Quat

@testset "Propagator Tests" begin
    
    """Defining derivative functions for RK4 tests"""
    function constant_deriv(t::Float64, x::Vector{Float64}, dummy_params)
        return zeros(Float64, length(x))
    end

    function slope_deriv(t::Float64, x::Vector{Float64}, dummy_params)
        return ones(Float64, length(x)) * 2
    end

    function exponential_deriv(t::Float64, x::Float64, dummy_params)
        return -x
    end

    function sine_deriv(t::Float64, x::Float64, dummy_params)
        return cos(t)
    end

    function logistic_deriv_3(t::Float64, x::Vector{Float64}, dummy_params)
        r = 1.0
        K = 10.0
        return r * x .* (ones(length(x)) - x/K)
    end

    dummy_params = 0.0

    @testset "RK4 Constant" begin
        vec = [1.0 ,0.0 ,0.0 ,2.0]

        # random values for dt and t
        result = rk4_step(1.0, 1.0, vec, constant_deriv, dummy_params)
        
        # Should be equal to original vector since slope of 0
        @test result == vec 
    end

    @testset "RK4 Constant Slope of 2" begin
        vec = [1.0, 0.0, 0.0, 2.0]

        # random values for dt and t
        result = rk4_step(2.0, 0.0, vec, slope_deriv, dummy_params)
        
        # Since dx/dt = 2, exact solution e.g.: x = 0.0 + 2*2 = 4.0
        @test result == [5.0, 4.0, 4.0, 6.0]
    end

    @testset "RK4 Exponential Decay Single Step" begin
        x = 4.0

        dt = 0.2

        # random values for dt and t
        result = rk4_step(dt, 0.0, x, exponential_deriv, dummy_params)
        
        # Exact solution: x = 4.0 * e^(-tmax) since x_dot = -x
        @test isapprox(result, 4*ℯ^(-dt), atol=1e-3)
    end

    @testset "RK4 Exponential Decay Multi Step" begin
        x = 4.0

        dt = 0.001;
        tmax = 10.0;  # Simulates from t = 0 to t = 1.0
        t = 0.0;

        while t < tmax
            x = rk4_step(dt, t, x, exponential_deriv, dummy_params)
            t += dt
        end

        # Exact solution: x = 4.0 * e^(-tmax) since x_dot = -x
        @test isapprox(x, 4*ℯ^(-tmax), rtol=1e-2)
    end

    @testset "RK4 Sine Single Step" begin
        x = 4.0

        dt = 0.01;

        result = rk4_step(dt, 0.0, x, sine_deriv, dummy_params)

        # Simply finding the result at the new time
        @test isapprox(result, 4.0 + sin(dt), rtol=1e-3)
    end

    @testset "RK4 Sine Multi Step 1 Second" begin
        x0 = 5.0

        dt = 0.01;
        tmax = 1.0;  # Simulates from t = 0 to t = 10.0
        t = 0.0;
        x = x0

        while t < tmax
            x = rk4_step(dt, t, x, sine_deriv, dummy_params)
            t += dt
        end

        # Simply finding the result at the new time
        @test isapprox(x, x0 + sin(tmax), rtol=1e-3)
    end

    @testset "RK4 Sine Multi Step 15 Second" begin
        x0 = 5.0

        dt = 0.005;
        tmax = 15.0;  # Simulates from t = 0 to t = 10.0
        t = 0.0;
        x = x0

        while t < tmax
            x = rk4_step(dt, t, x, sine_deriv, dummy_params)
            t += dt
        end

        # Simply finding the result at the new time
        @test isapprox(x, x0 + sin(tmax), rtol=1e-3)
    end

    @testset "RK4 Logistic Growth" begin
        x0 = [1.0, 2.0, 5.0]

        dt = 0.001;
        tmax = 15.0;  # Simulates from t = 0 to t = 10.0
        t = 0.0;
        x = x0

        while t < tmax
            x = rk4_step(dt, t, x, logistic_deriv_3, dummy_params)
            t += dt
        end

        # Should be near (10, 10, 10) eventually at fairly aggressive logistic growth
        @test isapprox(x, [10,10,10], atol=1e-4)
    end

    @testset "RK4 State Integration" begin
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
