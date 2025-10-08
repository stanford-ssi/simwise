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

        # Take one RK4 step
        dt = 1.0  # 1 second
        state1 = rk4_step(simple_dynamics, state0, dt, params)

        # Analytical solution: a(t) = a0 * exp(-k*t)
        k = 0.1
        a_expected = oe0[1] * exp(-k * dt)

        # Check if RK4 result is close to analytical solution (within 1e-6 relative error)
        @test isapprox(state1.orbital_elements[1], a_expected, rtol=1e-6)

        # Check that time advanced correctly
        @test isapprox(state1.t, t0 + dt/86400.0, rtol=1e-10)
    end
end
