using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: propagate_keplerian, state_from_oe
@testset "Orbital Dynamics Tests" begin
    @testset "Keplerian Propagation - Circular Orbit" begin
        # Test circular orbit (e=0)
        # For a circular orbit, the true anomaly should advance uniformly

        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0  # MJD

        # Circular LEO at 400 km
        a = (6378.137 + 400) * 1e3  # [m]
        e = 0.0
        i = 0.0  # equatorial
        Ω = 0.0
        ω = 0.0
        ν0 = 0.0
        oe0 = [a, e, i, Ω, ω, ν0]

        state = state_from_oe(q0, ω0, t0, oe0)

        # Propagate one orbit
        n = sqrt(μ_earth / a^3)  # mean motion [rad/s]
        T = 2π / n  # orbital period [s]

        dt = T / 100  # time step
        n_steps = 100

        times = Float64[]
        oe_history = []
        r_eci_history = []
        v_eci_history = []

        for i in 1:n_steps
            t = (i-1) * dt
            push!(times, t)
            push!(oe_history, copy(state.orbital_elements))
            push!(r_eci_history, copy(state.r_eci))
            push!(v_eci_history, copy(state.v_eci))

            oe_new = propagate_keplerian(state, dt)
            state = state_from_oe([state.q.w, state.q.x, state.q.y, state.q.z], state.ω, state.t + dt/86400.0, oe_new)
        end

        # Record final state after all propagations
        push!(times, n_steps * dt)
        push!(oe_history, copy(state.orbital_elements))
        push!(r_eci_history, copy(state.r_eci))
        push!(v_eci_history, copy(state.v_eci))

        # Test that we completed one full orbit
        final_anomaly = oe_history[end][6]
        @test isapprox(final_anomaly, 2π, atol=0.1)

        # Plot all orbital elements
        try
            a_vals = [oe[1] for oe in oe_history]
            e_vals = [oe[2] for oe in oe_history]
            i_vals = [oe[3] for oe in oe_history]
            Ω_vals = [oe[4] for oe in oe_history]
            ω_vals = [oe[5] for oe in oe_history]
            ν_vals = [oe[6] for oe in oe_history]

            analytical_ν = n .* times

            p1 = plot(times ./ 60, a_vals, ylabel="a [m]", title="Circular Orbit (e=0)", legend=false)
            p2 = plot(times ./ 60, e_vals, ylabel="e [-]", ylims=(0, 1), legend=false)
            p3 = plot(times ./ 60, i_vals, ylabel="i [rad]", ylims=(0, 2pi), legend=false)
            p4 = plot(times ./ 60, Ω_vals, ylabel="Ω [rad]", ylims=(0, 2pi), legend=false)
            p5 = plot(times ./ 60, ω_vals, ylabel="ω [rad]", ylims=(0, 2pi), xlabel="Time [min]", legend=false)
            p6 = plot(times ./ 60, ν_vals, label="ν [rad]", ylims=(0, 2pi), lw=2, ylabel="ν [rad]", xlabel="t [min]")
            plot!(p6, times ./ 60, analytical_ν, label="Mean anomaly [rad]", lw=2, ls=:dash)

            p = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(800, 900))
            savefig(p, "plots/orbit_circular.png")
            println("✓ Circular orbit plot saved to test/plots/orbit_circular.png")

            # Plot r_eci and v_eci
            r_x = [r[1] for r in r_eci_history] ./ 1e3  # Convert to km
            r_y = [r[2] for r in r_eci_history] ./ 1e3
            r_z = [r[3] for r in r_eci_history] ./ 1e3
            r_mag = [norm(r) for r in r_eci_history] ./ 1e3

            v_x = [v[1] for v in v_eci_history] ./ 1e3  # Convert to km/s
            v_y = [v[2] for v in v_eci_history] ./ 1e3
            v_z = [v[3] for v in v_eci_history] ./ 1e3
            v_mag = [norm(v) for v in v_eci_history] ./ 1e3

            pr1 = plot(times ./ 60, r_x, label="r_x", ylabel="Position [km]", title="ECI Position (Circular)", lw=2)
            plot!(pr1, times ./ 60, r_y, label="r_y", lw=2)
            plot!(pr1, times ./ 60, r_z, label="r_z", lw=2)

            pr2 = plot(times ./ 60, r_mag, ylims=(0, maximum(r_mag) * 1.1), label="||r||", ylabel="Magnitude [km]", xlabel="Time [min]", lw=2, legend=false)

            pv1 = plot(times ./ 60, v_x, label="v_x", ylabel="Velocity [km/s]", title="ECI Velocity (Circular)", lw=2)
            plot!(pv1, times ./ 60, v_y, label="v_y", lw=2)
            plot!(pv1, times ./ 60, v_z, label="v_z", lw=2)

            pv2 = plot(times ./ 60, v_mag, ylims=(0, maximum(v_mag) * 1.1), label="||v||", ylabel="Magnitude [km/s]", xlabel="Time [min]", lw=2, legend=false)

            # 2D trajectory plot
            ptraj = plot(r_x, r_y, aspect_ratio=:equal, label="Orbit", xlabel="X [km]", ylabel="Y [km]", title="Orbital Trajectory (XY)", lw=2)
            scatter!(ptraj, [0], [0], label="Earth", markersize=8, color=:blue)

            p_state = plot(pr1, pr2, pv1, pv2, ptraj, layout=(3,2), size=(800, 900))
            savefig(p_state, "plots/orbit_circular_state.png")
            println("✓ Circular orbit state plot saved to test/plots/orbit_circular_state.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end

    @testset "Keplerian Propagation - Eccentric Orbit" begin
        # Test eccentric orbit (e=0.2)
        # True anomaly should vary non-uniformly due to Kepler's second law

        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0  # MJD

        # Eccentric orbit
        a = (6378.137 + 600) * 1e3  # [m]
        e = 0.2
        i = 0.0
        Ω = 0.0
        ω = 0.0
        ν0 = 0.0
        oe0 = [a, e, i, Ω, ω, ν0]

        state = state_from_oe(q0, ω0, t0, oe0)

        # Propagate one orbit
        n = sqrt(μ_earth / a^3)  # mean motion [rad/s]
        T = 2π / n  # orbital period [s]

        dt = T / 100  # time step
        n_steps = 100

        times = Float64[]
        oe_history = []
        r_eci_history = []
        v_eci_history = []

        for i in 1:n_steps
            t = (i-1) * dt
            push!(times, t)
            push!(oe_history, copy(state.orbital_elements))
            push!(r_eci_history, copy(state.r_eci))
            push!(v_eci_history, copy(state.v_eci))

            oe_new = propagate_keplerian(state, dt)
            state = state_from_oe([state.q.w, state.q.x, state.q.y, state.q.z], state.ω, state.t + dt/86400.0, oe_new)
        end

        # Record final state after all propagations
        push!(times, n_steps * dt)
        push!(oe_history, copy(state.orbital_elements))
        push!(r_eci_history, copy(state.r_eci))
        push!(v_eci_history, copy(state.v_eci))

        # Test that we completed one full orbit
        final_anomaly = oe_history[end][6]
        @test isapprox(final_anomaly, 2π, atol=0.1)

        # Test that other elements stayed constant
        @test state.orbital_elements[1] ≈ a
        @test state.orbital_elements[2] ≈ e
        @test state.orbital_elements[3] ≈ i
        @test state.orbital_elements[4] ≈ Ω
        @test state.orbital_elements[5] ≈ ω

        # Plot all orbital elements
        try
            a_vals = [oe[1] for oe in oe_history]
            e_vals = [oe[2] for oe in oe_history]
            i_vals = [oe[3] for oe in oe_history]
            Ω_vals = [oe[4] for oe in oe_history]
            ω_vals = [oe[5] for oe in oe_history]
            ν_vals = [oe[6] for oe in oe_history]

            mean_anomalies = n .* times

            p1 = plot(times ./ 60, a_vals, ylabel="a [m]", title="Eccentric Orbit (e=0.2)", legend=false)
            p2 = plot(times ./ 60, e_vals, ylabel="e [-]", ylims=(0, 1), legend=false)
            p3 = plot(times ./ 60, i_vals, ylabel="i [rad]", ylims=(0, 2pi), legend=false)
            p4 = plot(times ./ 60, Ω_vals, ylabel="Ω [rad]", ylims=(0, 2pi), legend=false)
            p5 = plot(times ./ 60, ω_vals, ylabel="ω [rad]", ylims=(0, 2pi), xlabel="Time [min]", legend=false)
            p6 = plot(times ./ 60, ν_vals, label="ν [rad]", ylims=(0, 2pi), lw=2, ylabel="ν [rad]", xlabel="t [min]")
            plot!(p6, times ./ 60, mean_anomalies, label="Mean anomaly [rad]", lw=2, ls=:dash)

            p = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(800, 900))
            savefig(p, "plots/orbit_eccentric.png")
            println("✓ Eccentric orbit plot saved to test/plots/orbit_eccentric.png")

            # Plot r_eci and v_eci
            r_x = [r[1] for r in r_eci_history] ./ 1e3  # Convert to km
            r_y = [r[2] for r in r_eci_history] ./ 1e3
            r_z = [r[3] for r in r_eci_history] ./ 1e3
            r_mag = [norm(r) for r in r_eci_history] ./ 1e3

            v_x = [v[1] for v in v_eci_history] ./ 1e3  # Convert to km/s
            v_y = [v[2] for v in v_eci_history] ./ 1e3
            v_z = [v[3] for v in v_eci_history] ./ 1e3
            v_mag = [norm(v) for v in v_eci_history] ./ 1e3

            pr1 = plot(times ./ 60, r_x, label="r_x", ylabel="Position [km]", title="ECI Position (Eccentric)", lw=2)
            plot!(pr1, times ./ 60, r_y, label="r_y", lw=2)
            plot!(pr1, times ./ 60, r_z, label="r_z", lw=2)

            pr2 = plot(times ./ 60, r_mag, ylims=(0, maximum(r_mag) * 1.1), label="||r||", ylabel="Magnitude [km]", xlabel="Time [min]", lw=2, legend=false)

            pv1 = plot(times ./ 60, v_x, label="v_x", ylabel="Velocity [km/s]", title="ECI Velocity (Eccentric)", lw=2)
            plot!(pv1, times ./ 60, v_y, label="v_y", lw=2)
            plot!(pv1, times ./ 60, v_z, label="v_z", lw=2)

            pv2 = plot(times ./ 60, v_mag, ylims=(0, maximum(v_mag) * 1.1), label="||v||", ylabel="Magnitude [km/s]", xlabel="Time [min]", lw=2, legend=false)

            # 2D trajectory plot - shows elliptical shape
            ptraj = plot(r_x, r_y, aspect_ratio=:equal, label="Orbit", xlabel="X [km]", ylabel="Y [km]", title="Orbital Trajectory (XY)", lw=2)
            scatter!(ptraj, [0], [0], label="Earth", markersize=8, color=:blue)

            p_state = plot(pr1, pr2, pv1, pv2, ptraj, layout=(3,2), size=(800, 900))
            savefig(p_state, "plots/orbit_eccentric_state.png")
            println("✓ Eccentric orbit state plot saved to test/plots/orbit_eccentric_state.png")
        catch e
            println("⚠ Could not generate plot: $e")
        end
    end

    @testset "Orbital Period Validation" begin
        # Verify orbital period matches theoretical value

        q0 = [1.0, 0.0, 0.0, 0.0]
        ω0 = [0.0, 0.0, 0.0]
        t0 = 60000.0

        # ISS-like orbit
        a = (6378.137 + 420) * 1e3  # [m]
        e = 0.0001
        oe0 = [a, e, 0.0, 0.0, 0.0, 0.0]

        state = state_from_oe(q0, ω0, t0, oe0)

        # Theoretical period
        n = sqrt(μ_earth / a^3)
        T_theory = 2π / n

        # Propagate for one full orbital period
        dt = 10.0  # [s]
        n_steps = ceil(Int, T_theory / dt)

        for i in 1:n_steps
            oe_new = propagate_keplerian(state, dt)
            state = state_from_oe([state.q.w, state.q.x, state.q.y, state.q.z], state.ω, state.t + dt/86400.0, oe_new)
        end

        elapsed_time = n_steps * dt

        # After one period, true anomaly should be close to initial value (modulo 2π)
        final_ν = state.orbital_elements[6]
        initial_ν = oe0[6]

        # Check that we've completed approximately one orbit
        @test isapprox(mod(final_ν, 2π), mod(initial_ν, 2π), atol=0.1)

        println("  Theoretical period: $(T_theory/60) min")
        println("  Measured period: $(elapsed_time/60) min")
    end
end

nothing