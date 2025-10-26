# Interactive Satellite State Visualization with GLMakie
# This script creates interactive 3D visualizations with animation and time sliders

using Simwise
using Simwise.Satellite: State, Parameters
using Simwise.Math: q_to_dcm
using Simwise.Dynamics: state_from_oe
using Simwise.Simulation: propagate
using Simwise: R_earth
using GLMakie
using LinearAlgebra

"""
    visualize_satellite_interactive(trajectory::Vector{State})

Create interactive 3D visualizations with:
1. Animated orbit around Earth
2. Animated attitude (body frame rotation)
3. Time slider to scrub through the simulation

# Arguments
- `trajectory::Vector{State}`: Array of states from propagation
"""
function visualize_satellite_interactive(trajectory::Vector{State})
    # Create figure with layout - bigger window
    fig = Figure(size=(1000, 600))

    # Time slider - make it taller
    time_idx = Observable(1)
    max_idx = length(trajectory)

    sg = SliderGrid(
        fig[2, :],
        (label="Time Step", range=1:max_idx, startvalue=1),
        height=60  # Make slider taller
    )

    # Connect slider to time index
    time_slider = sg.sliders[1]
    connect!(time_idx, time_slider.value)

    # Extract position data
    x_pos = [state.r_eci[1]/1e3 for state in trajectory]  # Convert to km
    y_pos = [state.r_eci[2]/1e3 for state in trajectory]
    z_pos = [state.r_eci[3]/1e3 for state in trajectory]

    # ===== PLOT 1: ORBITAL TRAJECTORY =====
    ax_orbit = Axis3(fig[1, 1],
        title="Orbital Trajectory",
        xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)",
        aspect=:equal,
        width=400,
        height=300)

    # Draw Earth as a simple sphere using meshscatter
    earth_radius = R_earth / 1e3  # km
    meshscatter!(ax_orbit, [Point3f(0, 0, 0)],
                 markersize=earth_radius, color=:lightblue)

    # Full orbit path (faint)
    lines!(ax_orbit, x_pos, y_pos, z_pos, color=(:gray, 0.3), linewidth=1)

    # Observable orbit path (updates with slider)
    orbit_path_x = @lift(x_pos[1:$time_idx])
    orbit_path_y = @lift(y_pos[1:$time_idx])
    orbit_path_z = @lift(z_pos[1:$time_idx])

    lines!(ax_orbit, orbit_path_x, orbit_path_y, orbit_path_z,
           color=:cyan, linewidth=3)

    # Current satellite position
    sat_x = @lift([x_pos[$time_idx]])
    sat_y = @lift([y_pos[$time_idx]])
    sat_z = @lift([z_pos[$time_idx]])

    scatter!(ax_orbit, sat_x, sat_y, sat_z,
             color=:red, markersize=20)

    # Time label for orbit
    time_text_orbit = @lift("t = $(round(trajectory[$time_idx].t * 86400, digits=1)) s")
    text!(ax_orbit, 0, 0, maximum(z_pos) * 1.2,
          text=time_text_orbit, fontsize=20, align=(:center, :center))

    # ===== PLOT 2: ATTITUDE VISUALIZATION =====
    ax_attitude = Axis3(fig[1, 2],
        title="Attitude (Body Frame in ECI)",
        xlabel="X", ylabel="Y", zlabel="Z",
        aspect=:equal,
        width=400,
        height=300)

    # ECI reference frame (fixed) - draw as lines instead of arrows for simplicity
    eci_scale = 1.5
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(eci_scale, 0, 0)],
           color=(:gray, 0.5), linewidth=3)
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(0, eci_scale, 0)],
           color=(:gray, 0.5), linewidth=3)
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(0, 0, eci_scale)],
           color=(:gray, 0.5), linewidth=3)

    # Body frame axes (update with slider) - also use lines
    body_x_line = @lift begin
        dcm = q_to_dcm(trajectory[$time_idx].q)'
        axis = dcm[:, 1]
        [Point3f(0, 0, 0), Point3f(axis)]
    end

    body_y_line = @lift begin
        dcm = q_to_dcm(trajectory[$time_idx].q)'
        axis = dcm[:, 2]
        [Point3f(0, 0, 0), Point3f(axis)]
    end

    body_z_line = @lift begin
        dcm = q_to_dcm(trajectory[$time_idx].q)'
        axis = dcm[:, 3]
        [Point3f(0, 0, 0), Point3f(axis)]
    end

    # Plot body axes as thick colored lines
    lines!(ax_attitude, body_x_line, color=:red, linewidth=5)
    lines!(ax_attitude, body_y_line, color=:green, linewidth=5)
    lines!(ax_attitude, body_z_line, color=:blue, linewidth=5)

    # Time label for attitude
    time_text_att = @lift("t = $(round(trajectory[$time_idx].t * 86400, digits=1)) s")
    text!(ax_attitude, 0, 0, 2.0,
          text=time_text_att, fontsize=20, align=(:center, :center))

    # Set axis limits
    limits!(ax_attitude, -2, 2, -2, 2, -2, 2)

    # ===== INFO PANEL =====
    info_text = @lift begin
        state = trajectory[$time_idx]
        q = state.q
        ω = state.ω
        r = state.r_eci
        v = state.v_eci

        """
        Quaternion: [$(round(q.w, digits=3)), $(round(q.x, digits=3)), $(round(q.y, digits=3)), $(round(q.z, digits=3))]
        Angular Velocity: [$(round(ω[1], digits=3)), $(round(ω[2], digits=3)), $(round(ω[3], digits=3))] rad/s
        Position: [$(round(r[1]/1e3, digits=1)), $(round(r[2]/1e3, digits=1)), $(round(r[3]/1e3, digits=1))] km
        Velocity: [$(round(v[1]/1e3, digits=2)), $(round(v[2]/1e3, digits=2)), $(round(v[3]/1e3, digits=2))] km/s
        """
    end

    Label(fig[3, :], info_text, fontsize=16, tellwidth=false, padding=(5, 5, 5, 5))

    # Play button for animation
    is_playing = Observable(false)
    play_button = Button(fig[4, 1], label=@lift($is_playing ? "Pause" : "Play"), fontsize=16, height=40)

    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
    end

    # Animation loop
    @async while true
        if is_playing[] && time_idx[] < max_idx
            time_idx[] = time_idx[] + 1
            sleep(0.03)  # ~30 fps
        elseif time_idx[] >= max_idx
            is_playing[] = false
            time_idx[] = 1  # Loop back
        end
        sleep(0.01)
    end

    # Reset button
    reset_button = Button(fig[4, 2], label="Reset", fontsize=16, height=40)
    on(reset_button.clicks) do _
        time_idx[] = 1
        is_playing[] = false
    end

    # Adjust row and column sizes for better spacing - give most space to plots!
    rowsize!(fig.layout, 1, Relative(0.75))  # Plots take 75% of space
    rowsize!(fig.layout, 2, 60)  # Slider fixed height
    rowsize!(fig.layout, 3, Auto())  # Info text auto
    rowsize!(fig.layout, 4, 50)  # Buttons fixed height
    rowgap!(fig.layout, 10)  # Smaller gaps

    # Make columns equal width
    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    colgap!(fig.layout, 20)

    display(fig)
    return fig
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running interactive satellite visualization...")

    # Define spacecraft parameters
    params = Parameters(
        10.0,                           # mass (kg)
        diagm([0.1, 0.1, 0.05]),       # inertia (kg·m^2)
        0.01,                           # area (m^2)
        2.2                             # drag coefficient
    )

    # Initial orbital elements (circular LEO at 400 km)
    a = R_earth + 400e3     # semi-major axis (m)
    e = 0.0                 # eccentricity
    i = 0.0                 # inclination (rad)
    Ω = 0.0                 # RAAN (rad)
    ω = 0.0                 # argument of periapsis (rad)
    ν = 0.0                 # true anomaly (rad)
    oe0 = [a, e, i, Ω, ω, ν]

    # Initial attitude (small rotation + tumbling)
    q0 = [1.0, 0.0, 0.0, 0.0]  # Identity quaternion (scalar-first)
    ω0 = [0.1, 0.05, 0.15]     # rad/s - tumbling motion

    # Create initial state
    state0 = state_from_oe(q0, ω0, 0.0, oe0)

    # Simulation parameters
    dt = 0.1          # s
    n_steps = 54000    # ~ 1 orbit for 400 km LEO

    println("Propagating satellite state...")
    trajectory = State[]
    push!(trajectory, state0)

    for i in 1:n_steps
        new_state = propagate(trajectory[end], dt, 1, params)
        push!(trajectory, new_state)
    end

    println("Creating interactive visualization...")
    println("Use the slider to scrub through time, or press Play to animate!")
    visualize_satellite_interactive(trajectory)

    println("Visualization ready! Close the window when done.")
    wait()
end
