# Interactive satellite visualization using GLMakie
# This provides a simple interface for visualizing satellite trajectories

using ..Satellite: State
using ..Math: q_to_dcm
using GLMakie
using LinearAlgebra

# Import constants
include("../constants.jl")

"""
    visualize(trajectory::Vector{State}; interactive=true)

Visualize satellite trajectory with both orbital and attitude dynamics.

# Arguments
- `trajectory::Vector{State}`: Array of states from propagation

# Keyword Arguments
- `interactive::Bool`: If true, creates interactive plot with slider and animation controls (default: true)

# Returns
- Figure handle

# Example
```julia
using Simwise

# Run your simulation
trajectory = propagate(initial_state, dt, n_steps, params)

# Visualize interactively
Simwise.Visualization.visualize(trajectory)
```
"""
function visualize(trajectory::Vector{State}; interactive=true)
    if interactive
        return visualize_interactive(trajectory)
    else
        return visualize_static(trajectory)
    end
end

"""
    visualize_interactive(trajectory::Vector{State})

Create interactive 3D visualization with time slider and animation controls.
"""
function visualize_interactive(trajectory::Vector{State})
    # Create figure with layout
    fig = Figure(size=(1800, 1000))

    # Time slider
    time_idx = Observable(1)
    max_idx = length(trajectory)

    sg = SliderGrid(
        fig[2, :],
        (label="Time Step", range=1:max_idx, startvalue=1),
        height=60
    )

    time_slider = sg.sliders[1]
    connect!(time_idx, time_slider.value)

    # Extract position data
    x_pos = [state.r_eci[1]/1e3 for state in trajectory]
    y_pos = [state.r_eci[2]/1e3 for state in trajectory]
    z_pos = [state.r_eci[3]/1e3 for state in trajectory]

    # ===== PLOT 1: ORBITAL TRAJECTORY =====
    ax_orbit = Axis3(fig[1, 1],
        title="Orbital Trajectory",
        xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)",
        aspect=:equal,
        width=700,
        height=600)

    # Earth
    earth_radius = R_earth / 1e3
    meshscatter!(ax_orbit, [Point3f(0, 0, 0)],
                 markersize=earth_radius, color=:lightblue)

    # Full orbit path (faint)
    lines!(ax_orbit, x_pos, y_pos, z_pos, color=(:gray, 0.3), linewidth=1)

    # Observable orbit path
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

    # Time label
    time_text_orbit = @lift("t = $(round(trajectory[$time_idx].t * 86400, digits=1)) s")
    text!(ax_orbit, 0, 0, maximum(z_pos) * 1.2,
          text=time_text_orbit, fontsize=20, align=(:center, :center))

    # ===== PLOT 2: ATTITUDE VISUALIZATION =====
    ax_attitude = Axis3(fig[1, 2],
        title="Attitude (Body Frame in ECI)",
        xlabel="X", ylabel="Y", zlabel="Z",
        aspect=:equal,
        width=700,
        height=600)

    # ECI reference frame (fixed)
    eci_scale = 1.5
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(eci_scale, 0, 0)],
           color=(:gray, 0.5), linewidth=3)
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(0, eci_scale, 0)],
           color=(:gray, 0.5), linewidth=3)
    lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(0, 0, eci_scale)],
           color=(:gray, 0.5), linewidth=3)

    # Body frame axes (update with slider)
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

    lines!(ax_attitude, body_x_line, color=:red, linewidth=5)
    lines!(ax_attitude, body_y_line, color=:green, linewidth=5)
    lines!(ax_attitude, body_z_line, color=:blue, linewidth=5)

    # Time label
    time_text_att = @lift("t = $(round(trajectory[$time_idx].t * 86400, digits=1)) s")
    text!(ax_attitude, 0, 0, 2.0,
          text=time_text_att, fontsize=20, align=(:center, :center))

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

    # ===== CONTROLS =====
    is_playing = Observable(false)
    play_button = Button(fig[4, 1], label=@lift($is_playing ? "Pause" : "Play"), fontsize=16, height=40)

    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
    end

    # Animation loop
    @async while true
        if is_playing[] && time_idx[] < max_idx
            time_idx[] = time_idx[] + 1
            sleep(0.03)
        elseif time_idx[] >= max_idx
            is_playing[] = false
            time_idx[] = 1
        end
        sleep(0.01)
    end

    reset_button = Button(fig[4, 2], label="Reset", fontsize=16, height=40)
    on(reset_button.clicks) do _
        time_idx[] = 1
        is_playing[] = false
    end

    # Layout adjustments
    rowsize!(fig.layout, 1, Relative(0.75))
    rowsize!(fig.layout, 2, 60)
    rowsize!(fig.layout, 3, Auto())
    rowsize!(fig.layout, 4, 50)
    rowgap!(fig.layout, 10)

    colsize!(fig.layout, 1, Relative(0.5))
    colsize!(fig.layout, 2, Relative(0.5))
    colgap!(fig.layout, 20)

    display(fig)
    return fig
end

"""
    visualize_static(trajectory::Vector{State})

Create static visualization showing full trajectory.
"""
function visualize_static(trajectory::Vector{State})
    fig = Figure(size=(1600, 800))

    # Extract data
    x_pos = [state.r_eci[1]/1e3 for state in trajectory]
    y_pos = [state.r_eci[2]/1e3 for state in trajectory]
    z_pos = [state.r_eci[3]/1e3 for state in trajectory]

    # Orbit plot
    ax_orbit = Axis3(fig[1, 1],
        title="Orbital Trajectory",
        xlabel="X (km)", ylabel="Y (km)", zlabel="Z (km)",
        aspect=:equal)

    earth_radius = R_earth / 1e3
    meshscatter!(ax_orbit, [Point3f(0, 0, 0)],
                 markersize=earth_radius, color=:lightblue)

    lines!(ax_orbit, x_pos, y_pos, z_pos, color=:cyan, linewidth=2)
    scatter!(ax_orbit, [x_pos[1]], [y_pos[1]], [z_pos[1]], color=:green, markersize=10)
    scatter!(ax_orbit, [x_pos[end]], [y_pos[end]], [z_pos[end]], color=:red, markersize=10)

    # Attitude snapshots at key times
    ax_attitude = Axis3(fig[1, 2],
        title="Attitude Evolution",
        xlabel="X", ylabel="Y", zlabel="Z",
        aspect=:equal)

    # Show initial and final orientations
    for (idx, color_alpha) in [(1, 0.3), (length(trajectory), 1.0)]
        dcm = q_to_dcm(trajectory[idx].q)'
        for (col, color) in enumerate([:red, :green, :blue])
            axis = dcm[:, col]
            lines!(ax_attitude, [Point3f(0, 0, 0), Point3f(axis)],
                   color=(color, color_alpha), linewidth=3)
        end
    end

    limits!(ax_attitude, -1.5, 1.5, -1.5, 1.5, -1.5, 1.5)

    display(fig)
    return fig
end

export visualize
