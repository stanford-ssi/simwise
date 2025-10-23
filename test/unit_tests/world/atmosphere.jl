using Plots
using LinearAlgebra

using Simwise.World: low_solar_activity, moderate_solar_activity, high_solar_activity

# println(typeof(low_solar_activity))
println()
# println([v[1] for (_, v) in high_solar_activity])
p = plot([k for (k, _) in high_solar_activity], [v["temp"] for (_, v) in high_solar_activity],
        title="Temperature vs Altitude", xlabel="Altitude (km)", ylabel="Temperature (K)")
savefig(p, "src/world/atmosphere_graphs/temp.png")
p = plot([k for (k, _) in high_solar_activity], [log(v["density"]) for (_, v) in high_solar_activity],
        title="Density vs Altitude", xlabel="Altitude (km)", ylabel="log(density)")
savefig(p, "src/world/atmosphere_graphs/density.png")
p = plot([k for (k, _) in high_solar_activity], [log(v["pressure"]) for (_, v) in high_solar_activity],
        title="Pressure vs Altitude", xlabel="Altitude (km)", ylabel="log(pressure)")
savefig(p, "src/world/atmosphere_graphs/pressure.png")
p = plot([k for (k, _) in high_solar_activity], [v["mol_wt"] for (_, v) in high_solar_activity],
        title="Mole Weight vs Altitude", xlabel="Altitude (km)", ylabel="mol_wt (kg/kmol)")
savefig(p, "src/world/atmosphere_graphs/mol_wt.png")
