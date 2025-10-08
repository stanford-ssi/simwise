# Atmospheric density model

"""
    atmospheric_density(altitude)

Compute atmospheric density at given altitude using exponential model.

# Arguments
- `altitude::Float64`: Altitude above Earth surface [m]

# Returns
- `Float64`: Atmospheric density [kg/m^3]
"""
function atmospheric_density(altitude::Float64)
    # TODO: Implement atmospheric density model (exponential or NRLMSISE-00)
    # ρ = ρ0 * exp(-h/H)
end
