module World

using LinearAlgebra

# Call every non-exported name with Simwise.World.[function]

include("atmosphere.jl")
export low_solar_activity, moderate_solar_activity, high_solar_activity, atmosphere_characteristics, AtmosphereEntry

include("magneticField.jl")
export magnetic_field_eci

include("sunVector.jl")
export sun_vector_eci

end