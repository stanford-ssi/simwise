module World

using ..Simwise  # parent module, for shared constants etc.
using LinearAlgebra

# Call every non-exported name with Simwise.World.[function]

include("atmosphere.jl")

include("magneticField.jl")
export magnetic_field_eci

include("sunVector.jl")
export sun_vector_eci

end