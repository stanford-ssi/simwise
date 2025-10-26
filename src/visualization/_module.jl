module Visualization

using LinearAlgebra

# Call every non-exported name with Simwise.Visualization.[function]

include("attitudeViz.jl")
include("orbitViz.jl")
include("visualize.jl")

# Export the main visualization function
export visualize

end