using Test

include("../src/Simwise.jl")
using .Simwise

@testset "Simwise Tests" begin
    include("propagator.jl")
    include("orbit.jl")
    include("attitude.jl")
    include("dynamics.jl")
    include("environment.jl")
end
