using Test

include("../src/Simwise.jl")
using .Simwise

@testset "Simwise Tests" begin
    include("rk4.jl")
    include("orbit.jl")
    include("attitude.jl")
    include("magneticField.jl")
    include("sunVector.jl")
end
