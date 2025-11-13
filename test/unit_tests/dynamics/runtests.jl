using Test

@testset "Dynamics Tests" begin
    include("energy.jl")
    include("orbit.jl")
    include("attitude.jl")
end

nothing