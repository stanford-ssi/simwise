using Test

@testset "Dynamics Tests" begin
    include("orbit.jl")
    include("attitude.jl")
end

nothing