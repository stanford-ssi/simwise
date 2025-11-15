using Test

@testset "Dynamics Tests" begin
    include("energy.jl")
    include("orbit.jl")
    include("attitude.jl")
    include("rigidBody.jl")
end

nothing