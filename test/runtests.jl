using Test

@testset "Simwise Tests" begin
    include("dynamics/runtests.jl")
    include("math/runtests.jl")
    include("world/runtests.jl")
end
