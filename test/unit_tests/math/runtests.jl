using Test

@testset "Math Tests" begin
    include("rk4.jl")
    include("quaternions.jl")
    include("transforms.jl")
end
