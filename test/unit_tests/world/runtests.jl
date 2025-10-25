using Test

@testset "World Tests" begin
    include("atmosphere.jl")
    include("magneticField.jl")
    include("sunVector.jl")
end

nothing