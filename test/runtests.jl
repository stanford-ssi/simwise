using Test
using Simwise

@testset "Simwise Tests" begin
    include("propagatorTests.jl")
    include("dynamicsTests.jl")
    include("environmentTests.jl")
end
