using Test

include("../src/Simwise.jl")
using .Simwise

@testset "Quaternion Tests" begin
    
    @testset "Test Struct" begin
        
        @testset "Default" begin
            q = Quat()

            @test q == Quat(1,0,0,0)
        end

        @testset "From Vec" begin
            v = Float64[1, 2, 3, 4]
            q = Quat(v)

            @test q == Quat(1,2,3,4)
        end

        @testset "To Vec" begin
            q = Quat(1, 2, 3, 4)
            v = to_vector(q)

            @test v == [1,2,3,4]
        end
    end

    @testset "Conjugate" begin
        result = quat_conj(Quat(1,0,0,0))
        @test result == Quat(1,0,0,0)

        conj = quat_conj(Quat(1, 2, 3, 4))
        @test conj == Quat(1, -2, -3, -4)
    end

    @testset "Inverse" begin
        inv = quat_inv(Quat(1, 2, 3, 4))
        length = 1^2 + 2^2 + 3^2 + 4^2
        @test inv == Quat(1 / length, 2 / length, 3 / length, 4 / length)
    end

    @testset "Normalization" begin
        q = Quat(2, 2, 2, 2)

        conj = normalize(q)
        @test conj == Quat(0.5, 0.5, 0.5, 0.5)

        normalize!(q)

        @test q == Quat(0.5, 0.5, 0.5, 0.5)
    end

    @testset "Multiply" begin
        q1 = Quat(2, 2, 2, 2)
        q2 = quat_mult(q1, q1)
        @test q2 == Quat(-8,8,8,8)

        q1 = Quat(3, 1, -2, 1)
        q2 = Quat(2, -1, 2, 3)
        q3 = quat_mult(q1, q2)
        @test q3 == Quat(8,-9,-2,11)
    end


    @testset "Apply" begin
        
    end

end
