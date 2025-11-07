using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: rigid_body_derivative
using Simwise.Math: rk4_step
#propagate_keplerian, state_from_oe


@testset "Rigid Body Tests" begin
    @testset "Force No Torque" begin

        # Initialized with (mass, inertia, force, torque) constructor
        params = Parameters(
            2.0, # kg
            diagm([2.0, 2.0, 2.0]), # Inertia
            [1.0,0.0,0.0], # x direction force
            zeros(3) # no torque
        )
        
        state = [
            0.0, 0.0, 0.0, # r = 0
            1.1, 1.2, 1.3, # v = 0
            1.0, 0.0, 0.0, 0.0, # Identity quaternion
            0.0, 0.0, 0.0 # ω = 0
            ]

        dot = rigid_body_derivative(0.0, state, params)

        @test dot[1:3] == [1.1, 1.2, 1.3]
        @test dot[4:6] == [0.5, 0.0, 0.0]
        @test dot[7:10] == [0.0, 0.0, 0.0, 0.0]
        @test dot[11:13] == [0.0, 0.0, 0.0]
    end

    @testset "Torque No Force" begin

        # Initialized with (mass, inertia, force, torque) constructor
        params = Parameters(
            2.0, # kg
            diagm([2.0, 2.0, 2.0]), # Inertia
            zeros(3), # no force
            [1.0, 2.0, 4.0] # x direction torque
        )
        
        state = [
            0.0, 0.0, 0.0, # r = 0
            1.1, 1.2, 1.3, # v = 0
            1.0, 0.0, 0.0, 0.0, # Identity quaternion
            0.0, 0.0, 0.0 # ω = 0
            ]

        dot = rigid_body_derivative(0.0, state, params)

        @test dot[1:3] == [1.1, 1.2, 1.3]
        @test dot[4:6] == [0.0, 0.0, 0.0]
        @test dot[7:10] == [0.0, 0.0, 0.0, 0.0]
        @test dot[11:13] == [0.5, 1.0, 2.0]
    end

    @testset "Complete random nonsense testing" begin

        I = [
            15.0 0.5 0.5;
            0.1 15.0 0.0;
            0.0 3.0 4.0
        ]

        # Initialized with (mass, inertia, force, torque) constructor
        params = Parameters(
            17.0, # kg
            I,
            [12.0, 200.0, 14.0], # force
            [0.4, 0.2, 0.7] # torque
        )
        
        state = [
            0.0, 0.0, 0.0, # r = 0
            1.1, 1.2, 1.3, # v = 0
            0.5, 0.5, 0.5, 0.5, # Identity quaternion
            2.0, 3.14, 3.52 # ω = 0
            ]

        dot = rigid_body_derivative(0.0, state, params)

        # Generated using Julia as a simple calculator with the textbook functions
        @test isapprox(dot[1:3], [1.1, 1.2, 1.3], atol=1e-12)
        @test isapprox(dot[4:6], [0.7058823529411765, 11.764705882352942, 0.8235294117647058], atol=1e-12)
        @test isapprox(dot[7:10], [-2.165, 0.595, 0.40500000000000014, 1.165], atol=1e-12)
        @test isapprox(dot[11:13], [6.1567301516750925, -4.715818201011166, 6.225913650758373], atol=1e-12)
    end

end

@testset "Rigid Body Integration" begin
    @testset "0 Time" begin
        
        dt = 0.0
        t = 1.0

        I = [
            15.0 0.5 0.5;
            0.1 15.0 0.0;
            0.0 3.0 4.0
        ]

        # Initialized with (mass, inertia, force, torque) constructor
        params = Parameters(
            17.0, # kg
            I,
            [12.0, 200.0, 14.0], # force
            [0.4, 0.2, 0.7] # torque
        )

        state = [
            0.0, 0.0, 0.0, # r = 0
            1.1, 1.2, 1.3, # v = 0
            0.5, 0.5, 0.5, 0.5, # Identity quaternion
            2.0, 3.14, 3.52 # ω = 0
            ]
        
        next_step = rk4_step(dt, t, state, rigid_body_derivative, params)

        @test next_step == state # No change in state
    end

end

nothing