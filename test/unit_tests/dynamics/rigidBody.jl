using Plots
using LinearAlgebra

using Simwise.Constants: μ_earth
using Simwise.Dynamics: rigid_body_derivative
#propagate_keplerian, state_from_oe


@testset "Rigid Body Tests" begin
    @testset "Force No Torque" begin

        q = Quat() # Identity
        ω = [0,0,pi] # Rotating by π rad/s

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
end

nothing