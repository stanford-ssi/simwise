using Test

using Simwise.Math: Quat, quat_apply, quat_conj, quat_inv, quat_mult, to_vector, q_to_axis_angle, angle_axis_to_q, normalize, normalize!, dcm_to_q, q_to_dcm
using Simwise: RAD_TO_DEG, DEG_TO_RAD
@testset "Quaternion Tests" begin
    
    @testset "Struct and Initialization" begin
        
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
        @test inv == Quat(1 / length, -2 / length, -3 / length, -4 / length)
    end

    @testset "Normalization" begin
        q = Quat(1,0,0,0)
        @test normalize(q) == Quat(1,0,0,0)

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

    @testset "Axis Angle to Q" begin
        @testset "Identity" begin
            # Angle of 0 should not affect the quaternion
            @test angle_axis_to_q(0.0, [1,0,0]) == Quat(1,0,0,0)
            @test angle_axis_to_q(0.0, [1,1,0]) == Quat(1,0,0,0)
            @test angle_axis_to_q(0.0, [1,0,1]) == Quat(1,0,0,0)
        end

        @testset "180 Degrees" begin
            @test angle_axis_to_q(180.0, [1,0,0], true) == Quat(0,1,0,0) # Around X axis
            @test angle_axis_to_q(180.0, [0,1,0], true) == Quat(0,0,1,0) # Around Y axis
            @test angle_axis_to_q(180.0, [0,0,1], true) == Quat(0,0,0,1) # Around Z axis
        end

        @testset "90 Degrees" begin
            @test angle_axis_to_q(90.0, [1,0,0], true) == Quat(1/sqrt(2),1/sqrt(2),0,0) # Around X axis
            @test angle_axis_to_q(90.0, [0,1,0], true) == Quat(1/sqrt(2),0,1/sqrt(2),0) # Around Y axis
            @test angle_axis_to_q(90.0, [0,0,1], true) == Quat(1/sqrt(2),0,0,1/sqrt(2)) # Around Z axis
        end

        @testset "All Axes" begin
            @test angle_axis_to_q(90.0, [1,1,1], true) == Quat(1/sqrt(2), 1/sqrt(2*3), 1/sqrt(2*3), 1/sqrt(2*3))
            @test angle_axis_to_q(90.0, [-1,-1,-1], true) == Quat(1/sqrt(2), -1/sqrt(2*3), -1/sqrt(2*3), -1/sqrt(2*3))
        end

        @testset "Example in PDF" begin
            # https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
            @test angle_axis_to_q(2*pi/3, [1,1,1]) == Quat(0.5, 0.5, 0.5, 0.5)

        end
    end

    @testset "Q to Axis Angle" begin
        @testset "Identity" begin
            @test q_to_axis_angle(Quat(1,0,0,0), true) == (0, [0,0,0])
        end

        @testset "180 Degrees" begin
            @test q_to_axis_angle(Quat(0,1,0,0), true) == (180, [1,0,0])
            @test q_to_axis_angle(Quat(0,0,1,0), true) == (180, [0,1,0])
            @test q_to_axis_angle(Quat(0,0,0,1), true) == (180, [0,0,1])
        end

        @testset "90 Degrees" begin
            angle, axis = q_to_axis_angle(Quat(1/sqrt(2),1/sqrt(2),0,0), true)
            @test isapprox(angle, 90)
            @test isapprox(axis, [1,0,0])
            angle, axis = q_to_axis_angle(Quat(1/sqrt(2),0,1/sqrt(2),0), true)
            @test isapprox(angle, 90)
            @test isapprox(axis, [0,1,0])
            angle, axis = q_to_axis_angle(Quat(1/sqrt(2),0,0,1/sqrt(2)), true)
            @test isapprox(angle, 90)
            @test isapprox(axis, [0,0,1])
        end

        @testset "All Axes" begin
            angle, axis = q_to_axis_angle(Quat(1/sqrt(2), 1/sqrt(2*3), 1/sqrt(2*3), 1/sqrt(2*3)), true)
            @test isapprox(angle, 90.0)
            @test isapprox(axis, [1/sqrt(3), 1/sqrt(3), 1/sqrt(3)])

            angle, axis = q_to_axis_angle(Quat(1/sqrt(2), -1/sqrt(2*3), -1/sqrt(2*3), -1/sqrt(2*3)), true)
            @test isapprox(angle, 90.0)
            @test isapprox(axis, [-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)])
        end

        @testset "Example in PDF" begin
            # https://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
            angle, axis = q_to_axis_angle(Quat(0.5, 0.5, 0.5, 0.5))
            @test isapprox(angle, 2*pi/3)
            @test isapprox(axis, [1/sqrt(3),1/sqrt(3),1/sqrt(3)])
        end

    end

    @testset "DCM to Q" begin

        @testset "Identity" begin
            dcm = [
                1.0 0.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 1.0
                ]
            q = dcm_to_q(dcm)
                
            @test q == Quat(1,0,0,0)
        end

        # TODO: MAKE THIS WORK FOR NON-NORMALIZED MATRICES
        # @testset "Non Normalized" begin
        #     dcm = inv([
        #         0.0008382 -0.1452129  0.9894001;
        #         0.2022121  0.9689857  0.1420454;
        #         -0.9793414  0.1999496  0.0301760
        #     ]) * 10

        #     @test isapprox(to_vector(dcm_to_q(dcm)), to_vector(Quat(0.7071068, 0.0204722, 0.6960552, 0.1228333)), atol=1e-6)

        # end

        @testset "180 Degree Axes rotations" begin
            dcm_x = inv([
                1.0 0.0 0.0;
                0.0 -1.0 0.0;
                0.0 0.0 -1.0
                ])
            @test dcm_to_q(dcm_x) == Quat(0,1,0,0)

            dcm_y = inv([
                -1.0 0.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 -1.0
                ])
            @test dcm_to_q(dcm_y) == Quat(0,0,1,0)

            dcm_z = inv([
                -1.0 0.0 0.0;
                0.0 -1.0 0.0;
                0.0 0.0 1.0
                ])
            @test dcm_to_q(dcm_z) == Quat(0,0,0,1)
        end

        @testset "Random rotations" begin
            # https://www.andre-gaschler.com/rotationconverter/

            dcm = inv([
                0.9772839  -0.1380712  0.1607873;
                0.1607873   0.9772839 -0.1380712;
                -0.1380712   0.1607873  0.9772839
                ])
            @test isapprox(to_vector(dcm_to_q(dcm)), to_vector(Quat(0.9914449, 0.0753593, 0.0753593, 0.0753593)), atol=1e-6)

            dcm = inv([
                0.9205715 -0.0656093  0.3850241;
                0.0701405  0.9975345  0.0022810;
                -0.3842245  0.0249059  0.9229037
            ])
            @test isapprox(to_vector(dcm_to_q(dcm)), to_vector(Quat(0.9799247, 0.0057721, 0.196252, 0.0346327)), atol=1e-6)

            dcm = inv([
                0.0008382 -0.1452129  0.9894001;
                0.2022121  0.9689857  0.1420454;
                -0.9793414  0.1999496  0.0301760
            ])
            @test isapprox(to_vector(dcm_to_q(dcm)), to_vector(Quat(0.7071068, 0.0204722, 0.6960552, 0.1228333)), atol=1e-6)

        end
    end

    @testset "Q to DCM" begin

        @testset "Identity" begin
            q = Quat(1,0,0,0)
            dcm = [
                1.0 0.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 1.0
                ]
                
            @test q_to_dcm(q) == dcm
        end

        @testset "180 Degree Axes rotations" begin
            q_x = Quat(0,1,0,0)
            dcm_x = inv([
                1.0 0.0 0.0;
                0.0 -1.0 0.0;
                0.0 0.0 -1.0
                ])
            @test q_to_dcm(q_x) == dcm_x

            q_y = Quat(0,0,1,0)
            dcm_y = inv([
                -1.0 0.0 0.0;
                0.0 1.0 0.0;
                0.0 0.0 -1.0
                ])
            @test q_to_dcm(q_y) == dcm_y

            q_z = Quat(0,0,0,1)
            dcm_z = inv([
                -1.0 0.0 0.0;
                0.0 -1.0 0.0;
                0.0 0.0 1.0
                ])
            @test q_to_dcm(q_z) == dcm_z
        end

        @testset "Random rotations" begin
            # https://www.andre-gaschler.com/rotationconverter/

            q = Quat(0.9914449, 0.0753593, 0.0753593, 0.0753593)
            dcm = inv([
                0.9772839  -0.1380712  0.1607873;
                0.1607873   0.9772839 -0.1380712;
                -0.1380712   0.1607873  0.9772839
                ])
            @test isapprox(q_to_dcm(q), dcm, atol=1e-6)

            q = Quat(0.9799247, 0.0057721, 0.196252, 0.0346327)
            dcm = inv([
                0.9205715 -0.0656093  0.3850241;
                0.0701405  0.9975345  0.0022810;
                -0.3842245  0.0249059  0.9229037
            ])
            @test isapprox(q_to_dcm(q), dcm, atol=1e-6)

            q = Quat(0.7071068, 0.0204722, 0.6960552, 0.1228333)
            dcm = inv([
                0.0008382 -0.1452129  0.9894001;
                0.2022121  0.9689857  0.1420454;
                -0.9793414  0.1999496  0.0301760
            ])
            @test isapprox(q_to_dcm(q), dcm, atol=1e-6)
        end

        @testset "Non normalized" begin
            q = Quat(10.0 * 0.9799247, 10.0 * 0.0057721, 10.0 * 0.196252, 10.0 * 0.0346327)
            dcm = inv([
                0.9205715 -0.0656093  0.3850241;
                0.0701405  0.9975345  0.0022810;
                -0.3842245  0.0249059  0.9229037
            ])
            @test isapprox(q_to_dcm(q), dcm, atol=1e-6)
        end
    end

    @testset "DCM Round-trip Consistency" begin
        quats = [
            Quat(1, 0, 0, 0),
            Quat(0.9914449, 0.0753593, 0.0753593, 0.0753593),
            Quat(0.7071068, 0.0204722, 0.6960552, 0.1228333)
        ]
        for q in quats
            dcm = q_to_dcm(q)
            q2 = dcm_to_q(dcm)
            @test isapprox(to_vector(q2), to_vector(normalize(q)), atol=1e-6)
        end
    end

    @testset "DCM Quaternion Sign Equivalence" begin
        q = Quat(0.7071068, 0.0204722, 0.6960552, 0.1228333)
        dcm = q_to_dcm(q)
        @test q_to_dcm(-q) == dcm
    end

    @testset "Apply" begin

        @testset "Z axis 45 deg around X axis - Active" begin
            """
            Starts with quaternion representing a 90deg rotation around the x axis.
            The matrix expresses the vector in coordiante system A (`V_a`= `[0,0,1]`) 
            as the same vector from a rotated coordinate system B (`V_b` = `[0,1,0]`)
            """
            q = angle_axis_to_q(90.0, [1,0,0], true)
            v = Float64[0,0,1]
            v_result = quat_apply(q, v, false)
            @test isapprox(v_result, [0, -1, 0])
        end

        @testset "Z axis from frame rotated 90 deg around X axis - passive" begin
            """
            Starts with quaternion representing a 90deg rotation around the x axis.
            The matrix rotates the vector `V` `[0,0,1]` to produce `V'` `[0,1,0]`
            """
            q = angle_axis_to_q(90.0, [1,0,0], true)
            v = Float64[0,0,1]
            v_result = quat_apply(q, v)
            @test isapprox(v_result, [0, 1, 0])
        end

        @testset "Z axis 1e-6 deg around X axis - Active" begin
            """
            Starts with quaternion representing a 90deg rotation around the x axis.
            The matrix expresses the vector in coordiante system A (`V_a`= `[0,0,1]`) 
            as the same vector from a rotated coordinate system B (`V_b` = `[0,1,0]`)
            """
            q = angle_axis_to_q(1e-6, [1,0,0], true)
            v = Float64[0,0,1]
            v_result = quat_apply(q, v, false)
            @test isapprox(v_result, [0, -sin(1e-6 * DEG_TO_RAD), cos(1e-6 * DEG_TO_RAD)])
        end

        @testset "Z axis from frame rotated 1e-6 deg around X axis - passive" begin
            """
            Starts with quaternion representing a 90deg rotation around the x axis.
            The matrix rotates the vector `V` `[0,0,1]` to produce `V'` `[0,1,0]`
            """
            q = angle_axis_to_q(1e-6, [1,0,0], true)
            v = Float64[0,0,1]
            v_result = quat_apply(q, v)
            @test isapprox(v_result, [0, sin(1e-6 * DEG_TO_RAD), cos(1e-6 * DEG_TO_RAD)])
        end

        @testset "Weird rotations" begin
            """
            Starts with quaternion representing a 90deg rotation around the x axis.
            The matrix rotates the vector `V` `[0,0,1]` to produce `V'` `[0,1,0]`
            """
            q = angle_axis_to_q(90.0, [0,0,1], true)
            v = Float64[1,1,1]
            v_result = quat_apply(q, v)
            @test isapprox(v_result, [1, -1, 1])

            # From https://www.vcalc.com/wiki/vector-rotation
            q = angle_axis_to_q(13.0, [2,13,5.6], true)
            v = Float64[1,2,7]
            v_result = quat_apply(q, v, false)
            @test isapprox(v_result, [2.2469464036,1.9261221082,6.7261642474])
        end


    end

end

nothing