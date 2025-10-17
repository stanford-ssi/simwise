using Test

using Simwise.Simulation: SensorPacketData, ReceivePacketData, ADCSDriver, send_data!, receive_data, wait_for_ready, reset_adcs!

@testset "Serial Interface - Packet Structures" begin
    # Test SensorPacketData construction
    @testset "SensorPacketData" begin
        sensor_data = SensorPacketData(
            Float32(60000.0),                                           # sim_mjd
            (Float32(0.1), Float32(0.2), Float32(0.3)),                 # w_body_raw
            (Float32(1.0e-5), Float32(2.0e-5), Float32(3.0e-5)),        # b_field_local
            ntuple(i -> UInt16(1000 + i), 16),                          # sun_sensors
            Float32(37.4), Float32(-122.1), Float32(400e3),             # gps_lat, lon, alt
            Float32(120000.0),                                          # gps_time
            (Float32(2025.0), Float32(10.0), Float32(13.0)),            # UTC_date
            Float32(5.0), Float32(12.0), Float32(0.417)                 # power, voltage, current
        )

        @test sensor_data.sim_mjd == Float32(60000.0)
        @test sensor_data.w_body_raw == (Float32(0.1), Float32(0.2), Float32(0.3))
        @test sensor_data.b_field_local[1] ≈ Float32(1.0e-5)
        @test length(sensor_data.sun_sensors) == 16
        @test sensor_data.sun_sensors[1] == UInt16(1001)
        @test sensor_data.gps_lat ≈ Float32(37.4)
        @test sensor_data.adcs_power ≈ Float32(5.0)
    end

    # Test ReceivePacketData construction
    @testset "ReceivePacketData" begin
        receive_data = ReceivePacketData(
            (Float32(1.0), Float32(0.0), Float32(0.0), Float32(0.0)),          # q_eci_to_principal
            (Float32(0.01), Float32(0.02), Float32(0.03)),              # w_principal
            (Float32(0.577), Float32(0.577), Float32(0.577)),           # sun_vector_principal
            (Float32(0.001), Float32(0.002), Float32(0.003)),           # magdrv_requested
            (Float32(100.0), Float32(100.0), Float32(100.0), Float32(100.0))   # w_reaction_wheels
        )

        @test receive_data.q_eci_to_principal == (Float32(1.0), Float32(0.0), Float32(0.0), Float32(0.0))
        @test receive_data.w_principal == (Float32(0.01), Float32(0.02), Float32(0.03))
        @test receive_data.magdrv_requested[1] ≈ Float32(0.001)
        @test length(receive_data.w_reaction_wheels) == 4
    end
end

# Note: The following tests require actual hardware connection to pass.
# They are included here for completeness but should be run in an environment
# where the ADCS hardware is connected.
@testset "Serial Interface - ADCS Driver" begin
    adcs = ADCSDriver("/dev/tty.usbmodem1101"; baudrate=115200, wait_for_port=true)
    # Test ADCSDriver initialization
    @testset "ADCSDriver Initialization" begin
        @test adcs.port_name == "/dev/tty.usbmodem1101"
        @test adcs.port !== nothing
    end

    # Test resetting ADCS 
    @testset "Resetting ADCS" begin
        @test reset_adcs!(adcs) == true
    end

    # Test sending data
    @testset "ADCSDriver Send Data" begin
        sensor_data = SensorPacketData(
            Float32(60000.0),
            (Float32(0.1), Float32(0.2), Float32(0.3)),
            (Float32(1.0e-5), Float32(2.0e-5), Float32(3.0e-5)),
            ntuple(i -> UInt16(1000 + i), 16),
            Float32(37.4), Float32(-122.1), Float32(400e3),
            Float32(120000.0),
            (Float32(2025.0), Float32(10.0), Float32(13.0)),
            Float32(5.0), Float32(12.0), Float32(0.417)
        )
        send_data!(adcs, sensor_data)
        @test true  # If we got here, send succeeded
    end

    # Test receiving data
    @testset "ADCSDriver Receive Data" begin
        received_packet = receive_data(adcs)
        @test received_packet !== nothing
        @test length(received_packet.q_eci_to_principal) == 4
        @test length(received_packet.w_principal) == 3
        @test length(received_packet.sun_vector_principal) == 3
        @test length(received_packet.magdrv_requested) == 3
        @test length(received_packet.w_reaction_wheels) == 4
    end
end
