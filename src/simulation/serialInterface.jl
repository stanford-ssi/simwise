# Interface between ADCS flight software and physics-based simulator
# via serial connection. Designed for Hardware-in-the-Loop (HITL) testing.
"""
ADCS Simulation Interface Driver

This module provides a clean driver interface for communicating with the ADCS
flight software via serial. Designed to be integrated into physics-based simulations.

# Core Driver API:
- `ADCSDriver`: Stateful driver object managing serial connection
- `send_data!(driver, sensor_data)`: Send sensor data to ADCS
- `receive_data(driver)`: Receive actuator commands from ADCS
- `wait_for_ready(driver)`: Wait for ADCS initialization
- `reset_adcs!(driver)`: Reset the ADCS flight computer to initial state

# Example usage
driver = ADCSDriver("/dev/tty.usbmodem1101")
reset_adcs!(driver)
sensor_data = SensorPacketData(...)
send_data!(driver, sensor_data)
actuator_data = receive_data(driver)
```
"""

using LibSerialPort

# ============================================================================
# PACKET DEFINITIONS
# ============================================================================

const SENSOR_PACKET_SIZE = 108
const ACTUATOR_PACKET_SIZE = 76

"""
Sensor data sent from simulator to ADCS flight software.
Contains all sensor readings needed for GNC algorithms.
"""
struct SensorPacketData
    sim_mjd::Float32                    # Modified Julian Date
    w_body_raw::NTuple{3, Float32}      # Angular velocity [rad/s] in body frame
    b_field_local::NTuple{3, Float32}   # Magnetic field [T] in body frame
    sun_sensors::NTuple{16, UInt16}     # 16 sun sensor ADC readings
    gps_lat::Float32                    # GPS latitude [deg]
    gps_lon::Float32                    # GPS longitude [deg]
    gps_alt::Float32                    # GPS altitude [m]
    gps_time::Float32                   # GPS time [s]
    UTC_date::NTuple{3, Float32}        # (year, month, day)
    adcs_power::Float32                 # Power consumption [W]
    adcs_voltage::Float32               # Bus voltage [V]
    adcs_current::Float32               # Bus current [A]
end

"""
Actuator commands received from ADCS flight software.
Contains GNC state estimates and control outputs.
"""
struct ReceivePacketData
    q_eci_to_principal::NTuple{4, Float32}    # Attitude quaternion
    w_principal::NTuple{3, Float32}           # Angular velocity estimate [rad/s]
    sun_vector_principal::NTuple{3, Float32}  # Sun vector in body frame
    magdrv_requested::NTuple{3, Float32}      # Magnetorquer dipole [A�m�]
    w_reaction_wheels::NTuple{4, Float32}     # Reaction wheel speeds [rad/s]
end

# ============================================================================
# DRIVER OBJECT
# ============================================================================

"""
ADCS serial driver maintaining connection state and internal buffers.
"""
mutable struct ADCSDriver
    port::LibSerialPort.SerialPort
    port_name::String
    tx_buffer::Vector{UInt8}      # Pre-allocated transmission buffer
    rx_buffer::Vector{UInt8}      # Accumulation buffer for receiving
    header_search_pos::Int        # Current search position in rx_buffer

    function ADCSDriver(port_name::String; baudrate::Int=115200, wait_for_port::Bool=true)
        # Wait for port to become available (device may be rebooting)
        if wait_for_port
            println("� Waiting for serial port to become available...")
            max_attempts = 100
            port_ready = false
            for attempt in 1:max_attempts
                try
                    # Try to open and immediately close to test availability
                    test_port = LibSerialPort.open(port_name, baudrate)
                    close(test_port)
                    port_ready = true
                    println(" Serial port available!")
                    break
                catch e
                    if attempt == max_attempts
                        error("L Serial port not available after $(max_attempts) attempts: $e")
                    end
                    sleep(0.1)
                end
            end

            # Wait additional time for device to stabilize
            if port_ready
                println("�  Waiting 2 seconds for device to stabilize...")
                sleep(2.0)
            end
        end

        # Now open the port for real
        port = LibSerialPort.open(port_name, baudrate)

        tx_buffer = Vector{UInt8}(undef, max(SENSOR_PACKET_SIZE, ACTUATOR_PACKET_SIZE))
        rx_buffer = Vector{UInt8}()
        sizehint!(rx_buffer, ACTUATOR_PACKET_SIZE * 2)

        new(port, port_name, tx_buffer, rx_buffer, 1)
    end
end

Base.close(driver::ADCSDriver) = close(driver.port)

# ============================================================================
# LOW-LEVEL SERIALIZATION/DESERIALIZATION
# ============================================================================

"""Serialize sensor packet into buffer. Returns checksum."""
@inline function _serialize_sensor_packet!(buffer::Vector{UInt8}, pkt::SensorPacketData)
    idx = 1

    # Header + padding
    @inbounds buffer[idx] = UInt8('S'); idx += 1
    @inbounds buffer[idx] = UInt8('S'); idx += 1
    @inbounds buffer[idx] = 0x00; idx += 1
    @inbounds buffer[idx] = 0x00; idx += 1

    # Write Float32 values directly to buffer
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.sim_mjd); idx += 4

    @inbounds for v in pkt.w_body_raw
        unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), v)
        idx += 4
    end

    @inbounds for v in pkt.b_field_local
        unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), v)
        idx += 4
    end

    @inbounds for v in pkt.sun_sensors
        unsafe_store!(Ptr{UInt16}(pointer(buffer, idx)), v)
        idx += 2
    end

    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.gps_lat); idx += 4
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.gps_lon); idx += 4
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.gps_alt); idx += 4
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.gps_time); idx += 4

    @inbounds for v in pkt.UTC_date
        unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), v)
        idx += 4
    end

    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.adcs_power); idx += 4
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.adcs_voltage); idx += 4
    @inbounds unsafe_store!(Ptr{Float32}(pointer(buffer, idx)), pkt.adcs_current); idx += 4

    # Calculate checksum (idx-1 is last data byte)
    checksum = UInt16(0)
    @inbounds @simd for i in 1:(idx-1)
        checksum += UInt16(buffer[i])
    end

    # Write checksum and padding
    @inbounds buffer[idx] = UInt8(checksum & 0xFF); idx += 1
    @inbounds buffer[idx] = UInt8((checksum >> 8) & 0xFF); idx += 1
    @inbounds buffer[idx] = 0x00; idx += 1
    @inbounds buffer[idx] = 0x00

    return checksum
end

"""Parse and validate actuator packet. Returns ReceivePacketData or nothing."""
@inline function _parse_actuator_packet(data::Vector{UInt8})
    @inbounds length(data) != 76 && (return nothing)
    @inbounds (data[1] != UInt8('A') || data[2] != UInt8('A')) && (return nothing)

    # Fast checksum validation
    expected = UInt16(0)
    @inbounds @simd for i in 1:72
        expected += UInt16(data[i])
    end
    @inbounds received = UInt16(data[73]) | (UInt16(data[74]) << 8)
    expected != received && (return nothing)

    # Direct memory reads for much faster parsing
    idx = 5
    @inbounds q = ntuple(i -> (v = unsafe_load(Ptr{Float32}(pointer(data, idx))); idx += 4; v), 4)
    @inbounds w = ntuple(i -> (v = unsafe_load(Ptr{Float32}(pointer(data, idx))); idx += 4; v), 3)
    @inbounds sun = ntuple(i -> (v = unsafe_load(Ptr{Float32}(pointer(data, idx))); idx += 4; v), 3)
    @inbounds mag = ntuple(i -> (v = unsafe_load(Ptr{Float32}(pointer(data, idx))); idx += 4; v), 3)
    @inbounds wheels = ntuple(i -> (v = unsafe_load(Ptr{Float32}(pointer(data, idx))); idx += 4; v), 4)

    return ReceivePacketData(q, w, sun, mag, wheels)
end

# ============================================================================
# PUBLIC DRIVER API
# ============================================================================

"""
Wait for ADCS to complete initialization and send ready marker.
Blocks until "ADCS_READY" is received or timeout (default 10s).
"""
function wait_for_ready(driver::ADCSDriver; timeout_s::Float64=10.0)
    text_buffer = UInt8[]
    sizehint!(text_buffer, 256)
    ready_marker = b"ADCS_READY"
    start_time = time()

    while (time() - start_time) < timeout_s
        try
            byte = read(driver.port, 1)[1]
            if byte >= 0x20 && byte <= 0x7E || byte == UInt8('\n')  # Printable or newline
                push!(text_buffer, byte)

                # Check for ready marker - search the entire buffer
                if length(text_buffer) >= length(ready_marker)
                    for i in 1:(length(text_buffer)-length(ready_marker)+1)
                        if view(text_buffer, i:(i+length(ready_marker)-1)) == ready_marker
                            return true
                        end
                    end
                end

                # Prevent unbounded growth
                if length(text_buffer) > 200
                    text_buffer = text_buffer[end-100:end]
                end
            end
        catch e
            sleep(0.001)
        end
    end

    return false  # Timeout
end

"""
Reset the ADCS flight computer by sending a packet of all zeros.
This triggers the flight computer to return to its initial state.
After sending the reset packet, waits for the ADCS to reboot and send ADCS_READY.
"""
function reset_adcs!(driver::ADCSDriver; timeout_s::Float64=20.0)
    println("🔄 Resetting ADCS flight computer...")

    # Clear driver's receive buffer to avoid stale data
    empty!(driver.rx_buffer)
    driver.header_search_pos = 1

    # Send a packet of all zeros (SENSOR_PACKET_SIZE bytes)
    fill!(driver.tx_buffer, 0x00)
    packet_view = view(driver.tx_buffer, 1:SENSOR_PACKET_SIZE)
    write(driver.port, packet_view)

    println("   Reset packet sent, waiting for reboot...")

    # Wait for ADCS to reboot and send ready signal
    if wait_for_ready(driver, timeout_s=timeout_s)
        println("✅ ADCS reset complete and ready!\n")
        return true
    else
        println("⚠️  Warning: ADCS_READY not received after reset\n")
        return false
    end
end

"""
Send sensor data to ADCS flight software.
This is non-blocking - data is written to serial and returns immediately.
"""
function send_data!(driver::ADCSDriver, sensor_data::SensorPacketData)
    # Serialize into driver's tx buffer
    _serialize_sensor_packet!(driver.tx_buffer, sensor_data)

    # Write packet to serial
    packet_view = view(driver.tx_buffer, 1:SENSOR_PACKET_SIZE)
    write(driver.port, packet_view)

    return nothing
end

"""
Receive actuator data from ADCS flight software.
Blocks until a valid packet is received. Returns ReceivePacketData.

This function maintains internal state in driver.rx_buffer for streaming reception.
"""
function receive_data(driver::ADCSDriver)
    while true
        # Read available data
        bytes_avail = bytesavailable(driver.port)
        if bytes_avail > 0
            chunk_size = min(bytes_avail, 256)
            chunk = read(driver.port, chunk_size)
            append!(driver.rx_buffer, chunk)

            # Search for 'AA' header
            @inbounds for i in driver.header_search_pos:(length(driver.rx_buffer)-1)
                if driver.rx_buffer[i] == UInt8('A') && driver.rx_buffer[i+1] == UInt8('A')
                    # Found potential packet header
                    packet_end = i + ACTUATOR_PACKET_SIZE - 1
                    if packet_end <= length(driver.rx_buffer)
                        # Have full packet, try to parse
                        packet_view = view(driver.rx_buffer, i:packet_end)
                        actuator_pkt = _parse_actuator_packet(collect(packet_view))

                        if !isnothing(actuator_pkt)
                            # Valid packet! Remove processed data and return
                            deleteat!(driver.rx_buffer, 1:packet_end)
                            driver.header_search_pos = 1
                            return actuator_pkt
                        end
                    end
                end
            end

            # Update search position
            driver.header_search_pos = max(1, length(driver.rx_buffer) - 1)

            # Prevent buffer overflow
            if length(driver.rx_buffer) > 1024
                deleteat!(driver.rx_buffer, 1:(length(driver.rx_buffer)-512))
                driver.header_search_pos = 1
            end
        else
            sleep(0.0001)  # No data available, short sleep
        end
    end
end
