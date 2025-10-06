# src/hardware/serial_interface.jl

using LibSerialPort
using StaticArrays

"""
    SensorData

Sensor readings to send to flight computer.
"""
struct SensorData
    timestamp::Float32          # seconds
    gyro::SVector{3,Float32}    # rad/s (x, y, z)
    accel::SVector{3,Float32}   # m/s² (x, y, z)
    mag::SVector{3,Float32}     # μT (x, y, z)
    photodiodes::SVector{16,Float32}  # raw ADC-like values
    gps_lat::Float32            # degrees
    gps_lon::Float32            # degrees
    gps_alt::Float32            # meters
end

"""
    ActuatorCommand

Actuator commands received from flight computer.
"""
struct ActuatorCommand
    timestamp::Float32              # seconds
    magnetorquer::SVector{3,Float32}  # -1 to 1 (x, y, z)
    reaction_wheels::SVector{3,Float32}  # speeds (TBD units)
end

"""
    SerialConnection

Manages async serial communication with flight computer.
"""
mutable struct SerialConnection
    port::SerialPort
    port_name::String
    baudrate::Int
    is_connected::Bool
    read_buffer::String
    last_command::Union{ActuatorCommand, Nothing}
    command_lock::ReentrantLock
end

"""
    SerialConnection(port_name::String; baudrate::Int=115200)

Open serial connection to flight computer.
"""
function SerialConnection(port_name::String; baudrate::Int=115200)
    try
        port = LibSerialPort.open(port_name, baudrate)
        set_speed(port, baudrate)
        set_frame(port, ndatabits=8, parity=SP_PARITY_NONE, nstopbits=1)
        set_flow_control(port, flowcontrol=SP_FLOWCONTROL_NONE)
        
        conn = SerialConnection(
            port,
            port_name,
            baudrate,
            true,
            "",
            nothing,
            ReentrantLock()
        )
        
        println("✓ Connected to flight computer on $port_name @ $baudrate baud")
        return conn
    catch e
        error("Failed to open serial port $port_name: $e")
    end
end

"""
    close(conn::SerialConnection)

Close serial connection.
"""
function Base.close(conn::SerialConnection)
    if conn.is_connected
        close(conn.port)
        conn.is_connected = false
        println("✓ Serial connection closed")
    end
end

"""
    encode_sensor_data(data::SensorData) -> String

Encode sensor data as ASCII message.
Format: S,timestamp,gx,gy,gz,ax,ay,az,mx,my,mz,pd0,...,pd15,lat,lon,alt\\n
"""
function encode_sensor_data(data::SensorData)
    # Build CSV string
    parts = String[
        "S",
        string(data.timestamp),
        string(data.gyro[1]), string(data.gyro[2]), string(data.gyro[3]),
        string(data.accel[1]), string(data.accel[2]), string(data.accel[3]),
        string(data.mag[1]), string(data.mag[2]), string(data.mag[3])
    ]
    
    # Add all 16 photodiode values
    for i in 1:16
        push!(parts, string(data.photodiodes[i]))
    end
    
    # Add GPS
    push!(parts, string(data.gps_lat))
    push!(parts, string(data.gps_lon))
    push!(parts, string(data.gps_alt))
    
    return join(parts, ",") * "\n"
end

"""
    decode_actuator_command(line::String) -> Union{ActuatorCommand, Nothing}

Decode actuator command from ASCII message.
Format: A,timestamp,mtq_x,mtq_y,mtq_z,rw1,rw2,rw3\\n
"""
function decode_actuator_command(line::String)
    parts = split(strip(line), ',')
    
    # Validate format
    if length(parts) != 8 || parts[1] != "A"
        @warn "Invalid actuator command format: $line"
        return nothing
    end
    
    try
        timestamp = parse(Float32, parts[2])
        mtq = SVector{3,Float32}(
            parse(Float32, parts[3]),
            parse(Float32, parts[4]),
            parse(Float32, parts[5])
        )
        rw = SVector{3,Float32}(
            parse(Float32, parts[6]),
            parse(Float32, parts[7]),
            parse(Float32, parts[8])
        )
        
        # Validate magnetorquer range
        if any(abs.(mtq) .> 1.0)
            @warn "Magnetorquer values out of range [-1, 1]: $mtq"
        end
        
        return ActuatorCommand(timestamp, mtq, rw)
    catch e
        @warn "Failed to parse actuator command: $e"
        return nothing
    end
end

"""
    send_sensor_data(conn::SerialConnection, data::SensorData)

Send sensor data to flight computer (non-blocking).
"""
function send_sensor_data(conn::SerialConnection, data::SensorData)
    if !conn.is_connected
        @warn "Serial connection not open"
        return false
    end
    
    try
        message = encode_sensor_data(data)
        write(conn.port, message)
        return true
    catch e
        @warn "Failed to send sensor data: $e"
        return false
    end
end

"""
    read_actuator_command(conn::SerialConnection) -> Union{ActuatorCommand, Nothing}

Non-blocking read of actuator command from flight computer.
Returns latest command if available, otherwise returns nothing.
"""
function read_actuator_command(conn::SerialConnection)
    if !conn.is_connected
        return nothing
    end
    
    try
        # Non-blocking read
        nbytes = bytesavailable(conn.port)
        if nbytes > 0
            incoming = String(read(conn.port, nbytes))
            conn.read_buffer *= incoming
            
            # Look for complete lines (ending with \n)
            while contains(conn.read_buffer, '\n')
                line_end = findfirst('\n', conn.read_buffer)
                line = conn.read_buffer[1:line_end-1]
                conn.read_buffer = conn.read_buffer[line_end+1:end]
                
                # Parse the line
                cmd = decode_actuator_command(line)
                if !isnothing(cmd)
                    lock(conn.command_lock) do
                        conn.last_command = cmd
                    end
                end
            end
        end
        
        # Return latest command
        lock(conn.command_lock) do
            return conn.last_command
        end
    catch e
        @warn "Error reading from serial: $e"
        return nothing
    end
end

"""
    get_latest_command(conn::SerialConnection) -> Union{ActuatorCommand, Nothing}

Get the most recent actuator command (thread-safe).
"""
function get_latest_command(conn::SerialConnection)
    lock(conn.command_lock) do
        return conn.last_command
    end
end

"""
    clear_command(conn::SerialConnection)

Clear the stored command (useful after processing).
"""
function clear_command!(conn::SerialConnection)
    lock(conn.command_lock) do
        conn.last_command = nothing
    end
end

# Convenience function for list_ports
"""
    list_serial_ports()

List available serial ports.
"""
function list_serial_ports()
    ports = LibSerialPort.get_port_list()
    if isempty(ports)
        println("No serial ports found")
    else
        println("Available serial ports:")
        for port in ports
            println("  - $port")
        end
    end
    return ports
end