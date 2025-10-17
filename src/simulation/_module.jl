module Simulation

using LinearAlgebra
using LibSerialPort

# Call every non-exported name with Simwise.Simulation.[function]

include("serialInterface.jl")
include("simulator.jl")

# Serial interface exports
export propagate
export SensorPacketData, ReceivePacketData, ADCSDriver, send_data!, receive_data, wait_for_ready, reset_adcs!

end