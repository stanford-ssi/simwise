# Simwise.jl

6-DOF orbital and attitude propagator for SAMWISE ADCS hardware-in-the-loop (HITL) testing.

## Features

- **6-DOF Dynamics**: Orbital mechanics and attitude kinematics with environmental models (gravity, magnetic field, drag)
- **Sensor Models**: Gyroscope, magnetometer, sun sensor, GPS
- **Actuator Models**: Reaction wheels, magnetorquers, thrusters
- **Serial Interface**: Real-time communication with flight computers
- **Flexible Integrators**: RK4, RK45, and more via DifferentialEquations.jl

## Quick Start

```julia
using Simwise

# Initialize satellite state
state = SatelliteState(
    position = [500.0, 0.0, 0.0],      # km (ECI)
    velocity = [0.0, 7.5, 0.0],         # km/s
    quaternion = [0.0, 0.0, 0.0, 1.0],  # body to ECI
    angular_velocity = [0.0, 0.0, 0.0]  # rad/s
)
# TODO: check spacex-provided parameters

# Propagate 60 seconds
new_state = propagate!(state, 60.0, config)

# Hardware-in-the-loop
serial = SerialConnection("/dev/ttyUSB0")
run_hil_simulation(duration=3600.0, dt=0.01)
```

## Structure

```
src/
├── dynamics/      # Orbital & attitude dynamics
├── propagators/   # Numerical integrators
├── sensors/       # Sensor models
├── actuators/     # Actuator models
├── hardware/      # Serial communication & protocols
└── utils/         # Frames, quaternions, constants
```

### Key Resources
- [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl/tree/master?tab=readme-ov-file) for core math and modeling
    - [SatelliteToolboxGeomagneticField.jl](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl) for IGRF 14 modeling
    - [SatelliteToolboxAtmosphericModels.jl](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl) for atmospheric models
    - etc.
- [SatelliteDynamics.jl](https://sisl.github.io/SatelliteDynamics.jl/latest/) for attitude and orbital dynamics, integrators, and propagators
    - Made by postdoc Duncan Eddy @ Stanford