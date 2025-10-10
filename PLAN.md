## milestone #1 - dual attitude / orbit propagator with simulated sensor data

```
src/
├── Simwise.jl            # main module file
├── constants.jl
├── simulator.jl
├── satellite/
│   └── state.jl
│   └── parameters.jl
│   └── sensors/                    [TODO]
│       ├── sunSensors.jl
│       ├── magnetometer.jl
│       ├── imu.jl
│       └── gps.jl
├── math/
|   └── quaternions.jl
|   └── transforms.jl
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   └── orbit.jl
├── world/                          [TODO]
│   ├── magneticField.jl          
│   └── sunVector.jl              
├── visualization/                  [TODO]
│   ├── orbitViz.jl
│   └── attitudeViz.jl
```

## milestone #2 - send sensor data to flight computer

```
src/
├── Simwise.jl            # main module file
├── constants.jl
├── simulator.jl
├── serialInterface.jl        [TODO - and flight computer side]
├── satellite/
│   └── state.jl
│   └── parameters.jl
│   └── sensors/
│       ├── sunSensors.jl
│       ├── magnetometer.jl
│       ├── imu.jl
│       └── gps.jl
├── math/
|   └── quaternions.jl
|   └── transforms.jl
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   └── orbit.jl
├── world/                        
│   ├── magneticField.jl          
│   └── sunVector.jl              
├── visualization/                
│   ├── orbitViz.jl
│   └── attitudeViz.jl
```


## milestone #3 - perturbation torques + actuator inputs from flight computer HITL

```
src/
├── Simwise.jl
├── constants.jl
├── simulator.jl
├── satellite/
│   ├── state.jl
│   ├── parameters.jl
│   ├── serialInterface.jl          [TODO - add actuator receiving side]
│   └── sensors/                
│       ├── sunSensors.jl
│       ├── magnetometer.jl
│       ├── imu.jl
│       └── gps.jl
├── math/
|   └── quaternions.jl
|   └── transforms.jl
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   ├── orbit.jl
│   └── torques/                    [TODO]
│       ├── drag.jl
│       ├── gravityGradient.jl
│       └── actuator.jl           
├── world/
│   ├── atmosphere.jl               [TODO]
│   ├── magneticField.jl          
│   └── sunVector.jl              
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl
```