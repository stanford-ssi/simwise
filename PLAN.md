## milestone #1 - dual attitude / orbit propagator, no controls

```
src/
├── Simwise.jl            # main module file
├── constants.jl
├── simulator.jl
├── satellite/
│   └── state.jl
│   └── parameters.jl
├── math/
|   └── quaternions.jl
|   └── transforms.jl
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   ├── orbit.jl
│   └── torques/
│       ├── drag.jl
│       └── gravityGradient.jl
├── world/
│   └── atmosphere.jl
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl


test/
├── runtests.jl
├── attitude.jl
├── orbit.jl
├── rk4.jl
└── plots/

examples/
├── simpleOrbit.jl
└── tumblingCubesat.jl
```

## milestone #2 - compute sensor data

```
src/
├── Simwise.jl
├── constants.jl
├── simulator.jl
├── satellite/
│   ├── state.jl
│   ├── parameters.jl
│   ├── serialInterface.jl        [NEW]
│   └── sensors/                  [NEW]
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
│   ├── forces/                   [NEW]
│   │   └── j2.jl                 [NEW]
│   └── torques/
│       ├── drag.jl
│       ├── gravityGradient.jl
│       └── actuator.jl           [NEW]
├── world/
│   ├── atmosphere.jl
│   ├── magneticField.jl          [NEW]
│   └── j2.jl                     [NEW]
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl

test/
├── runtests.jl
├── attitude.jl
├── orbit.jl
├── rk4.jl
├── sensorTests.jl                [NEW]
└── plots/

examples/
├── simpleOrbit.jl
├── tumblingCubesat.jl
└── hilSimulation.jl              

```

## milestone #3 - hardware in the loop w/ actuator feedback

```
src/
├── Simwise.jl
├── constants.jl
├── simulator.jl
├── satellite/
│   ├── state.jl
│   ├── parameters.jl
│   ├── serialInterface.jl
│   ├── sensors/
│   │   ├── sunSensors.jl
│   │   ├── magnetometer.jl
│   │   ├── imu.jl
│   │   └── gps.jl
│   └── actuators/                [NEW]
│       ├── reactionWheels.jl
│       └── magnetorquers.jl
├── math/
|   └── quaternions.jl
|   └── transforms.jl
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   ├── orbit.jl
│   ├── forces/
│   │   └── j2.jl
│   └── torques/
│       ├── drag.jl
│       ├── gravityGradient.jl
│       └── actuator.jl
├── world/
│   ├── atmosphere.jl
│   ├── magneticField.jl
│   └── j2.jl
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl

test/
├── runtests.jl
├── attitude.jl
├── orbit.jl
├── rk4.jl
├── sensorTests.jl
├── actuatorTests.jl              [NEW]
└── plots/

examples/
├── simpleOrbit.jl
├── tumblingCubesat.jl
└── hilSimulation.jl              
```