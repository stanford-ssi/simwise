## milestone #1 - dual attitude / orbit propagator, no controls

```
src/
├── Simwise.jl            # main module file
├── constants.jl
├── state.jl
├── parameters.jl
├── propagator/
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   ├── orbit.jl
│   └── torques/
│       ├── drag.jl
│       └── gravityGradient.jl
├── environment/
│   └── atmosphere.jl
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl
└── utils/
    └── transforms.jl

test/
├── runtests.jl
├── propagatorTests.jl
├── dynamicsTests.jl
└── environmentTests.jl

examples/
├── simpleOrbit.jl
└── tumblingCubesat.jl
```

## milestone #2 - compute sensor data

```
src/
├── Simwise.jl
├── constants.jl
├── state.jl
├── parameters.jl
├── serialInterface.jl            [NEW]
├── propagator/
│   └── rk4.jl
├── dynamics/
│   ├── attitude.jl
│   ├── orbit.jl
│   ├── forces/
│   │   └── j2.jl                 [NEW]
│   └── torques/
│       ├── drag.jl
│       ├── gravityGradient.jl
│       └── actuator.jl           [NEW]
├── environment/
│   ├── atmosphere.jl
│   ├── magneticField.jl          [NEW]
│   └── j2.jl                     [NEW]
├── sensors/                      [NEW]
│   ├── sunSensors.jl
│   ├── magnetometer.jl
│   ├── imu.jl
│   └── gps.jl
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl
└── utils/
    └── transforms.jl

test/
├── runtests.jl
├── propagatorTests.jl
├── dynamicsTests.jl
├── environmentTests.jl
└── sensorTests.jl                [NEW]

examples/
├── simpleOrbit.jl
└── tumblingCubesat.jl
```

## milestone #3 - hardware in the loop w/ actuator feedback

```
src/
├── Simwise.jl
├── constants.jl
├── state.jl
├── parameters.jl
├── serialInterface.jl
├── propagator/
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
├── environment/
│   ├── atmosphere.jl
│   ├── magneticField.jl
│   └── j2.jl
├── sensors/
│   ├── sunSensors.jl
│   ├── magnetometer.jl
│   ├── imu.jl
│   └── gps.jl
├── actuators/                    [NEW]
│   ├── reactionWheels.jl
│   └── magnetorquers.jl
├── visualization/
│   ├── orbitViz.jl
│   └── attitudeViz.jl
└── utils/
    └── transforms.jl

test/
├── runtests.jl
├── propagatorTests.jl
├── dynamicsTests.jl
├── environmentTests.jl
├── sensorTests.jl
└── actuatorTests.jl              [NEW]

examples/
├── simpleOrbit.jl
├── tumblingCubesat.jl
└── hilSimulation.jl              [NEW]
```