# Physical constants

using SatelliteDynamics

# Earth
const μ_earth = 3.986004418e14  # m^3/s^2 - Standard gravitational parameter
const R_earth = 6378137.0       # m - Earth equatorial radius
const ω_earth = 7.2921159e-5    # rad/s - Earth rotation rate

# Universal
const G = 6.67430e-11           # m^3/(kg·s^2) - Gravitational constant

# Time
const seconds_per_day = 86400.0
const mjd_2000 = SatelliteDynamics.MJD2000
const mjd_0 = SatelliteDynamics.MJD_ZERO

# Sun Sensor Constants
const NUM_SUN_SENSORS = 16  # Total number of sun sensors
const SQRT_2_INV = 1 / sqrt(2.0)
