using SatelliteDynamics
using Dates
using DataStructures

# Atmospheric density model
using ..Constants: R_earth

# Define a simple struct
struct AtmosphereEntry
    alt_km::Float64
    temp::Float64
    density::Float64
    pressure::Float64
    mol_wt::Float64
end

"""
    atmosphere_characteristics(altitude_km)

http://www.braeunig.us/space/atmos.htm

# Arguments
- `altitude_km::Float64`: Altitude above Earth surface [km]
- `solar_activity::String`: "high", "moderate", "low" (Defaults to "moderate")

# Returns
- `density::Float64`: Atmospheric density [kg/m3]
- `temp::Float64`: Temperature [K]
- `pressure::Float64`: Pressure [Pa]
- `mol_wt::Float64`: Molecular Weight of air [kg/kmol]

"""
function atmosphere_characteristics(altitude_km::Float64, solar_activity::String = "moderate")

    activity_table =
        solar_activity == "low"   ? low_solar_activity :
        solar_activity == "high"  ? high_solar_activity :
                                    moderate_solar_activity

    # Extract sorted keys and values
    alts = collect(keys(activity_table))
    entries = collect(values(activity_table))

    # Find indices bounding the altitude
    upper_idx = searchsortedfirst(alts, altitude_km)
    lower_idx = clamp(upper_idx - 1, 1, length(alts))
    upper_idx = clamp(upper_idx, 1, length(alts))

    lower_entry = entries[lower_idx]
    upper_entry = entries[upper_idx]

    # Clamp
    if altitude_km <= alts[1]
        return AtmosphereEntry(alts[1], entries[1].temp, entries[1].density, entries[1].pressure, entries[1].mol_wt)
    elseif altitude_km >= alts[end]
        return AtmosphereEntry(alts[end], entries[end].temp, entries[end].density, entries[end].pressure, entries[end].mol_wt)
    end

    # Linear interpolation helper
    interp(x, x1, x2, y1, y2) = y1 + (x - x1) / (x2 - x1) * (y2 - y1)

    # Interpolated values
    temp = interp(altitude_km, alts[lower_idx], alts[upper_idx], lower_entry.temp, upper_entry.temp)
    density = interp(altitude_km, alts[lower_idx], alts[upper_idx], lower_entry.density, upper_entry.density)
    pressure = interp(altitude_km, alts[lower_idx], alts[upper_idx], lower_entry.pressure, upper_entry.pressure)
    mol_wt = interp(altitude_km, alts[lower_idx], alts[upper_idx], lower_entry.mol_wt, upper_entry.mol_wt)

    return AtmosphereEntry(altitude_km, temp, density, pressure, mol_wt)
end

"""Overload of atmosphere_characteristics for a vector"""
function atmosphere_characteristics(altitudes_km::AbstractVector, solar_activity::String = "moderate")
    return [atmosphere_characteristics(alt, solar_activity) for alt in altitudes_km]
end



"""
    atmosphere_vallado(altitude)

From Vallado 2013, Table 8-4.
Not tested.

# Arguments
- `altitude::Float64`: Altitude above Earth surface [m]

# Returns
- `density::Float64`: Atmospheric density [kg/m3]

"""
function atmosphere_vallado(h_km::Float64)
    # Ellipsoid Altitude [km], Base density [kg/m3], Scale height [km]
    
    # Find appropriate layer
    h_ellipsoid, base_density, h_scale = vallado_layers[1]
    for layer in vallado_layers
        if h_km >= layer[1]
            h_ellipsoid, base_density, h_scale = layer
        else
            break
        end
    end
    
    base_density * exp(-(h_km - h_ellipsoid) / h_scale)
    return density
end

##############################################################################################################
#           Atmospheric Models
##############################################################################################################



# Low Solar Activity
low_solar_activity = OrderedDict(
    #               alt (km),   temp (K),     density (km/m3),      pressure (Pa),      mol_wt (kg_kmol)
    0 => AtmosphereEntry(    0,      300.2511,       1.17E+00,           1.01E+05,           28.9502 ),
    20 => AtmosphereEntry(    20,     206.2085,       9.48E-02,           5.62E+03,           28.9502 ),
    40 => AtmosphereEntry(    40,     257.6979,       4.07E-03,           3.01E+02,           28.9502 ),
    60 => AtmosphereEntry(    60,     244.1212,       3.31E-04,           2.24E+01,           28.9502 ),
    80 => AtmosphereEntry(    80,     203.1065,       1.69E-05,           9.81E-01,           29.1353 ),
    100 => AtmosphereEntry(    100,    168.7219,       5.77E-07,           2.89E-02,           28.0036 ),
    120 => AtmosphereEntry(    120,    335.8684,       1.70E-08,           1.92E-03,           26.3548 ),
    140 => AtmosphereEntry(    140,    485.8594,       2.96E-09,           4.37E-04,           24.5645 ),
    160 => AtmosphereEntry(    160,    570.0652,       9.65E-10,           1.53E-04,           23.2784 ),
    180 => AtmosphereEntry(    180,    667.8662,       3.90E-10,           9.62E-05,           22.5037 ),
    200 => AtmosphereEntry(    200,    684.9187,       1.75E-10,           4.76E-05,           21.2516 ),
    220 => AtmosphereEntry(    220,    692.6487,       8.47E-11,           2.43E-05,           20.0935 ),
    240 => AtmosphereEntry(    240,    696.1697,       4.31E-11,           1.31E-05,           19.0789 ),
    260 => AtmosphereEntry(    260,    697.7811,       2.30E-11,           7.31E-06,           18.2300 ),
    280 => AtmosphereEntry(    280,    698.5220,       1.27E-11,           4.29E-06,           17.5402 ),
    300 => AtmosphereEntry(    300,    698.8649,       7.24E-12,           1.47E-06,           16.9831 ),
    320 => AtmosphereEntry(    320,    699.0233,       4.21E-12,           1.49E-06,           16.5214 ),
    340 => AtmosphereEntry(    340,    699.0973,       2.50E-12,           9.01E-07,           16.1147 ),
    360 => AtmosphereEntry(    360,    699.1320,       1.51E-12,           5.57E-07,           15.7719 ),
    380 => AtmosphereEntry(    380,    699.1483,       9.20E-13,           3.56E-07,           15.3028 ),
    400 => AtmosphereEntry(    400,    699.1561,       5.68E-13,           2.23E-07,           14.8185 ),
    420 => AtmosphereEntry(    420,    699.1597,       3.54E-13,           1.45E-07,           14.2332 ),
    440 => AtmosphereEntry(    440,    699.1615,       2.23E-13,           9.61E-08,           13.5181 ),
    460 => AtmosphereEntry(    460,    699.1623,       1.42E-13,           6.54E-08,           12.6581 ),
    480 => AtmosphereEntry(    480,    699.1627,       9.20E-14,           4.58E-08,           11.6591 ),
    500 => AtmosphereEntry(    500,    699.1629,       6.03E-14,           3.32E-08,           10.5547 ),
    520 => AtmosphereEntry(    520,    699.1630,       4.03E-14,           2.48E-08,           9.4006  ),
    540 => AtmosphereEntry(    540,    699.1630,       2.75E-14,           1.94E-08,           8.2657  ),
    560 => AtmosphereEntry(    560,    699.1631,       1.93E-14,           1.58E-08,           7.2141  ),
    580 => AtmosphereEntry(    580,    699.1631,       1.39E-14,           1.28E-08,           6.2904  ),
    600 => AtmosphereEntry(    600,    699.1631,       1.03E-14,           1.06E-08,           5.5149  ),
    620 => AtmosphereEntry(    620,    699.1631,       7.96E-15,           9.40E-09,           4.8664  ),
    640 => AtmosphereEntry(    640,    699.1631,       6.24E-15,           8.27E-09,           4.3891  ),
    660 => AtmosphereEntry(    660,    699.1631,       5.06E-15,           7.56E-09,           4.0012  ),
    680 => AtmosphereEntry(    680,    699.1631,       4.21E-15,           6.62E-09,           3.6999  ),
    700 => AtmosphereEntry(    700,    699.1631,       3.58E-15,           6.00E-09,           3.4648  ),
    720 => AtmosphereEntry(    720,    699.1631,       3.09E-15,           5.48E-09,           3.2789  ),
    740 => AtmosphereEntry(    740,    699.1631,       2.70E-15,           5.02E-09,           3.1289  ),
    760 => AtmosphereEntry(    760,    699.1631,       2.39E-15,           4.63E-09,           3.0049  ),
    780 => AtmosphereEntry(    780,    699.1631,       2.13E-15,           4.28E-09,           2.8996  ),
    800 => AtmosphereEntry(    800,    699.1631,       1.91E-15,           3.98E-09,           2.8075  ),
    820 => AtmosphereEntry(    820,    699.1631,       1.73E-15,           3.68E-09,           2.7249  ),
    840 => AtmosphereEntry(    840,    699.1631,       1.56E-15,           3.43E-09,           2.6492  ),
    860 => AtmosphereEntry(    860,    699.1631,       1.42E-15,           3.21E-09,           2.5784  ),
    880 => AtmosphereEntry(    880,    699.1631,       1.30E-15,           3.00E-09,           2.5113  ),
    900 => AtmosphereEntry(    900,    699.1631,       1.18E-15,           2.81E-09,           2.4470  )
)


# Moderate Solar Activity
moderate_solar_activity = OrderedDict(
    #               alt (km),   temp (K),     density (km/m3),      pressure (Pa),      mol_wt (kg_kmol)
    0 => AtmosphereEntry(    0,      300.2511,       1.17E+00,           1.01E+05,           28.9502 ),
    20 => AtmosphereEntry(    20,     206.2085,       9.49E-02,           5.62E+03,           28.9502 ),
    40 => AtmosphereEntry(    40,     257.6979,       4.07E-03,           3.02E+02,           28.9502 ),
    60 => AtmosphereEntry(    60,     244.1212,       3.31E-04,           2.32E+01,           28.9502 ),
    80 => AtmosphereEntry(    80,     196.3636,       1.68E-05,           9.45E-01,           29.0175 ),
    100 => AtmosphereEntry(    100,    194.0160,       5.08E-07,           2.77E-02,           27.7137 ),
    120 => AtmosphereEntry(    120,    374.9775,       1.88E-08,           2.17E-03,           25.8745 ),
    140 => AtmosphereEntry(    140,    565.5703,       3.28E-09,           5.73E-04,           24.3329 ),
    160 => AtmosphereEntry(    160,    767.5532,       1.18E-09,           2.31E-04,           23.1225 ),
    180 => AtmosphereEntry(    180,    877.6229,       5.51E-10,           1.80E-04,           22.4106 ),
    200 => AtmosphereEntry(    200,    931.2806,       2.91E-10,           1.05E-04,           21.4734 ),
    220 => AtmosphereEntry(    220,    963.2701,       1.66E-10,           6.44E-05,           20.6108 ),
    240 => AtmosphereEntry(    240,    982.4191,       9.91E-11,           4.09E-05,           19.8292 ),
    260 => AtmosphereEntry(    260,    993.9173,       6.16E-11,           2.68E-05,           19.1337 ),
    280 => AtmosphereEntry(    280,    1000.8427,      3.94E-11,           1.77E-05,           18.5256 ),
    300 => AtmosphereEntry(    300,    1005.2670,      2.58E-11,           1.20E-05,           18.0037 ),
    320 => AtmosphereEntry(    320,    1007.9620,      1.72E-11,           8.20E-06,           17.5537 ),
    340 => AtmosphereEntry(    340,    1009.0030,      1.16E-11,           5.69E-06,           17.1721 ),
    360 => AtmosphereEntry(    360,    1010.0423,      7.99E-12,           3.98E-06,           16.8449 ),
    380 => AtmosphereEntry(    380,    1010.6166,      5.55E-12,           2.81E-06,           16.5597 ),
    400 => AtmosphereEntry(    400,    1010.9688,      3.89E-12,           2.01E-06,           16.3044 ),
    420 => AtmosphereEntry(    420,    1011.1853,      2.75E-12,           1.44E-06,           16.0669 ),
    440 => AtmosphereEntry(    440,    1011.3190,      1.96E-12,           1.04E-06,           15.8360 ),
    460 => AtmosphereEntry(    460,    1011.4014,      1.40E-12,           7.55E-07,           15.6008 ),
    480 => AtmosphereEntry(    480,    1011.4519,      1.01E-12,           5.53E-07,           15.3509 ),
    500 => AtmosphereEntry(    500,    1011.4845,      7.30E-13,           4.07E-07,           15.0760 ),
    520 => AtmosphereEntry(    520,    1011.5043,      5.31E-13,           3.02E-07,           14.7660 ),
    540 => AtmosphereEntry(    540,    1011.5168,      3.88E-13,           2.27E-07,           14.4148 ),
    560 => AtmosphereEntry(    560,    1011.5245,      2.85E-13,           1.71E-07,           14.0125 ),
    580 => AtmosphereEntry(    580,    1011.5294,      2.11E-13,           1.31E-07,           13.5547 ),
    600 => AtmosphereEntry(    600,    1011.5325,      1.56E-13,           1.01E-07,           13.0389 ),
    620 => AtmosphereEntry(    620,    1011.5345,      1.17E-13,           7.86E-08,           12.4656 ),
    640 => AtmosphereEntry(    640,    1011.5357,      8.79E-14,           6.24E-08,           11.8428 ),
    660 => AtmosphereEntry(    660,    1011.5365,      6.65E-14,           5.01E-08,           11.1779 ),
    680 => AtmosphereEntry(    680,    1011.5370,      5.08E-14,           4.07E-08,           10.4854 ),
    700 => AtmosphereEntry(    700,    1011.5374,      3.91E-14,           3.36E-08,           9.7818  ),
    720 => AtmosphereEntry(    720,    1011.5375,      3.04E-14,           2.82E-08,           9.0847  ),
    740 => AtmosphereEntry(    740,    1011.5377,      2.39E-14,           2.39E-08,           8.4111  ),
    760 => AtmosphereEntry(    760,    1011.5377,      1.90E-14,           2.05E-08,           7.7753  ),
    780 => AtmosphereEntry(    780,    1011.5378,      1.53E-14,           1.78E-08,           7.1894  ),
    800 => AtmosphereEntry(    800,    1011.5378,      1.25E-14,           1.56E-08,           6.6631  ),
    820 => AtmosphereEntry(    820,    1011.5378,      1.03E-14,           1.40E-08,           6.1949  ),
    840 => AtmosphereEntry(    840,    1011.5378,      8.64E-15,           1.26E-08,           5.7711  ),
    860 => AtmosphereEntry(    860,    1011.5379,      7.32E-15,           1.14E-08,           5.4132  ),
    880 => AtmosphereEntry(    880,    1011.5379,      6.28E-15,           1.04E-08,           5.1066  ),
    900 => AtmosphereEntry(    900,    1011.5379,      5.46E-15,           9.47E-09,           4.8460  )
)


# High Solar Activity
high_solar_activity = OrderedDict(
    #               alt (km),   temp (K),     density (km/m3),      pressure (Pa),      mol_wt (kg_kmol)
    0 => AtmosphereEntry(    0,      300.2511,       1.16E+00,           9.98E+04,           28.9502 ),
    20 => AtmosphereEntry(    20,     206.2085,       9.41E-02,           5.57E+03,           28.9502 ),
    40 => AtmosphereEntry(    40,     257.6979,       4.04E-03,           2.99E+02,           28.9502 ),
    60 => AtmosphereEntry(    60,     244.1212,       3.28E-04,           2.30E+01,           28.9502 ),
    80 => AtmosphereEntry(    80,     172.2146,       1.68E-05,           8.42E-01,           28.5290 ),
    100 => AtmosphereEntry(    100,    297.3338,       2.78E-07,           2.47E-02,           26.1997 ),
    120 => AtmosphereEntry(    120,    493.4874,       1.24E-08,           1.96E-03,           24.3506 ),
    140 => AtmosphereEntry(    140,    879.9174,       3.34E-09,           1.01E-03,           22.5096 ),
    160 => AtmosphereEntry(    160,    1143.5236,      2.23E-09,           9.06E-04,           21.3129 ),
    180 => AtmosphereEntry(    180,    1314.3427,      1.28E-09,           6.76E-04,           20.7706 ),
    200 => AtmosphereEntry(    200,    1423.6469,      8.28E-10,           4.86E-04,           20.1836 ),
    220 => AtmosphereEntry(    220,    1493.7864,      5.69E-10,           3.60E-04,           19.6664 ),
    240 => AtmosphereEntry(    240,    1538.9154,      4.08E-10,           2.72E-04,           19.2046 ),
    260 => AtmosphereEntry(    260,    1568.0294,      3.00E-10,           2.08E-04,           18.7901 ),
    280 => AtmosphereEntry(    280,    1586.8613,      2.25E-10,           1.61E-04,           18.4178 ),
    300 => AtmosphereEntry(    300,    1599.0714,      1.71E-10,           1.26E-04,           18.0821 ),
    320 => AtmosphereEntry(    320,    1607.0154,      1.32E-10,           9.93E-05,           17.7852 ),
    340 => AtmosphereEntry(    340,    1612.1920,      1.03E-10,           7.86E-05,           17.5186 ),
    360 => AtmosphereEntry(    360,    1615.5731,      8.05E-11,           6.26E-05,           17.2812 ),
    380 => AtmosphereEntry(    380,    1617.7916,      6.35E-11,           5.01E-05,           17.0669 ),
    400 => AtmosphereEntry(    400,    1619.2476,      5.04E-11,           4.02E-05,           16.8818 ),
    420 => AtmosphereEntry(    420,    1620.2062,      4.02E-11,           3.25E-05,           16.7142 ),
    440 => AtmosphereEntry(    440,    1620.8390,      3.23E-11,           2.63E-05,           16.5643 ),
    460 => AtmosphereEntry(    460,    1621.2577,      2.60E-11,           2.13E-05,           16.4297 ),
    480 => AtmosphereEntry(    480,    1621.5356,      2.10E-11,           1.74E-05,           16.3079 ),
    500 => AtmosphereEntry(    500,    1621.7200,      1.70E-11,           1.42E-05,           16.1967 ),
    520 => AtmosphereEntry(    520,    1621.8420,      1.38E-11,           1.16E-05,           16.0940 ),
    540 => AtmosphereEntry(    540,    1621.9253,      1.13E-11,           9.50E-06,           15.9980 ),
    560 => AtmosphereEntry(    560,    1621.9803,      9.21E-12,           7.81E-06,           15.9067 ),
    580 => AtmosphereEntry(    580,    1622.0172,      7.55E-12,           6.44E-06,           15.8187 ),
    600 => AtmosphereEntry(    600,    1622.0421,      6.20E-12,           5.32E-06,           15.7321 ),
    620 => AtmosphereEntry(    620,    1622.0588,      5.10E-12,           4.40E-06,           15.6457 ),
    640 => AtmosphereEntry(    640,    1622.0702,      4.21E-12,           3.65E-06,           15.5578 ),
    660 => AtmosphereEntry(    660,    1622.0778,      3.47E-12,           3.03E-06,           15.4673 ),
    680 => AtmosphereEntry(    680,    1622.0830,      2.88E-12,           2.52E-06,           15.3725 ),
    700 => AtmosphereEntry(    700,    1622.0865,      2.38E-12,           2.11E-06,           15.2723 ),
    720 => AtmosphereEntry(    720,    1622.0890,      1.98E-12,           1.76E-06,           15.1653 ),
    740 => AtmosphereEntry(    740,    1622.0906,      1.65E-12,           1.48E-06,           15.0503 ),
    760 => AtmosphereEntry(    760,    1622.0918,      1.37E-12,           1.24E-06,           14.9260 ),
    780 => AtmosphereEntry(    780,    1622.0925,      1.15E-12,           1.05E-06,           14.7912 ),
    800 => AtmosphereEntry(    800,    1622.0930,      9.59E-13,           8.84E-07,           14.6446 ),
    820 => AtmosphereEntry(    820,    1622.0934,      8.05E-13,           7.48E-07,           14.4854 ),
    840 => AtmosphereEntry(    840,    1622.0936,      6.74E-13,           6.36E-07,           14.3131 ),
    860 => AtmosphereEntry(    860,    1622.0939,      5.67E-13,           5.42E-07,           14.1244 ),
    880 => AtmosphereEntry(    880,    1622.0940,      4.77E-13,           4.63E-07,           13.9210 ),
    900 => AtmosphereEntry(    900,    1622.0940,      4.03E-13,           3.97E-07,           13.7015 )
)

# Helpers
solar_activity_alts = collect(keys(low_solar_activity)) # Same for all tables


vallado_layers = [
    (0,    1.225,     7.249), # 0-25 km
    (25,   3.899e-2,  6.349), # 25-30 km
    (30,   1.774e-2,  6.682), # 30-40 km
    (40,   3.972e-3,  7.554), # 40-50 km
    (50,   1.057e-3,  8.382), # 50-60 km
    (60,   3.206e-4,  7.714), # 60-70 km
    (70,   8.770e-5,  6.549), # 70-80 km
    (80,   1.905e-5,  5.799), # 80-90 km
    (90,   3.396e-6,  5.382), # 90-100 km
    (100,  5.297e-7,  5.877), # 100-110 km
    (110,  9.661e-8,  7.263), # 110-120 km
    (120,  2.438e-8,  9.473), # 120-130 km
    (130,  8.484e-9,  12.636), # 130-140 km
    (140,  3.845e-9,  16.149), # 140-150 km
    (150,  2.070e-9,  22.523), # 150-180 km
    (180,  5.464e-10, 29.740), # 180-200 km
    (200,  2.789e-10, 37.105), # 200-250 km
    (250,  7.248e-11, 45.546), # 250-300 km
    (300,  2.418e-11, 53.628), # 300-350 km
    (350,  9.518e-12, 53.298), # 350-400 km
    (400,  3.725e-12, 58.515), # 400-450 km
    (450,  1.585e-12, 60.828), # 450-500 km
    (500,  6.967e-13, 63.822), # 500-600 km
    (600,  1.454e-13, 71.835), # 600-700 km
    (700,  3.614e-14, 88.667), # 700-800 km
    (800,  1.170e-14, 124.64), # 800-900 km
    (900,  5.245e-15, 181.05), # 900-1000 km
    (1000, 3.019e-15, 268.00), # 1000+ km
]


