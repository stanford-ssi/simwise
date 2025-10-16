# Atmospheric density model
using ..Constants: R_earth


"""
    atmosphere_characteristics(altitude)

From Vallado 2013, Table 8-4

# Arguments
- `altitude::Float64` => Altitude above Earth surface [m]

# Returns
- `density::Float64` => Atmospheric density [kg/m3]
- `T::Float64`: Temperature [K]

"""
function atmosphere_vallado(h_km::Float64)
    # Ellipsoid Altitude [km], Base density [kg/m3], Scale height [km]
    layers = [
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
    
    # Find appropriate layer
    h_ellipsoid, base_density, h_scale = layers[1]
    for layer in layers
        if h_km >= layer[1]
            h_ellipsoid, base_density, h_scale = layer
        else
            break
        end
    end
    
    density = base_density * exp(-(h_km - h_ellipsoid) / h_scale)
    return density
end