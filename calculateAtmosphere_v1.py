import math
from unitConverter import *

# Constants
R = 287.05  # Specific gas constant for dry air, J/(kg·K)
g0_m_s2 = 9.80665  # Standard gravitational acceleration, m/s²

# Define layers of the atmosphere
layers = [
    {'altitude_m': -1000, 'T_base_K': 288.15, 'p_base_Pa': 101325, 'L_K_m': -0.0065}, #made this layer up. HLV below 0 m
    {'altitude_m': 0, 'T_base_K': 288.15, 'p_base_Pa': 101325, 'L_K_m': -0.0065},
    {'altitude_m': 11000, 'T_base_K': 216.65, 'p_base_Pa': 22632.06, 'L_K_m': 0.0},
    {'altitude_m': 20000, 'T_base_K': 216.65, 'p_base_Pa': 5474.89, 'L_K_m': 0.001},
    {'altitude_m': 32000, 'T_base_K': 228.65, 'p_base_Pa': 868.02, 'L_K_m': 0.0028},
    {'altitude_m': 47000, 'T_base_K': 270.65, 'p_base_Pa': 110.91, 'L_K_m': 0.0},
    {'altitude_m': 51000, 'T_base_K': 270.65, 'p_base_Pa': 66.94, 'L_K_m': -0.0028},
    {'altitude_m': 71000, 'T_base_K': 214.65, 'p_base_Pa': 3.96, 'L_K_m': -0.002},
]

def calculateAtmosphere(altitude_ft, units='metric'):
    converter = unitConverter()
    altitude_m = converter.convert(altitude_ft, 'ft', 'm')  # Convert feet to meters

    for layer in layers:
        if altitude_m < layer['altitude_m']:
            break
        previous_layer = layer

    T_base_K = previous_layer['T_base_K']
    p_base_Pa = previous_layer['p_base_Pa']
    L_K_m = previous_layer['L_K_m']
    h_base_m = previous_layer['altitude_m']

    if L_K_m == 0:
        p_Pa = p_base_Pa * math.exp(-g0_m_s2 * (altitude_m - h_base_m) / (R * T_base_K))
    else:
        p_Pa = p_base_Pa * (T_base_K / (T_base_K + L_K_m * (altitude_m - h_base_m))) ** (g0_m_s2 / (R * L_K_m))

    T_K = T_base_K + L_K_m * (altitude_m - h_base_m)
    density_kgpm3 = p_Pa / (R * T_K)

    if units == 'imperial':
        T_R = converter.convert(T_K, 'K', 'R')
        p_psf = converter.convert(p_Pa, 'Pa', 'psf')
        density_slugpft3 = converter.convert(density_kgpm3, 'kgpm3', 'slugpft3') 
        return T_R, p_psf, density_slugpft3
    else:
        return T_K, p_Pa, density_kgpm3


