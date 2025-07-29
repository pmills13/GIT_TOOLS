# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 18:03:46 2024

@author: mlfun
"""

import math
from unitConverter import *

# Constants
R = 287.05  # Specific gas constant for dry air, J/(kg·K)
g0_mps2 = 9.80665  # Standard gravitational acceleration, m/s²
gamma = 1.4 #ratio of specific heats

class Airspeed:
    def __init__(self, altitude_ft, airspeed_value=None, airspeed_type=None, units='metric'):
        self.altitude_ft = altitude_ft
        self.airspeed_value = airspeed_value
        self.airspeed_type = airspeed_type
        self.units = units
        self.converter = unitConverter()
        self.atmosphere = self.get_atmospheric_properties()

        if airspeed_value is not None and airspeed_type is not None:
            self.mach = self.calculate_mach()
            self.tas_mps = self.calculate_tas_mps()
            self.ktas = self.calculate_ktas()
            self.keas = self.calculate_keas()
            self.dynamicPressure_psf = self.calculate_dynamicPressure_psf()
            self.dynamicPressure_Pa = self.calculate_dynamicPressure_Pa()

        self.totals = self.get_total_properties()

    def get_atmospheric_properties(self):
        from calculateAtmosphere_v1 import calculateAtmosphere
        T_K, p_Pa, density_kgpm3 = calculateAtmosphere(self.altitude_ft, units='metric')
        return {'T_K': T_K, 'p_Pa': p_Pa, 'density_kgpm3': density_kgpm3}

    def get_total_properties(self):
        Tt_K = self.atmosphere['T_K']*(1+(gamma-1)/2*self.mach**2)
        pt_Pa = self.atmosphere['p_Pa']*(1+(gamma-1)/2*self.mach**2)**(gamma/(gamma-1))
        return {'Tt_K': Tt_K, 'pt_Pa': pt_Pa}

    def calculate_a(self):
        a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
        return a

    def calculate_mach(self):
        if self.airspeed_type == 'Mach':
            return self.airspeed_value
        elif self.airspeed_type == 'trueAirspeed_mps':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return (self.airspeed_value) / a
        elif self.airspeed_type == 'KTAS':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return (self.airspeed_value * 0.514444) / a
        elif self.airspeed_type in ['KEAS']:
            return self.calculate_ktas() / self.calculate_a()
        elif self.airspeed_type == 'dynamicPressure_psf':
            return math.sqrt(2 * self.airspeed_value * 47.8803/self.atmosphere['density_kgpm3']) /  math.sqrt(gamma * R * self.atmosphere['T_K'])
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return math.sqrt(2 * self.airspeed_value/self.atmosphere['density_kgpm3']) / math.sqrt(gamma * R * self.atmosphere['T_K'])
        return None

    def calculate_tas_mps(self):
        if self.airspeed_type == 'Mach':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return self.airspeed_value * a
        elif self.airspeed_type == 'trueAirspeed_mps':
            return self.airspeed_value
        elif self.airspeed_type == 'KTAS':
            return self.airspeed_value*0.514444
        elif self.airspeed_type in ['KEAS']:
            return self.calculate_keas() / math.sqrt(self.atmosphere['density_kgpm3'] / self.calculate_density_sea_level())*0.514444
        elif self.airspeed_type == 'dynamicPressure_psf':
            return math.sqrt((2 * self.airspeed_value * 47.8803) / self.atmosphere['density_kgpm3'])
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return math.sqrt((2 * self.airspeed_value) / self.atmosphere['density_kgpm3'])
        return None

    def calculate_ktas(self):
        if self.airspeed_type == 'Mach':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return self.airspeed_value * a / 0.514444
        elif self.airspeed_type == 'trueAirspeed_mps':
            return self.airspeed_value / 0.514444
        elif self.airspeed_type == 'KTAS':
            return self.airspeed_value
        elif self.airspeed_type in ['KEAS']:
            return self.calculate_keas() / math.sqrt(self.atmosphere['density_kgpm3'] / self.calculate_density_sea_level())
        elif self.airspeed_type == 'dynamicPressure_psf':
            return math.sqrt((2 * self.airspeed_value * 47.8803) / self.atmosphere['density_kgpm3']) / 0.514444
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return math.sqrt((2 * self.airspeed_value) / self.atmosphere['density_kgpm3']) / 0.514444
        return None

    def calculate_keas(self):
        if self.airspeed_type == 'Mach':
            return self.calculate_ktas() * math.sqrt(self.atmosphere['density_kgpm3'] / self.calculate_density_sea_level())
        elif self.airspeed_type == 'KEAS':
            return self.airspeed_value
        elif self.airspeed_type == 'trueAirspeed_mps':
            return self.airspeed_value/0.514444 / math.sqrt(self.calculate_density_sea_level() / self.atmosphere['density_kgpm3'])
        elif self.airspeed_type == 'KTAS':
            return self.airspeed_value / math.sqrt(self.calculate_density_sea_level() / self.atmosphere['density_kgpm3'])
        elif self.airspeed_type == 'dynamicPressure_psf':
            return math.sqrt((2 * self.airspeed_value * 47.8803) / self.calculate_density_sea_level()) / 0.514444
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return math.sqrt((2 * self.airspeed_value) / self.calculate_density_sea_level()) / 0.514444
        return None

    def calculate_dynamicPressure_psf(self):
        if self.airspeed_type == 'Mach':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value * a) ** 2 / 47.8803
        elif self.airspeed_type == 'KEAS':
            return 0.5 * self.calculate_density_sea_level() * (self.airspeed_value * 0.514444) ** 2 / 47.8803
        elif self.airspeed_type == 'trueAirspeed_mps':
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value ) ** 2 / 47.8803
        elif self.airspeed_type == 'KTAS':
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value * 0.514444) ** 2 / 47.8803
        elif self.airspeed_type == 'dynamicPressure_psf':
            return self.airspeed_value
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return self.airspeed_value / 47.8803
        return None

    def calculate_dynamicPressure_Pa(self):
        if self.airspeed_type == 'Mach':
            a = math.sqrt(1.4 * R * self.atmosphere['T_K'])  # Speed of sound
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value * a) ** 2
        elif self.airspeed_type == 'KEAS':
            return 0.5 * self.calculate_density_sea_level() * (self.airspeed_value * 0.514444) ** 2
        elif self.airspeed_type == 'trueAirspeed_mps':
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value ) ** 2
        elif self.airspeed_type == 'KTAS':
            return 0.5 * self.atmosphere['density_kgpm3'] * (self.airspeed_value * 0.514444) ** 2
        elif self.airspeed_type == 'dynamicPressure_psf':
            return self.airspeed_value * 47.8803
        elif self.airspeed_type == 'dynamicPressure_Pa':
            return self.airspeed_value
        return None

    def calculate_density_sea_level(self):
        return 1.225  # kg/m³, standard sea level density

    def convert_to_imperial(self):
        self.atmosphere['T_R'] = self.converter.convert(self.atmosphere['T_K'], 'K', 'R')
        self.atmosphere['p_psf'] = self.converter.convert(self.atmosphere['p_Pa'], 'Pa', 'psf')
        self.atmosphere['density_slugpft3'] = self.converter.convert(self.atmosphere['density_kgpm3'], 'kgpm3', 'slugpft3')
        self.totals['Tt_R'] = self.converter.convert(self.totals['Tt_K'], 'K', 'R')
        self.totals['pt_psf'] = self.converter.convert(self.totals['pt_Pa'], 'Pa', 'psf')

    def get_dynamic_viscosity_kgpms(self,T_K):
        Tref_K = 273.15
        S_K = 110.4
        mu_ref_kgpms = 1.716e-5

        mu_kgpms = mu_ref_kgpms*(T_K/Tref_K)**1.5*(Tref_K+S_K)/(T_K+S_K)

        return mu_kgpms

    def get_results(self):
        if self.units == 'imperial':
            self.convert_to_imperial()
            return {
                'Mach': self.mach,
                'KTAS': self.ktas,
                'trueAirspeed_mps': self.tas_mps,
                'KEAS': self.keas,
                'dynamicPressure_psf': self.dynamicPressure_psf,
                'dynamicPressure_Pa': self.dynamicPressure_Pa,
                'temperature': self.atmosphere['T_R'],
                'pressure': self.atmosphere['p_psf'],
                'density': self.atmosphere['density_slugpft3'],
                'totalTemperature': self.totals['Tt_R'],
                'totalPressure': self.totals['pt_psf'],
                'temperature_units': 'R',
                'pressure_units': 'psf',
                'density_units': 'slug/ft³'
            }
        else:
            self.mu_kgpms = self.get_dynamic_viscosity_kgpms(self.atmosphere['T_K'])

            return {
                'Mach': self.mach,
                'KTAS': self.ktas,
                'trueAirspeed_mps': self.tas_mps,
                'KEAS': self.keas,
                'dynamicPressure_psf': self.dynamicPressure_psf,
                'dynamicPressure_Pa': self.dynamicPressure_Pa,
                'temperature': self.atmosphere['T_K'],
                'pressure': self.atmosphere['p_Pa'],
                'density': self.atmosphere['density_kgpm3'],
                'totalTemperature': self.totals['Tt_K'],
                'totalPressure': self.totals['pt_Pa'],
                'temperature_units': 'K',
                'pressure_units': 'Pa',
                'density_units': 'kg/m³',
                'dynamicViscosity_kgpms':self.mu_kgpms,
                'ReynoldsNumber': self.atmosphere['density_kgpm3']*self.tas_mps/self.mu_kgpms,
                'ReynoldsNumber_units' : '1/m'
            }

## Example usage:
# altitude_ft = 0
# airspeed_value = 340.29228686527705 # Mach
# airspeed_type = 'trueAirspeed_mps'
# units = 'imperial'  # or 'metric'

# results = Airspeed(altitude_ft, airspeed_value, airspeed_type).get_results()

# print(f"Mach: {results['Mach']}")
# print(f"KTAS: {results['KTAS']}")
# print(f"V (m/s): {results['trueAirspeed_mps']}")
# print(f"KEAS: {results['KEAS']}")
# print(f"Dynamic Pressure (psf): {results['dynamicPressure_psf']}")
# print(f"Dynamic Pressure (Pa): {results['dynamicPressure_Pa']}")
# print(f"Temperature: {results['Temperature']} {results['Temperature_units']}")
# print(f"Pressure: {results['Pressure']} {results['Pressure_units']}")
# print(f"Total Temperature: {results['totalTemperature']} {results['Temperature_units']}")
# print(f"Total Pressure: {results['totalPressure']} {results['Pressure_units']}")
