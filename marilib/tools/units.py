#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

from marilib import numpy
from copy import deepcopy

      
def s_min(min): return min*60.   # Translate minutes into seconds

def min_s(s): return s/60.   # Translate seconds into minutes

def s_h(h): return h*3600.   # Translate hours into seconds

def h_s(s): return s/3600.   # Translate seconds into hours

def m_ft(ft): return ft*0.3048   # Translate feet into metres

def ft_m(m): return m/0.3048   # Translate metres into feet

def m_NM(NM): return NM*1852.   # Translate nautical miles into metres

def NM_m(m): return m/1852.   # Translate metres into nautical miles

def mps_kmph(kmph): return kmph*1000./3600.   # Translate knots into meters per second

def kmph_mps(mps): return mps*3600./1000.   # Translate knots into meters per second

def mps_kt(kt): return kt*1852/3600   # Translate knots into meters per second

def kt_mps(mps): return mps*3600./1852.   # Translate meters per second into knots

def mps_ftpmin(ftpmin): return ftpmin*0.3048/60.   # Translate feet per minutes into meters per second

def ftpmin_mps(mps): return mps/0.3048*60.   # Translate meters per second into feet per minutes

def liter_usgal(usgal): return usgal*3.7853982   # Translate US gallons into liters

def usgal_liter(liter): return liter/3.7853982   # Translate liters into US gallons

def rad_deg(deg): return deg*numpy.pi/180.   # Translate degrees into radians

def deg_rad(rad): return rad*180./numpy.pi   # Translate radians into degrees

def J_kWh(kWh): return kWh*3.6e6   # Translate kWh into J

def J_MWh(MWh): return MWh*3.6e9   # Translate MWh into J

def kWh_J(J): return J/3.6e6   # Translate J into kWh

def MWh_J(J): return J/3.6e9   # Translate J into MWh


def smart_round(X, S):
    Fac = (10 * numpy.ones(S))**numpy.min(4, max(0, 4 - round(numpy.log10(S))))
    return round(X * Fac)  # Fac


#=========================================================================
#
#	Generic unit converter
#
#=========================================================================
UNIT = {}

# dim = "Distance"
UNIT["m"] = 1.
UNIT["cm"] = 0.01
UNIT["mm"] = 0.001
UNIT["inch"] = 0.0254
UNIT["in"] = 0.0254
UNIT["km"] = 1000.
UNIT["ft"] = 0.3048
UNIT["NM"] = 1852.

# dim = "YearlyDistance"
UNIT["km/year"] = 1.
UNIT["1e12.km/year"] = 1.e-12

# dim = "Area"
UNIT["m2"] = 1.
UNIT["cm2"] = 0.0001
UNIT["ft2"] = 0.0929030
UNIT["inch2"] = 0.00064516
UNIT["in2"] = 0.00064516

# dim = "Duration"
UNIT["s"] = 1.
UNIT["min"] = 60.
UNIT["h"] = 3600.
UNIT["an"] = 31557600.
UNIT["year"] = 31557600.

# dim = "Velocity"
UNIT["m/s"] = 1.
UNIT["ft/s"] = 0.3048
UNIT["ft/min"] = 0.00508
UNIT["km/h"] = 0.2777777778
UNIT["kt"] = 0.5144444444
UNIT["mph"] = 0.4469444

# dim = "Acceleration"
UNIT["m/s2"] = 1.
UNIT["km/s2"] = 1000.
UNIT["ft/s2"] = 0.3048
UNIT["kt/s"] = 0.5144444444

# dim = "AbsoluteTemperature"
UNIT["KELVIN"] = 1.
UNIT["Kelvin"] = 1.
UNIT["K"] = 1.
UNIT["RANKINE"] = 0.5555556
UNIT["Rankine"] = 0.5555556
UNIT["R"] = 0.5555556

# dim = "AbsoluteTemperatureinCelsius"
UNIT["CELSIUS"] = 1.
UNIT["Celsius"] = 1.
UNIT["C"] = 1.

# dim = "DeltaofTemperature"
UNIT["degK"] = 1.
UNIT["degC"] = 1.
UNIT["degF"] = 0.5555556
UNIT["degR"] = 0.5555556

# dim = "Temparaturevariationrate"
UNIT["degK/s"] = 1.
UNIT["degF/s"] = 0.5555556

# dim = "Temparaturegradiant"
UNIT["degK/m"] = 1.
UNIT["degK/km"] = 0.001
UNIT["degF/m"] = 0.5555556
UNIT["degF/km"] = 5.555556e-4

# dim = "DeciBell"
UNIT["dB"] = 1.

# dim = "EffectivePerceivedDeciBell"
UNIT["EPNdB"] = 1.

# dim = "Mass"
UNIT["kg"] = 1.
UNIT["g"] = 0.001
UNIT["lb"] = 0.4535924
UNIT["lbm"] = 0.4535924
UNIT["t"] = 1000.

# dim = "MassIndex"
UNIT["g/kg"] = 1.

# dim = "MasstoForceratio"
UNIT["kg/N"] = 1.
UNIT["g/N"] = 0.001
UNIT["g/kN"] = 0.000001

# dim = "MassperSeat"
UNIT["kg/seat"] = 1.
UNIT["g/seat"] = 0.001

# dim = "MassperSeatandperDistance"
UNIT["kg/m/seat"] = 1.
UNIT["kg/NM/seat"] = 0.000540

# dim = "SpecificConsumptionvsThrust"
UNIT["kg/N/s"] = 1.
UNIT["kg/daN/h"] = 2.77778e-05
UNIT["lb/lbf/h"] = 0.000028327

# dim = "SpecificEnergyConsumption"
UNIT["J/N/s"] = 1.
UNIT["kJ/daN/h"] = 1.e3
UNIT["MJ/lbf/h"] = 1.e6

# dim = "SpecificConsumptionvsPower"
UNIT["kg/W/s"] = 1.
UNIT["kg/kW/h"] = 2.77778e-07
UNIT["lb/shp/h"] = 1.68969e-07

# dim = "Force"
UNIT["N"] = 1.
UNIT["kN"] = 1000.
UNIT["lbf"] = 4.4482198
UNIT["klbf"] = 4448.2198
UNIT["daN"] = 10.
UNIT["kgf"] = 9.8066502

# dim = "Pressure"
UNIT["Pa"] = 1.
UNIT["kPa"] = 1000.
UNIT["MPa"] = 1000000.
UNIT["kgf/m2"] = 9.8066502
UNIT["atm"] = 101325.
UNIT["bar"] = 100000.
UNIT["mbar"] = 100.
UNIT["psi"] = 6895.
UNIT["N/m2"] = 1.
UNIT["daN/m2"] = 10.

# dim = "Pressurevariationrate"
UNIT["Pa/s"] = 1.
UNIT["atm/s"] = 101325.
UNIT["bar/s"] = 100000.

# dim = "VolumetricMass"
UNIT["kg/m3"] = 1.
UNIT["kg/l"] = 1000.
UNIT["lb/ft3"] = 16.018499

# dim = "MassSensitivity"
UNIT["1/kg"] = 1.
UNIT["%/kg"] = 0.01
UNIT["%/ton"] = 0.01 * 0.001

# dim = "VolumetricMassFlow"
UNIT["kg/m3/s"] = 1.
UNIT["lb/ft3/s"] = 16.018499

# dim = "StandardUnit"
UNIT["si"] = 1.
UNIT["std"] = 1.
UNIT["uc"] = 1.
UNIT["cu"] = 1.

# dim = "Angle"
UNIT["rad"] = 1.
UNIT["deg"] = 0.0174533

# dim = "Volume"
UNIT["m3"] = 1.
UNIT["dm3"] = 0.001
UNIT["cm3"] = 0.000001
UNIT["litres"] = 0.001
UNIT["litre"] = 0.001
UNIT["l"] = 0.001
UNIT["ft3"] = 0.0283168

# dim = "VolumeFlow"
UNIT["m3/s"] = 1.
UNIT["litre/s"] = 0.001
UNIT["l/s"] = 0.001
UNIT["ft3/s"] = 0.0283168
UNIT["m3/min"] = 60.
UNIT["litre/min"] = 0.06
UNIT["l/min"] = 0.06
UNIT["ft3/min"] = 1.699008
UNIT["m3/h"] = 3600.
UNIT["litre/h"] = 3.6
UNIT["l/h"] = 3.6
UNIT["ft3/h"] = 101.94048

# dim = "VolumeCoefficient"
UNIT["m2/kN"] = 1.

# dim = "MachNumber"
UNIT["Mach"] = 1.
UNIT["mach"] = 1.

# dim = "DragCount"
UNIT["cx"] = 1.
UNIT["dc"] = 0.0001

# dim = "DragSensitivity"
UNIT["1/cx"] = 1.
UNIT["%/cx"] = 0.01
UNIT["%/dc"] = 0.01 * 10000.

# dim = "MachNumbervariationrate"
UNIT["Mach/s"] = 1.
UNIT["mach/s"] = 1.

# dim = "MassFlow"
UNIT["kg/s"] = 1.
UNIT["kg/h"] = 0.0002778
UNIT["lb/s"] = 0.4535924
UNIT["kg/min"] = 0.0166667
UNIT["lb/min"] = 0.00756
UNIT["lb/h"] = 0.000126

# dim = "Power"
UNIT["Watt"] = 1.
UNIT["W"] = 1.
UNIT["kW"] = 1.e3
UNIT["MW"] = 1.e6
UNIT["GW"] = 1.e9
UNIT["TW"] = 1.e12
UNIT["shp"] = 745.70001

# dim = "PowerDensity"
UNIT["W/kg"] = 1.
UNIT["kW/kg"] = 1.e3

# dim = "PowerDensityPerTime"
UNIT["kW/daN/h"] = 1 / 36.

# dim = "Euro"
UNIT["E"] = 1.
UNIT["Ec."] = 0.01
UNIT["kE"] = 1000.
UNIT["ME"] = 1000000.

# dim = "Cost"
UNIT["$"] = 1.
UNIT["$c."] = 0.01
UNIT["k$"] = 1000.
UNIT["M$"] = 1000000.

# dim = "HourlyCost"
UNIT["$/h"] = 1.

# dim = "Utilisation"
UNIT["trip/year"] = 1.

# dim = "TripCost"
UNIT["$/vol"] = 1.
UNIT["$/trip"] = 1.

# dim = "CosttoDistanceratio"
UNIT["$/km"] = 1.
UNIT["$/NM"] = 0.5399568

# dim = "CostDistancePax"
UNIT["$/km/pax"] = 1.
UNIT["$/NM/pax"] = 0.5399568

# dim = "Lineic"
UNIT["1/m"] = 1.

# dim = "Viscosity"
UNIT["Poises"] = 1.
UNIT["Pl.10e6"] = 0.000001

# dim = "SpecificCost"
UNIT["$/pax/km"] = 1.
UNIT["$/pax/NM"] = 0.5399568

# dim = "InverseAnglular"
UNIT["1/rad"] = 1.
UNIT["1/deg"] = 57.29578

# dim = "InversesquaredAngular"
UNIT["1/rad2"] = 1.

# dim = "MassicDistance"
UNIT["m/kg"] = 1.
UNIT["km/t"] = 1.
UNIT["km/kg"] = 1000.
UNIT["NM/t"] = 1.852
UNIT["NM/kg"] = 1852.
UNIT["NM/lb"] = 4082.8923

# dim = "SurfacicMass"
UNIT["kg/m2"] = 1.
UNIT["lb/ft2"] = 4.8825102

# dim = "LineicMass"
UNIT["kg/m"] = 1.
UNIT["kg/km"] = 0.001
UNIT["lb/m"] = 0.4535924

# dim = "Momentum"
UNIT["N.m"] = 1.
UNIT["daN.m"] = 10.
UNIT["kgf.m"] = 9.8066502
UNIT["lbf.ft"] = 1.3558174

# dim = "InertiaMomentum"
UNIT["kg.m2"] = 1.
UNIT["lb.m2"] = 0.4535924

# dim = "AngularVelocity"
UNIT["rad/s"] = 1.
UNIT["deg/s"] = 0.0174533
UNIT["rpm"] = 0.1047198

# dim = "AngularAcceleration"
UNIT["rad/s2"] = 1.
UNIT["deg/s2"] = 0.0174533
UNIT["rpm/s"] = 0.1047198

# dim = "Energy"
UNIT["J"] = 1.
UNIT["kJ"] = 1.e3
UNIT["MJ"] = 1.e6
UNIT["GJ"] = 1.e9
UNIT["TJ"] = 1.e12
UNIT["Wh"] = 3600.
UNIT["kWh"] = 3600.e3
UNIT["MWh"] = 3600.e6
UNIT["GWh"] = 3600.e9
UNIT["TWh"] = 3600.e12

# dim = "EnergyDensity"
UNIT["J/kg"] = 1.
UNIT["kJ/kg"] = 1.e3
UNIT["MJ/kg"] = 1.e6
UNIT["GJ/kg"] = 1.e9
UNIT["TJ/kg"] = 1.e12
UNIT["Wh/kg"] = 3600.
UNIT["kWh/kg"] = 3600.e3
UNIT["MWh/kg"] = 3600.e6
UNIT["GWh/kg"] = 3600.e9
UNIT["TWh/kg"] = 3600.e12

# dim = "MassicEnergy"
UNIT["J/kg"] = 1.
UNIT["kJ/kg"] = 1000.
UNIT["MJ/kg"] = 1000000.
UNIT["btu/lb"] = 2325.9612

# dim = "FuelCost"
UNIT["$/l"] = 1.
UNIT["$/gal"] = 0.264173
UNIT["$/USgal"] = 0.264173
UNIT["$/USbrl"] = 0.00838644

# dim = "BatteryMassCost"
UNIT["$/kg"] = 1.

# dim = "BatteryEnergyCost"
UNIT["$/kWh"] = 1. / UNIT['kWh']

# dim = "nodimension"
UNIT["sd"] = 1
UNIT["no_dim"] = 1
UNIT["%"] = 0.01
UNIT["%/%"] = 1.

# dim = "integer"
UNIT["integer"] = 1
UNIT["int"] = 1
UNIT["entier"] = 1
UNIT["numeric"] = 1

# dim = "variouscounts"
UNIT["aircraft"] = 1
UNIT["engine"] = 1
UNIT["pilot"] = 1
UNIT["attendant"] = 1
UNIT["trolley"] = 1
UNIT["toilet"] = 1
UNIT["seat"] = 1
UNIT["door"] = 1
UNIT["wheel"] = 1

# dim = "string"
UNIT["string"] = 1
UNIT["text"] = 1

# dim = "textdate"
UNIT["text_date"] = 1

# dim = "GlobalWarmingEnergy"
UNIT["W/m2/km/year"] = 1.
UNIT["1e-6.W/m2/km/year"] = 1.e-6
UNIT["1e-12.W/m2/km/year"] = 1.e-12

# dim = "SpecificMassicEmission"
UNIT["g/seat/m"] = 1.
UNIT["g/seat/km"] = 0.001

# dim = "SpecificVolumicConsumption"
UNIT["m3/seat/m"] = 1.
UNIT["l/seat/100km"] = 0.01

# dim = "CO2metric"
UNIT["kg/NM/m^0.48"] = 1.
UNIT["kg/km/m^0.48"] = 1.852
UNIT["kg/km/m0.48"] = 1.852
UNIT["kg/m/m^0.48"] = 1852.

# dim = "GlobalWarmingTemperature"
UNIT["K/m2/km/year"] = 1.
UNIT["1e-6.K/m2/km/year"] = 1.e-6
UNIT["1e-12.K/m2/km/year"] = 1.e-12

# dim = "DataStructure"
UNIT["structure"] = 1
UNIT["dict"] = 1
UNIT["array"] = 1


# Conversion functions
#-------------------------------------------------------------------------


def convert_from(ulab, val):
    # Convert val expressed in ulab to corresponding standard unit
    if isinstance(val, (type(None), str)):
        return val
    if isinstance(val, list):
        return [convert_from(ulab, v) for v in val]
    if isinstance(val, tuple):
        return (convert_from(ulab, v) for v in val)
    if isinstance(val, numpy.ndarray):
        return numpy.array([convert_from(ulab, v) for v in val])
    if isinstance(val, dict):
        dic_val = deepcopy(val)
        for k, v in dic_val.items():
            dic_val[k] = convert_from(ulab, v)
        return dic_val
    return val * UNIT[ulab]


def convert_to(ulab, val):
    # Convert val expressed in standard unit to ulab
    if isinstance(val, (type(None), str)):
        return val
    if isinstance(val, list):
        return [convert_to(ulab, v) for v in val]
    if isinstance(val, tuple):
        return tuple([convert_to(ulab, v) for v in val])
    if isinstance(val, numpy.ndarray):
        return numpy.array([convert_to(ulab, v) for v in val])
    if isinstance(val, dict):
        dic_val = deepcopy(val)
        for k, v in dic_val.items():
            dic_val[k] = convert_to(ulab, v)
        return dic_val
    return val / UNIT[ulab]
