#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class DesignDriver(object):
    """
    Top level design drivers
    """
    INFO = {\
    "design_range":{"unit":"NM", "om":1.e3, "txt":"Range of design mission"},
    "cruise_mach":{"unit":"mach", "om":1.e0, "txt":"Nominal cruise Mach number"},
    "ref_cruise_altp":{"unit":"ft", "om":1.e4, "txt":"Reference cruise altitude (generally 35000ft)"},
    "top_of_climb_altp":{"unit":"ft", "om":1.e4, "txt":"Top of climb altitude (may be lower or equal to reference cruise altitude"}
    }
    def __init__(self, design_range = None,
                       cruise_mach = None,
                       ref_cruise_altp = None,
                       top_of_climb_altp = None):
        self.design_range = design_range
        self.cruise_mach = cruise_mach
        self.ref_cruise_altp = ref_cruise_altp
        self.top_of_climb_altp = top_of_climb_altp


#--------------------------------------------------------------------------------------------------------------------------------
class LowSpeed(object):
    """
    Low speed performance data
    """
    INFO = {\
    "disa_tofl":{"unit":"degK", "om":1.e1, "txt":"Temperature shift for take off field length computation"},
    "altp_tofl":{"unit":"ft", "om":1.e4, "txt":"Altitude for take off field length computation"},
    "kvs1g_tofl":{"unit":"no_dim", "om":1.e0, "txt":"Minimum allowed stall speed margin at take off"},
    "req_tofl":{"unit":"m", "om":1.e3, "txt":"Maximum take off field length at MTOW and given conditions"},
    "eff_tofl":{"unit":"m", "om":1.e3, "txt":"Effective take off field length at MTOW and given condition"},
    "eff_kvs1g":{"unit":"no_dim", "om":1.e0, "txt":"Effective stall speed margin at take off"},
    "seg2_path":{"unit":"no_dim", "om":1.e-1, "txt":"Air path at 35 ft at take off"},
    "limitation":{"unit":"int", "om":1.e0, "txt":"Active limitation, 0: error, 1: field length, 2: min climb path"},
    "perfo_constraint_1":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on Take Off Field Length, must be kept positive"},
    "disa_app_speed":{"unit":"degK", "om":1.e1, "txt":"Temperature shift for approach speed computation"},
    "altp_app_speed":{"unit":"ft", "om":1.e3, "txt":"Altitude for approach speed computation"},
    "kvs1g_app_speed":{"unit":"no_dim", "om":1.e0, "txt":"Minimum allowed stall speed margin at landing"},
    "req_app_speed":{"unit":"kt", "om":1.e2, "txt":"Maximum approach speed at MLW and given conditions"},
    "eff_app_speed":{"unit":"kt", "om":1.e2, "txt":"Effective approach speed at MLW and given condition"},
    "perfo_constraint_2":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on Approach Speed, must be kept positive"},
    "disa_oei":{"unit":"degK", "om":1.e1, "txt":"Temperature shift for One Engine Inoperative (OEI)"},
    "req_oei_altp":{"unit":"ft", "om":1.e4, "txt":"Required One Engine Inoperative (OEI) minimum altitude"},
    "req_oei_path":{"unit":"%", "om":1.e-1, "txt":"Required minimum slope OEI at 95%MTOW, required altitude and MCN rating"},
    "eff_oei_path":{"unit":"%", "om":1.e-1, "txt":"Effective slope OEI at 95%MTOW, required altitude and MCN rating"},
    "oei_best_speed":{"unit":"kt", "om":1.e2, "txt":"Calibrated Air Speed (CAS) at which slope is maximum in given conditions"},
    "perfo_constraint_3":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on One Engine Inoperative performance, must be kept positive"}
    }
    def __init__(self, disa_tofl = None,
                       altp_tofl = None,
                       kvs1g_tofl = None,
                       req_tofl = None,
                       eff_tofl = None,
                       eff_kvs1g = None,
                       seg2_path = None,
                       limitation = None,
                       perfo_constraint_1 = None,
                       disa_app_speed = None,
                       altp_app_speed = None,
                       kvs1g_app_speed = None,
                       req_app_speed = None,
                       eff_app_speed = None,
                       perfo_constraint_2 = None,
                       disa_oei = None,
                       req_oei_altp = None,
                       req_oei_path = None,
                       eff_oei_path = None,
                       oei_best_speed = None,
                       perfo_constraint_3 = None):
        self.disa_tofl = disa_tofl
        self.altp_tofl = altp_tofl
        self.kvs1g_tofl = kvs1g_tofl
        self.req_tofl = req_tofl
        self.eff_tofl = eff_tofl
        self.eff_kvs1g = eff_kvs1g
        self.seg2_path = seg2_path
        self.limitation = limitation
        self.perfo_constraint_1 = perfo_constraint_1
        self.disa_app_speed = disa_app_speed
        self.altp_app_speed = altp_app_speed
        self.kvs1g_app_speed = kvs1g_app_speed
        self.req_app_speed = req_app_speed
        self.eff_app_speed = eff_app_speed
        self.perfo_constraint_2 = perfo_constraint_2
        self.disa_oei = disa_oei
        self.req_oei_altp = req_oei_altp
        self.req_oei_path = req_oei_path
        self.eff_oei_path = eff_oei_path
        self.oei_best_speed = oei_best_speed
        self.perfo_constraint_3 = perfo_constraint_3

#--------------------------------------------------------------------------------------------------------------------------------
class HighSpeed(object):
    """
    High speed performance data
    """
    INFO = {\
    "disa_climb":{"unit":"degK", "om":1.e1, "txt":"Temperature shift for Maximum climb speed computation"},
    "req_vz_climb":{"unit":"ft/min", "om":1.e2, "txt":"Required minimum climb speed at 97%MTOW, nominal initial cruise altitude and MCL rating"},
    "eff_vz_climb":{"unit":"ft/min", "om":1.e2, "txt":"Effective climb speed at 97%MTOW, nominal initial cruise altitude and MCL rating"},
    "perfo_constraint_1":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on climb performance with MCL rating, must be kept positive"},
    "req_vz_cruise":{"unit":"ft/min", "om":1.e2, "txt":"Required minimum climb speed at 97%MTOW, nominal initial cruise altitude and MCR rating"},
    "eff_vz_cruise":{"unit":"ft/min", "om":1.e2, "txt":"Effective climb speed at 97%MTOW, nominal initial cruise altitude and MCR rating"},
    "perfo_constraint_2":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on climb performance with MCR rating, must be kept positive"},
    "req_toc_altp":{"unit":"ft", "om":1.e4, "txt":"Targeted Top Of Climb Altitude (TOC) for Time To Climb (TTC) computation"},
    "cas1_ttc":{"unit":"kt", "om":1.e2, "txt":"Calibrated Air Speed (CAS) below 10000ft for TTC computation"},
    "cas2_ttc":{"unit":"kt", "om":1.e2, "txt":"Calibrated Air Speed (CAS) above 10000ft for TTC computation"},
    "req_ttc":{"unit":"min", "om":1.e1, "txt":"Required maximum Time To Climb"},
    "eff_ttc":{"unit":"min", "om":1.e1, "txt":"Effective Time To Climb"},
    "perfo_constraint_3":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on time to climb, must be kept positive"},
    }
    def __init__(self, disa_climb = None,
                       req_vz_climb = None,
                       eff_vz_climb = None,
                       perfo_constraint_1 = None,
                       req_vz_cruise = None,
                       eff_vz_cruise = None,
                       perfo_constraint_2 = None,
                       req_toc_altp = None,
                       cas1_ttc = None,
                       cas2_ttc = None,
                       req_ttc = None,
                       eff_ttc = None,
                       perfo_constraint_3 = None):
        self.disa_climb = disa_climb
        self.req_vz_climb = req_vz_climb
        self.eff_vz_climb = eff_vz_climb
        self.perfo_constraint_1 = perfo_constraint_1
        self.req_vz_cruise = req_vz_cruise
        self.eff_vz_cruise = eff_vz_cruise
        self.perfo_constraint_2 = perfo_constraint_2
        self.req_toc_altp = req_toc_altp
        self.cas1_ttc = cas1_ttc
        self.cas2_ttc = cas2_ttc
        self.req_ttc = req_ttc
        self.eff_ttc = eff_ttc
        self.perfo_constraint_3 = perfo_constraint_3

#--------------------------------------------------------------------------------------------------------------------------------
class MaxPayloadMission(object):
    """
    Max payload mission data
    """
    INFO = {\
    "range":{"unit":"NM", "om":1.e3, "txt":"Range of the max payload mission"},
    "payload":{"unit":"kg", "om":1.e4, "txt":"Payload of the max payload mission"},
    "tow":{"unit":"kg", "om":1.e4, "txt":"Take off weight of the max payload mission"},
    "total_fuel":{"unit":"kg", "om":1.e4, "txt":"Total fuel of the max payload mission"},
    "block_fuel":{"unit":"kg", "om":1.e4, "txt":"Block fuel of the max payload mission"},
    "block_time":{"unit":"h", "om":1.e1, "txt":"Block time of the max payload mission"},
    "block_enrg":{"unit":"MWh", "om":1.e1, "txt":"Block energy of the max payload mission"},
    "total_enrg":{"unit":"MWh", "om":1.e1, "txt":"Total energy of the max payload mission"},
    "req_battery_mass":{"unit":"kg", "om":1.e3, "txt":"Required battery mass of the max payload mission"}
    }
    def __init__(self, range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_enrg = None,
                       total_enrg = None,
                       req_battery_mass = None):
        self.range = range
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 
        self.block_enrg = block_enrg
        self.total_enrg = total_enrg
        self.req_battery_mass = req_battery_mass

#--------------------------------------------------------------------------------------------------------------------------------
class NominalMission(object):
    """
    Nominal mission data
    """
    INFO = {\
    "range":{"unit":"NM", "om":1.e3, "txt":"Range of the nominal mission"},
    "payload":{"unit":"kg", "om":1.e4, "txt":"Payload of the nominal mission"},
    "nominal_cruise_mach":{"unit":"mach", "om":1.e0, "txt":"Cruise mach of nominal mission"},
    "nominal_cruise_altp":{"unit":"ft", "om":1.e0, "txt":"Cruise altitude of nominal mission"},
    "tow":{"unit":"kg", "om":1.e4, "txt":"Take off weight of the nominal mission"},
    "total_fuel":{"unit":"kg", "om":1.e4, "txt":"Total fuel of the nominal mission"},
    "block_fuel":{"unit":"kg", "om":1.e4, "txt":"Block fuel of the nominal mission"},
    "block_time":{"unit":"h", "om":1.e1, "txt":"Block time of the nominal mission"},
    "block_enrg":{"unit":"MWh", "om":1.e1, "txt":"Block energy of the nominal mission"},
    "total_enrg":{"unit":"MWh", "om":1.e1, "txt":"Total energy of the nominal mission"},
    "req_battery_mass":{"unit":"kg", "om":1.e3, "txt":"Required battery mass of the nominal mission"}
    }
    def __init__(self, range = None,
                       payload = None,
                       nominal_cruise_mach = None,
                       nominal_cruise_altp = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_enrg = None,
                       total_enrg = None,
                       req_battery_mass = None):
        self.range = range
        self.payload = payload
        self.nominal_cruise_mach = nominal_cruise_mach
        self.nominal_cruise_altp = nominal_cruise_altp
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 
        self.block_enrg = block_enrg
        self.total_enrg = total_enrg
        self.req_battery_mass = req_battery_mass

#--------------------------------------------------------------------------------------------------------------------------------
class MaxFuelMission(object):
    """
    Max fuel mission data
    """
    INFO = {\
    "range":{"unit":"NM", "om":1.e3, "txt":"Range of the max fuel mission"},
    "payload":{"unit":"kg", "om":1.e4, "txt":"Payload of the max fuel mission"},
    "tow":{"unit":"kg", "om":1.e4, "txt":"Take off weight of the max fuel mission"},
    "total_fuel":{"unit":"kg", "om":1.e4, "txt":"Total fuel of the max fuel mission"},
    "block_fuel":{"unit":"kg", "om":1.e4, "txt":"Block fuel of the max fuel mission"},
    "block_time":{"unit":"h", "om":1.e1, "txt":"Block time of the max fuel mission"},
    "block_enrg":{"unit":"MWh", "om":1.e1, "txt":"Block energy of the max fuel mission"},
    "total_enrg":{"unit":"MWh", "om":1.e1, "txt":"Total energy of the max fuel mission"},
    "req_battery_mass":{"unit":"kg", "om":1.e3, "txt":"Required battery mass of the max fuel mission"}
    }
    def __init__(self, range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_enrg = None,
                       total_enrg = None,
                       req_battery_mass = None):
        self.range = range
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 
        self.block_enrg = block_enrg
        self.total_enrg = total_enrg
        self.req_battery_mass = req_battery_mass

#--------------------------------------------------------------------------------------------------------------------------------
class ZeroPayloadMission(object):
    """
    Zero payload mission data
    """
    INFO = {\
    "range":{"unit":"NM", "om":1.e3, "txt":"Range of the zero payload mission"},
    "tow":{"unit":"kg", "om":1.e4, "txt":"Take off weight of the zero payload mission"},
    "total_fuel":{"unit":"kg", "om":1.e4, "txt":"Total fuel of the zero payload mission"},
    "block_fuel":{"unit":"kg", "om":1.e4, "txt":"Block fuel of the zero payload mission"},
    "block_time":{"unit":"h", "om":1.e1, "txt":"Block time of the zero payload mission"},
    "block_enrg":{"unit":"MWh", "om":1.e1, "txt":"Block energy of the zero payload mission"},
    "total_enrg":{"unit":"MWh", "om":1.e1, "txt":"Total energy of the zero payload mission"},
    "req_battery_mass":{"unit":"kg", "om":1.e3, "txt":"Required battery mass of the zero payload mission"}
    }
    def __init__(self, range = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_enrg = None,
                       total_enrg = None,
                       req_battery_mass = None):
        self.range = range
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 
        self.block_enrg = block_enrg
        self.total_enrg = total_enrg
        self.req_battery_mass = req_battery_mass

#--------------------------------------------------------------------------------------------------------------------------------
class CostMission(object):
    """
    Mission data for cost evaluation
    """
    INFO = {\
    "disa":{"unit":"degK", "om":1.e1, "txt":"Temperature shift of the cost evaluation mission"},
    "range":{"unit":"NM", "om":1.e3, "txt":"Range of the cost evaluation mission"},
    "payload":{"unit":"kg", "om":1.e4, "txt":"Payload of the cost evaluation mission"},
    "tow":{"unit":"kg", "om":1.e4, "txt":"Take off weight of the cost evaluation mission"},
    "total_fuel":{"unit":"kg", "om":1.e4, "txt":"Total fuel of the cost evaluation mission"},
    "block_fuel":{"unit":"kg", "om":1.e4, "txt":"Block fuel of the cost evaluation mission"},
    "block_time":{"unit":"h", "om":1.e1, "txt":"Block time of the cost evaluation mission"},
    "block_enrg":{"unit":"MWh", "om":1.e1, "txt":"Block energy of the cost evaluation mission"},
    "total_enrg":{"unit":"MWh", "om":1.e1, "txt":"Total energy of the cost evaluation mission"},
    "req_battery_mass":{"unit":"kg", "om":1.e3, "txt":"Required battery mass of the cost evaluation mission"},
    "block_CO2":{"unit":"kg", "om":1.e4, "txt":"Mass of carbon dioxide emitted during the mission"}
    }
    def __init__(self, disa = None,
                       range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_enrg = None,
                       total_enrg = None,
                       req_battery_mass = None,
                       block_CO2 = None):
        self.disa = disa
        self.range = range 
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time
        self.block_enrg = block_enrg
        self.total_enrg = total_enrg
        self.req_battery_mass = req_battery_mass
        self.block_CO2 = block_CO2

#--------------------------------------------------------------------------------------------------------------------------------
class Economics(object):
    """
    Cost data
    """
    INFO = {\
    "gear_price":{"unit":"M$", "om":1.e1, "txt":"Price of landing gears"},
    "engine_price":{"unit":"M$", "om":1.e1, "txt":"Price of one engine"},
    "airplane_price":{"unit":"M$", "om":1.e1, "txt":"Price of the airplane"},
    "battery_price":{"unit":"k$", "om":1.e1, "txt":"Total price of the battery (eventual)"},
    "battery_mass_price":{"unit":"$/kg", "om":1.e1, "txt":"Mass price of battery (eventual)"},
    "fuel_price":{"unit":"$/gal", "om":1.e1, "txt":"Fuel price"},
    "elec_price":{"unit":"$/kWh", "om":1.e-1, "txt":"Price of electricity"},
    "labor_cost":{"unit":"$/h", "om":1.e1, "txt":"Labor cost"},
    "irp":{"unit":"year", "om":1.e1, "txt":"Interest recovery period"},
    "period":{"unit":"year", "om":1.e1, "txt":"Utilisation period"},
    "interest_rate":{"unit":"%", "om":1.e1, "txt":"Interest rate"},
    "utilisation":{"unit":"int", "om":1.e3, "txt":"Number of flights per year"},
    "cockpit_crew_cost":{"unit":"$/trip", "om":1.e3, "txt":"Cockpit crew cost"},
    "cabin_crew_cost":{"unit":"$/trip", "om":1.e3, "txt":"Cabin crew cost"},
    "fuel_cost":{"unit":"$/trip", "om":1.e3, "txt":"Fuel cost"},
    "elec_cost":{"unit":"$/trip", "om":1.e3, "txt":"Cost of electricity"},
    "landing_fees":{"unit":"$/trip", "om":1.e3, "txt":"Landing fees"},
    "navigation_fees":{"unit":"$/trip", "om":1.e3, "txt":"Navigation fees"},
    "catering_cost":{"unit":"$/trip", "om":1.e3, "txt":"Catering cost"},
    "pax_handling_cost":{"unit":"$/trip", "om":1.e3, "txt":"Pax handling cost"},
    "ramp_handling_cost":{"unit":"$/trip", "om":1.e3, "txt":"Ramp handling cost"},
    "standard_operating_cost":{"unit":"$/trip", "om":1.e4, "txt":"Standard operating cost"},
    "cash_operating_cost":{"unit":"$/trip", "om":1.e4, "txt":"Cash operating cost"},
    "total_investment":{"unit":"$/trip", "om":1.e3, "txt":"Total investment"},
    "interest":{"unit":"$/trip", "om":1.e3, "txt":"Interest"},
    "insurance":{"unit":"$/trip", "om":1.e3, "txt":"Insurance"},
    "depreciation":{"unit":"$/trip", "om":1.e3, "txt":"Depreciation"},
    "direct_operating_cost":{"unit":"$/trip", "om":1.e4, "txt":"Direct operating cost"}
    }
    def __init__(self, gear_price = None,
                       engine_price = None,
                       airplane_price = None,
                       battery_price = None,
                       battery_mass_price = None,
                       fuel_price = None,
                       elec_price = None,
                       labor_cost = None,
                       irp = None,
                       period = None,
                       interest_rate = None,
                       utilisation = None,
                       cockpit_crew_cost = None,
                       cabin_crew_cost = None,
                       fuel_cost = None,
                       elec_cost = None,
                       landing_fees = None,
                       navigation_fees = None,
                       catering_cost = None,
                       pax_handling_cost = None,
                       ramp_handling_cost = None,
                       standard_operating_cost = None,
                       cash_operating_cost = None,
                       total_investment = None,
                       interest = None,
                       insurance = None,
                       depreciation = None,
                       direct_operating_cost = None):
        self.gear_price = gear_price
        self.engine_price = engine_price
        self.airplane_price = airplane_price
        self.battery_price = battery_price
        self.battery_mass_price = battery_mass_price
        self.fuel_price = fuel_price
        self.elec_price = elec_price
        self.labor_cost = labor_cost
        self.irp = irp 
        self.period = period 
        self.interest_rate = interest_rate 
        self.utilisation = utilisation
        self.cockpit_crew_cost = cockpit_crew_cost
        self.cabin_crew_cost = cabin_crew_cost 
        self.fuel_cost = fuel_cost
        self.elec_cost = elec_cost
        self.landing_fees = landing_fees
        self.navigation_fees = navigation_fees 
        self.catering_cost = catering_cost 
        self.pax_handling_cost = pax_handling_cost 
        self.ramp_handling_cost = ramp_handling_cost 
        self.standard_operating_cost = standard_operating_cost
        self.cash_operating_cost = cash_operating_cost
        self.total_investment = total_investment
        self.interest = interest
        self.insurance = insurance
        self.depreciation = depreciation
        self.direct_operating_cost = direct_operating_cost 

#--------------------------------------------------------------------------------------------------------------------------------
class Environmental_Impact(object):
    """
    Environmental impact data
    """
    INFO = {\
    "rgf":{"unit":"m2", "om":1.e2, "txt":"Reference Geometric Factor, close to cabin floor pressurized area (but higher)"},
    "CO2_metric":{"unit":"kg/km/m0.48", "om":1.e0, "txt":"Fuel efficiency metric"},
    "CO2_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of carbon dioxide emitted per kg of fuel"},
    "H2O_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of water emitted per kg of fuel"},
    "SO2_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of sulfur dioxide emitted per kg of fuel"},
    "NOx_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of nitrogen oxide emitted per kg of fuel"},
    "CO_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of carbon monoxide emitted per kg of fuel"},
    "HC_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of unburnt hydrocarbon emitted per kg of fuel"},
    "sulfuric_acid_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of sulfuric acid emitted per kg of fuel"},
    "nitrous_acid_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of nitrous acid emitted per kg of fuel"},
    "nitric_acid_index":{"unit":"g/kg", "om":1.e3, "txt":"Mass of nitric acid emitted per kg of fuel"},
    "soot_index":{"unit":"int", "om":1.e12, "txt":"Number of soot particle emitted per kg of fuel"}
    }
    def __init__(self, rgf = None,
                       CO2_metric = None,
                       CO2_index = None,
                       H2O_index = None,
                       SO2_index = None,
                       NOx_index = None,
                       CO_index = None,
                       HC_index = None,
                       sulfuric_acid_index = None,
                       nitrous_acid_index = None,
                       nitric_acid_index = None,
                       soot_index = None):
        self.rgf = rgf
        self.CO2_metric = CO2_metric
        self.CO2_index = CO2_index
        self.H2O_index = H2O_index
        self.SO2_index = SO2_index
        self.NOx_index = NOx_index
        self.CO_index = CO_index
        self.HC_index = HC_index
        self.sulfuric_acid_index = sulfuric_acid_index
        self.nitrous_acid_index = nitrous_acid_index
        self.nitric_acid_index = nitric_acid_index
        self.soot_index = soot_index

