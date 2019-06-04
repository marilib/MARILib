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
    def __init__(self, design_range = None,
                       cruise_mach = None,
                       ref_cruise_altp = None,
                       top_of_climb_altp = None):
        """
        Constructor :
            :param design_range: Uu NM - OoM 10^3 - Range of design mission
            :param cruise_mach: Uu mach - OoM 10^0 - Nominal cruise Mach number
            :param ref_cruise_altp: Uu ft - OoM 10^4 - Reference cruise altitude (generally 35000ft)
            :param top_of_climb_altp: Uu ft - OoM 10^4 - Top of climb altitude (may be lower or equal to reference cruise altitude
        """
        self.design_range = design_range 
        self.cruise_mach = cruise_mach
        self.ref_cruise_altp = ref_cruise_altp
        self.top_of_climb_altp = top_of_climb_altp

#--------------------------------------------------------------------------------------------------------------------------------
class LowSpeed(object):
    """
    Low speed performance data
    """
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
        """
        Constructor :
            :param disa_tofl: Uu degK - OoM 10^1 - Temperature shift for take off field length computation
            :param altp_tofl: Uu ft - OoM 10^4 - Altitude for take off field length computation
            :param kvs1g_tofl: Uu no_dim - OoM 10^0 - Minimum allowed stall speed margin at take off
            :param req_tofl: Uu m - OoM 10^3 - Maximum take off field length at MTOW and given conditions
            :param eff_tofl: Uu m - OoM 10^3 - Effective take off field length at MTOW and given condition
            :param eff_kvs1g: Uu no_dim - OoM 10^0 - Effective stall speed margin at take off
            :param seg2_path: Uu no_dim - OoM 10^-1 - Air path at 35 ft at take off
            :param limitation: Uu int - OoM 10^0 - Active limitation, 0: error, 1: field length, 2: min climb path
            :param perfo_constraint_1: Uu m - OoM 10^0 - Constraint on Take Off Field Length, must be kept positive
            :param disa_app_speed: Uu degK - OoM 10^1 - Temperature shift for approach speed computation
            :param altp_app_speed: Uu ft - OoM 10^3 - Altitude for approach speed computation
            :param kvs1g_app_speed: Uu no_dim - OoM 10^0 - Minimum allowed stall speed margin at landing
            :param req_app_speed: Uu kt - OoM 10^2 - Maximum approach speed at MLW and given conditions
            :param eff_app_speed: Uu kt - OoM 10^2 - Effective approach speed at MLW and given condition
            :param perfo_constraint_2: Uu m - OoM 10^0 - Constraint on Approach Speed, must be kept positive
            :param disa_oei: Uu degK - OoM 10^1 - Temperature shift for One Engine Inoperative (OEI)
            :param req_oei_altp: Uu ft - OoM 10^4 - Required One Engine Inoperative (OEI) minimum altitude
            :param req_oei_path: Uu % - OoM 10^-1 - Required minimum slope OEI at 95%MTOW, required altitude and MCN rating
            :param eff_oei_path: Uu % - OoM 10^-1 - Effective slope OEI at 95%MTOW, required altitude and MCN rating
            :param oei_best_speed: Uu kt - OoM 10^2 - Calibrated Air Speed (CAS) at which slope is maximum in given conditions
            :param perfo_constraint_3: Uu % - OoM 10^-1 - Constraint on One Engine Inoperative performance, must be kept positive
        """
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
                       perfo_constraint_3 = None,
                       cruise_sfc = None,
                       cruise_lod = None):
        """
        Constructor :
            :param disa_climb: Uu degK - OoM 10^1 - Temperature shift for Maximum climb speed computation
            :param req_vz_climb: Uu ft/min - OoM 10^2 - Required minimum climb speed at 97%MTOW, nominal initial cruise altitude and MCL rating
            :param eff_vz_climb: Uu ft/min - OoM 10^2 - Effective climb speed at 97%MTOW, nominal initial cruise altitude and MCL rating
            :param perfo_constraint_1: Uu ft/min - OoM 10^0 - Constraint on climb performance with MCL rating, must be kept positive
            :param req_vz_cruise: Uu ft/min - OoM 10^2 - Required minimum climb speed at 97%MTOW, nominal initial cruise altitude and MCR rating
            :param eff_vz_cruise: Uu ft/min - OoM 10^2 - Effective climb speed at 97%MTOW, nominal initial cruise altitude and MCR rating
            :param perfo_constraint_2: Uu ft/min - OoM 10^0 - Constraint on climb performance with MCR rating, must be kept positive
            :param req_toc_altp: Uu ft - OoM 10^4 - Targeted Top Of Climb Altitude (TOC) for Time To Climb (TTC) computation
            :param cas1_ttc: Uu kt - OoM 10^2 - Calibrated Air Speed (CAS) below 10000ft for TTC computation
            :param cas2_ttc: Uu kt - OoM 10^2 - Calibrated Air Speed (CAS) above 10000ft for TTC computation
            :param req_ttc: Uu min - OoM 10^1 - Required maximum Time To Climb
            :param eff_ttc: Uu min - OoM 10^1 - Effective Time To Climb
            :param perfo_constraint_3: Uu min - OoM 10^0 - Constraint on time to climb, must be kept positive
            :param cruise_sfc: Uu kg/daN/h - OoM 10^0 - Specific fuel consumption for nominal mission cruise
            :param cruise_lod: Uu no_dim - OoM 10^1 - Lift over drag ratio for nominal mission cruise
        """
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
        self.cruise_sfc = cruise_sfc
        self.cruise_lod = cruise_lod

#--------------------------------------------------------------------------------------------------------------------------------
class MaxPayloadMission(object):
    """
    Max payload mission data
    """
    def __init__(self, range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None):
        """
        Constructor :
            :param range: Uu NM - OoM 10^3 - Range of the max payload mission
            :param payload: Uu kg - OoM 10^4 - Payload of the max payload mission
            :param tow: Uu kg - OoM 10^4 - Take off weight of the max payload mission
            :param total_fuel: Uu kg - OoM 10^4 - Total fuel of the max payload mission
            :param block_fuel: Uu kg - OoM 10^4 - Block fuel of the max payload mission
            :param block_time: Uu h - OoM 10^1 - Block time of the max payload mission
        """
        self.range = range 
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 

#--------------------------------------------------------------------------------------------------------------------------------
class NominalMission(object):
    """
    Nominal mission data
    """
    def __init__(self, range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None):
        """
        Constructor :
            :param range: Uu NM - OoM 10^3 - Range of the nominal mission
            :param payload: Uu kg - OoM 10^4 - Payload of the nominal mission
            :param tow: Uu kg - OoM 10^4 - Take off weight of the nominal mission
            :param total_fuel: Uu kg - OoM 10^4 - Total fuel of the nominal mission
            :param block_fuel: Uu kg - OoM 10^4 - Block fuel of the nominal mission
            :param block_time: Uu h - OoM 10^1 - Block time of the nominal mission
        """
        self.range = range 
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 

#--------------------------------------------------------------------------------------------------------------------------------
class MaxFuelMission(object):
    """
    Max fuel mission data
    """
    def __init__(self, range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None):
        """
        Constructor :
            :param range: Uu NM - OoM 10^3 - Range of the max fuel mission
            :param payload: Uu kg - OoM 10^4 - Payload of the max fuel mission
            :param tow: Uu kg - OoM 10^4 - Take off weight of the max fuel mission
            :param total_fuel: Uu kg - OoM 10^4 - Total fuel of the max fuel mission
            :param block_fuel: Uu kg - OoM 10^4 - Block fuel of the max fuel mission
            :param block_time: Uu h - OoM 10^1 - Block time of the max fuel mission
        """
        self.range = range 
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 

#--------------------------------------------------------------------------------------------------------------------------------
class ZeroPayloadMission(object):
    """
    Zero payload mission data
    """
    def __init__(self, range = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None):
        """
        Constructor :
            :param range: Uu NM - OoM 10^3 - Range of the zero payload mission
            :param tow: Uu kg - OoM 10^4 - Take off weight of the zero payload mission
            :param total_fuel: Uu kg - OoM 10^4 - Total fuel of the zero payload mission
            :param block_fuel: Uu kg - OoM 10^4 - Block fuel of the zero payload mission
            :param block_time: Uu h - OoM 10^1 - Block time of the zero payload mission
        """
        self.range = range 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time 

#--------------------------------------------------------------------------------------------------------------------------------
class CostMission(object):
    """
    Mission data for cost evaluation
    """
    def __init__(self, disa = None,
                       range = None,
                       payload = None,
                       tow = None,
                       total_fuel = None,
                       block_fuel = None,
                       block_time = None,
                       block_CO2 = None):
        """
        Constructor :
            :param disa: Uu degK - OoM 10^1 - Temperature shift of the cost evaluation mission
            :param range: Uu NM - OoM 10^3 - Range of the cost evaluation mission
            :param payload: Uu kg - OoM 10^4 - Payload of the cost evaluation mission
            :param tow: Uu kg - OoM 10^4 - Take off weight of the cost evaluation mission
            :param total_fuel: Uu kg - OoM 10^4 - Total fuel of the cost evaluation mission
            :param block_fuel: Uu kg - OoM 10^4 - Block fuel of the cost evaluation mission
            :param block_time: Uu h - OoM 10^1 - Block time of the cost evaluation mission
            :param block_CO2: Uu kg - OoM 10^4 - Mass of carbon dioxide emitted during the mission
        """
        self.disa = disa 
        self.range = range 
        self.payload = payload 
        self.tow = tow  
        self.total_fuel = total_fuel 
        self.block_fuel = block_fuel 
        self.block_time = block_time
        self.block_CO2 = block_CO2

#--------------------------------------------------------------------------------------------------------------------------------
class Economics(object):
    """
    Cost data
    """
    def __init__(self, gear_price = None,
                       engine_price = None,
                       battery_price = None,
                       airplane_price = None,
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
        """
        Constructor :
            :param gear_price: Uu M$ - OoM 10^1 - Price of landing gears 
            :param engine_price: Uu M$ - OoM 10^1 - Price of one engine 
            :param battery_price: Uu $/kg - OoM 10^1 - Mass price of battery (eventual)
            :param airplane_price: Uu M$ - OoM 10^1 - Price of the airplane
            :param fuel_price: Uu $/gal - OoM 10^1 - Fuel price 
            :param elec_price: Uu $/kWh - OoM 10^-1 - Price of electricity
            :param labor_cost: Uu $/h - OoM 10^1 - Labor cost
            :param irp: Uu year - OoM 10^1 - Interest recovery period 
            :param period: Uu year - OoM 10^1 - Utilisation period 
            :param interest_rate: Uu % - OoM 10^1 - Interest rate 
            :param utilisation: Uu int - OoM 10^3 - Number of flights per year
            :param cockpit_crew_cost: Uu $/trip - OoM 10^3 - Cockpit crew cost
            :param cabin_crew_cost: Uu $/trip - OoM 10^3 - Cabin crew cost
            :param fuel_cost: Uu $/trip - OoM 10^3 - Fuel cost
            :param landing_fees: Uu $/trip - OoM 10^3 - Landing fees
            :param navigation_fees: Uu $/trip - OoM 10^3 - Navigation fees
            :param catering_cost: Uu $/trip - OoM 10^3 - Catering cost
            :param pax_handling_cost: Uu $/trip - OoM 10^3 - Pax handling cost
            :param ramp_handling_cost: Uu $/trip - OoM 10^3 - Ramp handling cost
            :param standard_operating_cost: Uu $/trip - OoM 10^4 - Standard operating cost
            :param cash_operating_cost: Uu $/trip - OoM 10^4 - Cash operating cost
            :param total_investment: Uu $/trip - OoM 10^3 - Total investment
            :param interest: Uu $/trip - OoM 10^3 - Interest
            :param insurance: Uu $/trip - OoM 10^3 - Insurance
            :param depreciation: Uu $/trip - OoM 10^3 - Depreciation
            :param direct_operating_cost: Uu $/trip - OoM 10^4 - Direct operating cost
        """
        self.gear_price = gear_price 
        self.engine_price = engine_price 
        self.battery_price = battery_price
        self.airplane_price = airplane_price
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
    def __init__(self, rgf = None,
                       CO2_metric = None,
                       CO2_index = None,
                       H2O_index = None,
                       SO2_index = None,
                       NOx_index = None,
                       CO_index = None,
                       HC_index = None,
                       sulphuric_acid_index = None,
                       nitrous_acid_index = None,
                       nitric_acid_index = None,
                       soot_index = None):
        """
        Constructor :
            :param rgf: Uu m2 - OoM 10^2 - Reference Geometric Factor, close to cabin floor pressurized area (but higher)
            :param CO2_metric: Uu kg/km/m0.48 - OoM 10^0 - Fuel efficiency metric
            :param CO2_index: Uu g/kg - OoM 10^3 - Mass of carbon dioxide emitted per kg of fuel
            :param H2O_index: Uu g/kg - OoM 10^3 - Mass of water emitted per kg of fuel
            :param SO2_index: Uu g/kg - OoM 10^3 - Mass of sulfur dioxide emitted per kg of fuel
            :param NOx_index: Uu g/kg - OoM 10^3 - Mass of nitrogen oxide emitted per kg of fuel
            :param CO_index: Uu g/kg - OoM 10^3 - Mass of carbon monoxide emitted per kg of fuel
            :param HC_index: Uu g/kg - OoM 10^3 - Mass of unburnt hydrocarbon emitted per kg of fuel
            :param sulphuric_acid_index: Uu g/kg - OoM 10^3 - Mass of sulfuric acid emitted per kg of fuel
            :param nitrous_acid_index: Uu g/kg - OoM 10^3 - Mass of nitrous acid emitted per kg of fuel
            :param nitric_acid_index: Uu g/kg - OoM 10^3 - Mass of nitric acid emitted per kg of fuel
            :param soot_index: Uu int - OoM 10^12 - Number of soot particle emitted per kg of fuel
        """
        self.rgf = rgf
        self.CO2_metric = CO2_metric
        self.H2O_index = H2O_index
        self.SO2_index = SO2_index
        self.NOx_index = NOx_index
        self.CO_index = CO_index
        self.HC_index = HC_index
        self.sulphuric_acid_index = sulphuric_acid_index
        self.nitrous_acid_index = nitrous_acid_index
        self.nitric_acid_index = nitric_acid_index
        self.soot_index = soot_index

