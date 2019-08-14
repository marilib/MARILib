#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

This file contains specific mission function when the airplane is fully powered with batteries,
thus, it does not change its mass along its flight.
In function eval_e_nominal_mission, constraint_3 has been computed considering that the airplane
takes the exact amount of batteries for the required range, eval_e_mission_coupling is written so far.
Functions e_mission_payload, e_mission_tow, e_mission_range, e_mission_energy_limited and eval_cost_e_mission
have been writen with the same approach (embarked battery mass corresponds to required battery mass).
Inportant remark : to keep consistency with fuelled airplanes, battery mass is book kept above OWE
"""


import numpy
from scipy.optimize import fsolve
from marilib.tools import units as unit

from marilib.earth import environment as earth
from marilib.aircraft_model.airplane import aerodynamics as airplane_aero, \
                                            regulation as regul
from marilib.airplane.propulsion import propulsion_models as propu
from marilib.aircraft_model.operations import flight_mechanics as flight


#===========================================================================================================
def b_mission(aircraft,dist_range,tow,altp,mach,disa):
    """
    Mission computation using breguet equation, fixed L/D and fixed sfc
    """

    engine = aircraft.turbofan_engine
    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    g = earth.gravity()

    # Proportion of fuel heat used for propulsion propulsion in ground phases
    # This factor is used to convert into energy the amount of fuel coming from allowance formulas
    e_factor = 0.25 * earth.fuel_heat("Kerosene")

    # Departure ground phases
    #-----------------------------------------------------------------------------------------------------------
    enrg_taxi_out = e_factor*(34. + 2.3e-4*engine.reference_thrust)*propulsion.n_engine
    time_taxi_out = 540.

    enrg_take_off = e_factor*1e-4*(2.8+2.3/engine.bpr)*tow
    time_take_off = 220.*tow/(engine.reference_thrust*propulsion.n_engine)

    # Mission leg
    #-----------------------------------------------------------------------------------------------------------
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    vsnd = earth.sound_speed(tamb)
    tas = vsnd*mach

    throttle = 1.
    nei = 0

    #WARNING : for EF1 architecture SFC is returned in Energy unit per Thrust unit
    if (propulsion.architecture=="EF1"):
        fn,sec,data = propu.ef1_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
    else:
        raise Exception("mission_b, propulsive architecture not allowed")

    mass = 0.95 * tow
    c_z = flight.lift_from_speed(aircraft,pamb,mach,mass)
    c_x,lod_cruise = airplane_aero.drag(aircraft, pamb, tamb, mach, c_z)

    if (propulsion.architecture=="EF1"):
        enrg_mission = (tow*sec*g*dist_range)/(tas*lod_cruise)
    else:
        raise Exception("propulsion.architecture index is not allowed")

    time_mission = 1.09*(dist_range/tas)

    # Arrival ground phases
    #-----------------------------------------------------------------------------------------------------------
    enrg_landing = e_factor*1e-4*(0.5+2.3/engine.bpr)*tow
    time_landing = 180.

    enrg_taxi_in = e_factor*(26. + 1.8e-4*engine.reference_thrust)*propulsion.n_engine
    time_taxi_in = 420.

    # Block fuel and time
    #-----------------------------------------------------------------------------------------------------------
    block_enrg = enrg_taxi_out + enrg_take_off + enrg_mission + enrg_landing + enrg_taxi_in
    block_time = time_taxi_out + time_take_off + time_mission + time_landing + time_taxi_in

    # Diversion fuel
    #-----------------------------------------------------------------------------------------------------------
    enrg_diversion = (tow*sec*g*regul.diversion_range())/(tas*lod_cruise)

    # Holding fuel
    #-----------------------------------------------------------------------------------------------------------
    altp_holding = unit.m_ft(1500.)
    mach_holding = 0.50 * mach
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp_holding,disa)
    lod_max, cz_lod_max = airplane_aero.lod_max(aircraft, pamb, tamb, mach_holding)
    enrg_holding = sec*(tow*g/lod_max)*regul.holding_time()

    # Total
    #-----------------------------------------------------------------------------------------------------------
    design_range = aircraft.design_driver.design_range

    total_enrg = enrg_mission*(1.+regul.reserve_fuel_ratio(design_range)) + enrg_diversion + enrg_holding

    #-----------------------------------------------------------------------------------------------------------
    return block_enrg,block_time,total_enrg


#===========================================================================================================
def b_specific_air_range(aircraft,altp,mass,mach,disa):

    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    g = earth.gravity()

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)

    vsnd = earth.sound_speed(tamb)

    Cz = flight.lift_from_speed(aircraft,pamb,mach,mass)

    [Cx,LoD] = airplane_aero.drag(aircraft,pamb,tamb,mach,Cz)

    thrust = mass*g / LoD

    nei = 0

     #WARNING : for EF1 architecture SFC is returned in Energy unit per Thrust unit
    if (propulsion.architecture=="EF1"):
        sec = propu.ef1_sec(aircraft,pamb,tamb,mach,MCR,thrust,nei)
    else:
        raise Exception("b_specific_air_range, propulsive architecture not allowed")

    sar = (vsnd*mach*LoD)/(mass*g*sec)

    return sar


#===========================================================================================================
def eval_nominal_b_mission(aircraft):
    """
    Compute nominal mission with range as input
    """

    disa = 0.
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0

    aircraft.nominal_mission.payload = aircraft.payload.nominal
    aircraft.nominal_mission.range = aircraft.design_driver.design_range
    aircraft.nominal_mission.tow = aircraft.weights.mtow

    range = aircraft.nominal_mission.range
    tow = aircraft.nominal_mission.tow

    block_enrg,block_time,total_enrg = b_mission(aircraft,range,tow,altp,mach,disa)

    aircraft.nominal_mission.total_enrg = total_enrg
    aircraft.nominal_mission.block_enrg = block_enrg
    aircraft.nominal_mission.block_time = block_time

    payload = aircraft.nominal_mission.payload
    owe = aircraft.weights.owe

    aircraft.nominal_mission.req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    if (aircraft.ef1_battery.stacking=="Variable"):
        mtow = owe + payload + aircraft.nominal_mission.req_battery_mass
    else:
        mtow = owe + payload

    aircraft.weights.mass_constraint_3 = aircraft.weights.mtow - mtow

    return


#===========================================================================================================
def eval_b_mission_coupling(aircraft):
    """
    Mass-Mission coupling
    This relation is put apart from nominal_mission because GEMS does not manage functions that compute their own input
    """

    if (aircraft.ef1_battery.stacking=="Variable"):
        aircraft.weights.mtow = aircraft.weights.owe + aircraft.nominal_mission.payload + aircraft.nominal_mission.req_battery_mass
    else:
        aircraft.weights.mtow = aircraft.weights.owe + aircraft.nominal_mission.payload

    return


#===========================================================================================================
def eval_payload_range_b_missions(aircraft):
    """
    Compute Payload - Range diagram corner points
    """
    disa = 0.
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    # Max payload mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    payload = aircraft.payload.maximum

    aircraft.max_payload_mission.tow = tow
    aircraft.max_payload_mission.payload = payload

    range,block_enrg,block_time,total_enrg = b_mission_range(aircraft,tow,payload,altp,mach,disa)

    aircraft.max_payload_mission.range = range
    aircraft.max_payload_mission.block_enrg = block_enrg
    aircraft.max_payload_mission.block_time = block_time
    aircraft.max_payload_mission.total_enrg = total_enrg

    aircraft.max_payload_mission.req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    # Max fuel mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    total_enrg = aircraft.weights.mfw * aircraft.propulsion.battery_energy_density

    aircraft.max_fuel_mission.tow = tow
    aircraft.max_fuel_mission.total_enrg = total_enrg

    range,payload,block_enrg,block_time = b_mission_energy_limited(aircraft,tow,total_enrg,altp,mach,disa)

    aircraft.max_fuel_mission.payload = payload
    aircraft.max_fuel_mission.range = range
    aircraft.max_fuel_mission.block_enrg = block_enrg
    aircraft.max_fuel_mission.block_time = block_time

    aircraft.max_fuel_mission.req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    # zero payload mission
    #------------------------------------------------------------------------------------------------------
    if (aircraft.ef1_battery.stacking=="Variable"):
        req_battery_mass = min(tow-aircraft.weights.owe, aircraft.weights.mfw)
        total_enrg = req_battery_mass * aircraft.propulsion.battery_energy_density
        tow = aircraft.weights.owe + req_battery_mass
    else:
        total_enrg = aircraft.weights.battery_in_owe * aircraft.propulsion.battery_energy_density
        tow = aircraft.weights.owe

    aircraft.zero_payload_mission.tow = tow
    aircraft.zero_payload_mission.total_enrg = total_enrg

    range,payload,block_enrg,block_time = b_mission_energy_limited(aircraft,tow,total_enrg,altp,mach,disa)

    aircraft.zero_payload_mission.range = range
    aircraft.zero_payload_mission.block_enrg = block_enrg
    aircraft.zero_payload_mission.block_time = block_time

    aircraft.zero_payload_mission.req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    return


#===========================================================================================================
def eval_cost_b_mission(aircraft):
    """
    Compute climb performances
    """

    # Cost mission
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    disa = aircraft.cost_mission.disa
    range = aircraft.cost_mission.range

    payload = aircraft.payload.nominal

    aircraft.cost_mission.payload = payload

    tow,block_enrg,block_time,total_enrg = b_mission_tow(aircraft,payload,range,altp,mach,disa)

    aircraft.cost_mission.tow = tow
    aircraft.cost_mission.total_fuel = 0.
    aircraft.cost_mission.block_fuel = 0.
    aircraft.cost_mission.total_enrg = total_enrg
    aircraft.cost_mission.block_enrg = block_enrg
    aircraft.cost_mission.block_time = block_time

    aircraft.cost_mission.req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    aircraft.cost_mission.block_CO2 = 0.

    return


#===========================================================================================================
def b_mission_payload(aircraft,tow,range,altp,mach,disa):
    """
    Mission simulation (payload weight is output)
    """

    weights = aircraft.weights

    [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)

    if (aircraft.ef1_battery.stacking=="Variable"):
        req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density
        payload = tow - weights.owe - req_battery_mass
    else:
        payload = tow - weights.owe

    return payload,block_enrg,block_time,total_enrg


#===========================================================================================================
def b_mission_tow(aircraft,payload,range,altp,mach,disa):
    """
    Mission simulation (take off weight is output)
    """

    weights = aircraft.weights

    def fct_mission(tow,aircraft,payload,range,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)
        req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density
        Y = tow - weights.owe - payload - req_battery_mass
        return Y
    #---------------------------------------------------------------------------------------

    if (aircraft.ef1_battery.stacking=="Variable"):

        tow_ini = weights.owe + payload + 2000.

        fct_arg = (aircraft,payload,range,altp,mach,disa)

        output_dict = fsolve(fct_mission, x0 = tow_ini, args=fct_arg, full_output=True)

        tow = output_dict[0][0]

        [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)

    else:

        tow = weights.owe + payload

        [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)

    return tow,block_enrg,block_time,total_enrg


#===========================================================================================================
def b_mission_range(aircraft,tow,payload,altp,mach,disa):
    """
    Mission simulation (range is output)
    """

    design_driver = aircraft.design_driver

    def fct_mission(range,aircraft,tow,payload,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)
        req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density
        if (aircraft.ef1_battery.stacking=="Variable"):
            Y = tow - weights.owe - payload - req_battery_mass
        else:
            Y = weights.battery_in_owe - req_battery_mass
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,payload,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    [block_enrg,block_time,total_enrg] = b_mission(aircraft,range,tow,altp,mach,disa)

    return range,block_enrg,block_time,total_enrg


#===========================================================================================================
def b_mission_energy_limited(aircraft,tow,total_enrg,altp,mach,disa):
    """
    Mission fuel limited (range & payload are output)
    """

    design_driver = aircraft.design_driver
    weights = aircraft.weights

    def fct_mission(range,aircraft,tow,total_enrg,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_enrg,block_time,enrg] = b_mission(aircraft,range,tow,altp,mach,disa)
        Y = total_enrg - enrg
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,total_enrg,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    block_enrg,block_time,total_enrg = b_mission(aircraft,range,tow,altp,mach,disa)

    req_battery_mass = total_enrg / aircraft.propulsion.battery_energy_density

    payload = tow - weights.owe - req_battery_mass

    return range,payload,block_enrg,block_time


