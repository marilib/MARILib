#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
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
def f_mission(aircraft,dist_range,tow,altp,mach,disa):
    """
    Mission computation using breguet equation, fixed L/D and fixed sfc
    """

    engine = aircraft.turbofan_engine
    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    g = earth.gravity()

    # Departure ground phases
    #-----------------------------------------------------------------------------------------------------------
    fuel_taxi_out = (34. + 2.3e-4*engine.reference_thrust)*propulsion.n_engine
    time_taxi_out = 540.

    fuel_take_off = 1e-4*(2.8+2.3/engine.bpr)*tow
    time_take_off = 220.*tow/(engine.reference_thrust*propulsion.n_engine)

    # Mission leg
    #-----------------------------------------------------------------------------------------------------------
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    vsnd = earth.sound_speed(tamb)
    tas = vsnd*mach

    throttle = 1.
    nei = 0

    if (propulsion.architecture=="TF"):
        fn,sfc,data = propu.turbofan_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
    elif (propulsion.architecture=="PTE1"):
        fn,sfc,sec,data = propu.pte1_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
    else:
        raise Exception("mission_f, propulsive architecture not allowed")

    mass = 0.95 * tow
    c_z = flight.lift_from_speed(aircraft,pamb,mach,mass)
    c_x,lod_cruise = airplane_aero.drag(aircraft, pamb, tamb, mach, c_z)

    if (propulsion.architecture=="TF"):
        fuel_mission = tow*(1-numpy.exp(-(sfc*g*dist_range)/(tas*lod_cruise)))
    elif (propulsion.architecture=="PTE1"):
        fuel_mission = tow*(1-numpy.exp(-(sfc*g*dist_range)/(tas*lod_cruise))) \
                        - (sfc/sec)*aircraft.pte1_battery.energy_cruise
    else:
        raise Exception("propulsion.architecture index is out of range")

    time_mission = 1.09*(dist_range/tas)

    mass = tow - (fuel_taxi_out + fuel_take_off + fuel_mission)

    # Arrival ground phases
    #-----------------------------------------------------------------------------------------------------------
    fuel_landing = 1e-4*(0.5+2.3/engine.bpr)*mass
    time_landing = 180.

    fuel_taxi_in = (26. + 1.8e-4*engine.reference_thrust)*propulsion.n_engine
    time_taxi_in = 420.

    # Block fuel and time
    #-----------------------------------------------------------------------------------------------------------
    block_fuel = fuel_taxi_out + fuel_take_off + fuel_mission + fuel_landing + fuel_taxi_in
    time_block = time_taxi_out + time_take_off + time_mission + time_landing + time_taxi_in

    # Diversion fuel
    #-----------------------------------------------------------------------------------------------------------
    fuel_diversion = mass*(1.-numpy.exp(-(sfc*g*regul.diversion_range())/(tas*lod_cruise)))

    # Holding fuel
    #-----------------------------------------------------------------------------------------------------------
    altp_holding = unit.m_ft(1500.)
    mach_holding = 0.50 * mach
    pamb,tamb,tstd,dtodz = earth.atmosphere(altp_holding,disa)
    lod_max, cz_lod_max = airplane_aero.lod_max(aircraft, pamb, tamb, mach_holding)
    fuel_holding = sfc*(mass*g/lod_max)*regul.holding_time()

    # Total
    #-----------------------------------------------------------------------------------------------------------
    design_range = aircraft.design_driver.design_range

    fuel_total = fuel_mission*(1.+regul.reserve_fuel_ratio(design_range)) + fuel_diversion + fuel_holding

    #-----------------------------------------------------------------------------------------------------------
    return block_fuel,time_block,fuel_total


#===========================================================================================================
def f_specific_air_range(aircraft,altp,mass,mach,disa):

    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    g = earth.gravity()

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)

    vsnd = earth.sound_speed(tamb)

    Cz = flight.lift_from_speed(aircraft,pamb,mach,mass)

    [Cx,LoD] = airplane_aero.drag(aircraft,pamb,tamb,mach,Cz)

    thrust = mass*g / LoD

    nei = 0

    sfc = propu.sfc(aircraft,pamb,tamb,mach,MCR,thrust,nei)

    sar = (vsnd*mach*LoD)/(mass*g*sfc)

    return sar


#===========================================================================================================
def eval_nominal_f_mission(aircraft):
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

    block_fuel,block_time,total_fuel = f_mission(aircraft,range,tow,altp,mach,disa)

    aircraft.nominal_mission.total_fuel = total_fuel
    aircraft.nominal_mission.block_fuel = block_fuel
    aircraft.nominal_mission.block_time = block_time

    payload = aircraft.nominal_mission.payload
    owe = aircraft.weights.owe

    mtow = owe + payload + total_fuel

    aircraft.weights.mass_constraint_3 = aircraft.weights.mtow - mtow

    return


#===========================================================================================================
def eval_f_mission_coupling(aircraft):
    """
    Mass-Mission coupling
    This relation is put apart from nominal_mission because GEMS does not manage functions that compute their own input
    """

    aircraft.weights.mtow = aircraft.weights.owe + aircraft.nominal_mission.payload + aircraft.nominal_mission.total_fuel

    return


#===========================================================================================================
def eval_payload_range_f_missions(aircraft):
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

    range,block_fuel,block_time,total_fuel = f_mission_range(aircraft,tow,payload,altp,mach,disa)

    aircraft.max_payload_mission.range = range
    aircraft.max_payload_mission.block_fuel = block_fuel
    aircraft.max_payload_mission.block_time = block_time
    aircraft.max_payload_mission.total_fuel = total_fuel

    # Max fuel mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    total_fuel = aircraft.weights.mfw

    aircraft.max_fuel_mission.tow = tow
    aircraft.max_fuel_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.max_fuel_mission.payload = payload
    aircraft.max_fuel_mission.range = range
    aircraft.max_fuel_mission.block_fuel = block_fuel
    aircraft.max_fuel_mission.block_time = block_time

    # zero payload mission
    #------------------------------------------------------------------------------------------------------
    total_fuel = aircraft.weights.mfw
    tow = aircraft.weights.owe + total_fuel

    aircraft.zero_payload_mission.tow = tow
    aircraft.zero_payload_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.zero_payload_mission.range = range
    aircraft.zero_payload_mission.block_fuel = block_fuel
    aircraft.zero_payload_mission.block_time = block_time

    return


#===========================================================================================================
def eval_cost_f_mission(aircraft):
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

    tow,block_fuel,block_time,total_fuel = f_mission_tow(aircraft,payload,range,altp,mach,disa)

    aircraft.cost_mission.tow = tow

    aircraft.cost_mission.total_fuel = total_fuel
    aircraft.cost_mission.block_fuel = block_fuel
    aircraft.cost_mission.block_time = block_time

    aircraft.cost_mission.block_CO2 = block_fuel * aircraft.environmental_impact.CO2_index

    return


#===========================================================================================================
def f_mission_payload(aircraft,tow,range,altp,mach,disa):
    """
    Mission simulation (payload weight is output)
    """

    weights = aircraft.weights

    [block_fuel,block_time,total_fuel] = f_mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return payload,block_fuel,block_time,total_fuel


#===========================================================================================================
def f_mission_tow(aircraft,payload,range,altp,mach,disa):
    """
    Mission simulation (take off weight is output)
    """

    weights = aircraft.weights

    def fct_mission(tow,aircraft,payload,range,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = f_mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    tow_ini = weights.owe + payload + 2000.

    fct_arg = (aircraft,payload,range,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = tow_ini, args=fct_arg, full_output=True)

    tow = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = f_mission(aircraft,range,tow,altp,mach,disa)

    return tow,block_fuel,block_time,total_fuel


#===========================================================================================================
def f_mission_range(aircraft,tow,payload,altp,mach,disa):
    """
    Mission simulation (range is output)
    """

    design_driver = aircraft.design_driver

    def fct_mission(range,aircraft,tow,payload,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = f_mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,payload,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = f_mission(aircraft,range,tow,altp,mach,disa)

    return range,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa):
    """
    Mission fuel limited (range & payload are output)
    """

    design_driver = aircraft.design_driver
    weights = aircraft.weights

    def fct_mission(range,aircraft,tow,total_fuel,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,fuel] = f_mission(aircraft,range,tow,altp,mach,disa)
        Y = total_fuel - fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,total_fuel,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    block_fuel,block_time,total_fuel = f_mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return range,payload,block_fuel,block_time


