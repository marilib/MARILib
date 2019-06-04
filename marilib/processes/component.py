#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python

Changed name from "solvers.py" to "component.py" on 21:05:2019
"""


from scipy.optimize import fsolve

from marilib.aircraft_model.operations import mission as perfo, \
                                              environmental_impact as environ, \
                                              pricing_and_costing as costing


#===========================================================================================================
def eval_nominal_mission(aircraft):
    """
    Compute nominal mission with range as input
    """

    disa = 0
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0

    aircraft.nominal_mission.payload = aircraft.payload.nominal
    aircraft.nominal_mission.range = aircraft.design_driver.design_range
    aircraft.nominal_mission.tow = aircraft.weights.mtow

    range = aircraft.nominal_mission.range
    tow = aircraft.nominal_mission.tow

    block_fuel,block_time,total_fuel = perfo.mission(aircraft,range,tow,altp,mach,disa)

    aircraft.nominal_mission.block_fuel = block_fuel
    aircraft.nominal_mission.block_time = block_time
    aircraft.nominal_mission.total_fuel = total_fuel

    payload = aircraft.nominal_mission.payload
    owe = aircraft.weights.owe

    mtow = owe + payload + total_fuel

    aircraft.weights.mass_constraint_3 = aircraft.weights.mtow - mtow

#    aircraft.weights.mtow = mtow    # MTOW is overwritten here

    return


#===========================================================================================================
def mission_payload(aircraft,tow,range,altp,mach,disa):
    """
    Mission simulation (payload weight is output)
    """

    weights = aircraft.weights

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return payload,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_tow(aircraft,payload,range,altp,mach,disa):
    """
    Mission simulation (take off weight is output)
    """

    weights = aircraft.weights

    def fct_mission(tow,aircraft,payload,range,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    tow_ini = weights.owe + payload + 2000

    fct_arg = (aircraft,payload,range,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = tow_ini, args=fct_arg, full_output=True)

    tow = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

    return tow,block_fuel,block_time,total_fuel


#===========================================================================================================
def mission_range(aircraft,tow,payload,altp,mach,disa):
    """
    Mission simulation (range is output)
    """

    design_driver = aircraft.design_driver

    def fct_mission(range,aircraft,tow,payload,altp,mach,disa):
    #=======================================================================================
        weights = aircraft.weights
        [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = tow - weights.owe - payload - total_fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,payload,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    [block_fuel,block_time,total_fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)

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
        [block_fuel,block_time,fuel] = perfo.mission(aircraft,range,tow,altp,mach,disa)
        Y = total_fuel - fuel
        return Y
    #---------------------------------------------------------------------------------------

    range_ini = design_driver.design_range

    fct_arg = (aircraft,tow,total_fuel,altp,mach,disa)

    output_dict = fsolve(fct_mission, x0 = range_ini, args=fct_arg, full_output=True)

    range = output_dict[0][0]

    block_fuel,block_time,total_fuel = perfo.mission(aircraft,range,tow,altp,mach,disa)

    payload = tow - weights.owe - total_fuel

    return range,payload,block_fuel,block_time


#===========================================================================================================
def eval_cost_mission(aircraft):
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

    tow,block_fuel,block_time,total_fuel = mission_tow(aircraft,payload,range,altp,mach,disa)

    aircraft.cost_mission.block_fuel = block_fuel
    aircraft.cost_mission.block_time = block_time
    aircraft.cost_mission.total_fuel = total_fuel

    aircraft.cost_mission.block_CO2 = block_fuel * aircraft.environmental_impact.CO2_index

    return


#===========================================================================================================
def eval_economics(aircraft):
    """
    Compute climb performances
    """

    # Economics
    #------------------------------------------------------------------------------------------------------
    block_fuel = aircraft.cost_mission.block_fuel
    block_time = aircraft.cost_mission.block_time

    direct_op_cost,cash_op_cost,battery_price,gear_price,engine_price,airplane_price = costing.eval_operating_costs(aircraft,block_fuel,block_time)

    aircraft.economics.battery_price = battery_price
    aircraft.economics.gear_price = gear_price
    aircraft.economics.engine_price = engine_price
    aircraft.economics.airplane_price = airplane_price

    aircraft.economics.direct_operating_cost = direct_op_cost
    aircraft.economics.cash_operating_cost = cash_op_cost

    return


#===========================================================================================================
def eval_take_off_performances(aircraft):
    """
    Compute take off field length
    """

    # Take off field length
    #------------------------------------------------------------------------------------------------------
    altp = aircraft.low_speed.altp_tofl
    disa = aircraft.low_speed.disa_tofl
    mass = aircraft.weights.mtow
    hld_conf_to = aircraft.aerodynamics.hld_conf_to

    tofl,seg2_path,eff_kvs1g,limitation = perfo.take_off_field_length(aircraft,altp,disa,mass,hld_conf_to)

    aircraft.low_speed.eff_tofl = tofl
    aircraft.low_speed.eff_kvs1g = eff_kvs1g
    aircraft.low_speed.seg2_path = seg2_path
    aircraft.low_speed.limitation = limitation

    aircraft.low_speed.perfo_constraint_1 = (aircraft.low_speed.req_tofl - tofl) / aircraft.low_speed.req_tofl

    return


#===========================================================================================================
def eval_landing_performances(aircraft):
    """
    Compute approach speed
    """

    # Approach speed
    #------------------------------------------------------------------------------------------------------
    altp = aircraft.low_speed.altp_app_speed
    disa = aircraft.low_speed.disa_app_speed
    mass = aircraft.weights.mlw
    hld_conf_ld = aircraft.aerodynamics.hld_conf_to

    app_speed = perfo.approach_speed(aircraft,altp,disa,mass,hld_conf_ld)

    aircraft.low_speed.eff_app_speed = app_speed

    aircraft.low_speed.perfo_constraint_2 = (aircraft.low_speed.req_app_speed - app_speed) / aircraft.low_speed.req_app_speed

    return


#===========================================================================================================
def eval_co2_metric(aircraft):
    """
    Compute climb performances
    """

    # Environment
    #------------------------------------------------------------------------------------------------------
    CO2_metric,rgf = environ.fuel_efficiency_metric(aircraft)

    aircraft.environmental_impact.rgf = rgf
    aircraft.environmental_impact.CO2_metric = CO2_metric

    return



