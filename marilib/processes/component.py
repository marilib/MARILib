#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python

Changed name from "solvers.py" to "component.py" on 21:05:2019
"""


from marilib.tools import units as unit

from marilib.aircraft_model.operations import other_performances as perfo, \
                                              environmental_impact as environ, \
                                              pricing_and_costing as costing


#===========================================================================================================
def eval_climb_performances(aircraft):
    """
    Compute climb performances
    """

    # Ceilings
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.design_driver.top_of_climb_altp
    oei_ceil_req = aircraft.low_speed.req_oei_altp

    vz_clb,vz_crz,oei_path,oei_speed = perfo.ceilings(aircraft,toc,oei_ceil_req)

    aircraft.low_speed.eff_oei_path = oei_path
    aircraft.low_speed.oei_best_speed = oei_speed
    aircraft.high_speed.eff_vz_climb = vz_clb
    aircraft.high_speed.eff_vz_cruise = vz_crz

    aircraft.low_speed.perfo_constraint_3 = (oei_path - aircraft.low_speed.req_oei_path) / aircraft.low_speed.req_oei_path

    aircraft.high_speed.perfo_constraint_1 = (vz_clb - aircraft.high_speed.req_vz_climb) / unit.mps_ftpmin(300.)
    aircraft.high_speed.perfo_constraint_2 = (vz_crz - aircraft.high_speed.req_vz_cruise) / unit.mps_ftpmin(300.)

    # Time to climb to requested altitude
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.high_speed.req_toc_altp
    disa = 0.
    mass = aircraft.weights.mtow
    vcas1 = aircraft.high_speed.cas1_ttc
    vcas2 = aircraft.high_speed.cas2_ttc
    mach = aircraft.design_driver.cruise_mach

    ttc = perfo.time_to_climb(aircraft,toc,disa,mass,vcas1,vcas2,mach)

    aircraft.high_speed.eff_ttc = ttc

    aircraft.high_speed.perfo_constraint_3 = (aircraft.high_speed.req_ttc - ttc) / aircraft.high_speed.req_ttc

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

    standard_op_cost,direct_op_cost,cash_op_cost,battery_price,gear_price,engine_price,airplane_price = costing.eval_operating_costs(aircraft,block_fuel,block_time)

    aircraft.economics.battery_price = battery_price
    aircraft.economics.gear_price = gear_price
    aircraft.economics.engine_price = engine_price
    aircraft.economics.airplane_price = airplane_price

    aircraft.economics.standard_operating_cost = standard_op_cost
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

    cabin = aircraft.cabin

    rgf = cabin.projected_area      # Reference Geometric Factor (Pressurized floor area)

    # Environment
    #------------------------------------------------------------------------------------------------------
    if (aircraft.propulsion.fuel_type=="Battery"):
        CO2_metric = 0.
    else:
        CO2_metric = environ.fuel_efficiency_metric(aircraft)

    aircraft.environmental_impact.rgf = rgf
    aircraft.environmental_impact.CO2_metric = CO2_metric

    return



