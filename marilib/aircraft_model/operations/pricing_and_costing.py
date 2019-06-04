#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


import numpy

from marilib.earth import environment as earth


#===========================================================================================================
def landing_gear_price(aircraft):
    """
    Typical value
    """

    ldg = aircraft.landing_gears

    gear_price = 720 * ldg.mass

    return gear_price


#===========================================================================================================
def one_engine_price(aircraft):
    """
    Regression on catalog prices
    """

    enigne = aircraft.turbofan_engine

    engine_price = ((2.115e-4*enigne.reference_thrust + 78.85)*enigne.reference_thrust)

    return engine_price


#===========================================================================================================
def one_airframe_price(aircraft):
    """
    Regression on catalog prices corrected with engine prices
    """

    weights = aircraft.weights
    
    airframe_price = 0.7e3*(9e4 + 1.15*weights.mwe - 1.8e9/(2e4 + weights.mwe**0.94))

    return airframe_price


#===========================================================================================================
def eval_operating_costs(aircraft,block_fuel,block_time):
    """
    Computes Cash and Direct Operating Costs per flight (based on AAE 451 Spring 2004)
    """

    cabin = aircraft.cabin
    propulsion = aircraft.propulsion
    battery = aircraft.battery
    engine = aircraft.turbofan_engine
    weights = aircraft.weights
    cost_mission = aircraft.cost_mission

    eco = aircraft.economics

    labor_cost = eco.labor_cost
    irp = eco.irp
    period = eco.period
    interest_rate = eco.interest_rate
    utilisation = eco.utilisation

    fuel_density = earth.fuel_density(propulsion.fuel_type)

    # Cash Operating Cost
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    eco.fuel_cost =   (block_fuel*(eco.fuel_price*1e3)/fuel_density) \
                    + battery.energy_cruise*eco.elec_price

    b_h = block_time/3600
    t_t = b_h + 0.25
    w_f = (10000 + weights.mwe-propulsion.mass)*1e-5

    labor_frame = ((1.26+1.774*w_f-0.1071*w_f**2)*t_t + (1.614+0.7227*w_f+0.1204*w_f**2))*labor_cost
    matrl_frame = (12.39+29.8*w_f+0.1806*w_f**2)*t_t + (15.20+97.330*w_f-2.8620*w_f**2)
    frame_cost = labor_frame + matrl_frame

    t_h = 0.05*(propulsion.reference_thrust_effective/4.4482198)*1e-4

    labor_engine = engine.n_engine*(0.645*t_t+t_h*(0.566*t_t+0.434))*labor_cost
    matrl_engine = engine.n_engine*(25*t_t+t_h*(0.62*t_t+0.38))
    engine_cost = labor_engine + matrl_engine

    w_g = weights.mtow*1e-3

    eco.cockpit_crew_cost = b_h*2*(440-0.532*w_g)

    eco.cabin_crew_cost = b_h*numpy.ceil(cabin.n_pax_ref/50)*labor_cost

    eco.landing_fees = 8.66*(weights.mtow*1e-3)

    eco.navigation_fees = 57*(cost_mission.range/185200)*numpy.sqrt((weights.mtow/1000)/50)

    eco.catering_cost = 3.07 * cabin.n_pax_ref

    eco.pax_handling_cost = 2 * cabin.n_pax_ref

    eco.ramp_handling_cost = 8.70 * cabin.n_pax_ref

    std_op_cost = eco.fuel_cost + frame_cost + engine_cost + eco.cockpit_crew_cost + eco.landing_fees + eco.navigation_fees

    cash_op_cost = std_op_cost + eco.cabin_crew_cost + eco.catering_cost + eco.pax_handling_cost + eco.ramp_handling_cost

    # DirectOperating Cost
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    engine_price = one_engine_price(aircraft)

    gear_price = landing_gear_price(aircraft)

    frame_price = one_airframe_price(aircraft)

    battery_price = eco.battery_price*battery.mass

    aircraft_price = frame_price + engine_price * engine.n_engine + gear_price + battery_price

    eco.total_investment = frame_price * 1.06 + engine.n_engine * engine_price * 1.025

    eco.interest = (eco.total_investment/(utilisation*period)) * (irp * 0.04 * (((1 + interest_rate)**irp)/((1 + interest_rate)**irp - 1)) - 1)

    eco.insurance = 0.0035 * aircraft_price/utilisation

    eco.depreciation = 0.99 * (eco.total_investment / (utilisation * period))     # Depreciation

    direct_op_cost = cash_op_cost + eco.interest + eco.depreciation + eco.insurance

    return direct_op_cost,cash_op_cost,battery_price,gear_price,engine_price,aircraft_price



