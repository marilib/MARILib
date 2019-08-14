#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
"""


from scipy.optimize import fsolve

from marilib.tools.math import maximize_1d

from marilib.aircraft_model.operations import mission_b, mission_f


#===========================================================================================================
def eval_payload_range_missions(aircraft):
    """
    Compute Payload - Range diagram corner points
    """

    if (aircraft.propulsion.fuel_type=="Battery"):
        mission_b.eval_payload_range_b_missions(aircraft)
    else:
        mission_f.eval_payload_range_f_missions(aircraft)

    return


#===========================================================================================================
def eval_nominal_mission(aircraft):
    """
    Nominal mission evaluation
    """

    if (aircraft.propulsion.fuel_type=="Battery"):
        mission_b.eval_nominal_b_mission(aircraft)
    else:
        mission_f.eval_nominal_f_mission(aircraft)


#===========================================================================================================
def eval_cost_mission(aircraft):
    """
    Cost mission evaluation
    """

    if (aircraft.propulsion.fuel_type=="Battery"):
        mission_b.eval_cost_b_mission(aircraft)
    else:
        mission_f.eval_cost_f_mission(aircraft)


#===========================================================================================================
def mission_tow(aircraft,payload,range,altp,mach,disa):
    """
    Mission simulation (take off weight is output)
    """

    if (aircraft.propulsion.fuel_type=="Battery"):
        tow,block_enrg,block_time,total_enrg = mission_b.b_mission_tow(aircraft,payload,range,altp,mach,disa)
        return tow,block_enrg,block_time,total_enrg
    else:
        tow,block_fuel,block_time,total_fuel = mission_f.f_mission_tow(aircraft,payload,range,altp,mach,disa)
        return tow,block_fuel,block_time,total_fuel


#===========================================================================================================
def specific_air_range(aircraft,altp,mass,mach,disa):

    if (aircraft.propulsion.fuel_type=="Battery"):
        sar = mission_b.b_specific_air_range(aircraft,altp,mass,mach,disa)
    else:
        sar = mission_f.f_specific_air_range(aircraft,altp,mass,mach,disa)
    return sar


#===========================================================================================================
def sar_max(aircraft,mass,mach,disa):

    def fct_sar_max(altp,mass,mach,disa,aircraft):
    #=======================================================================================
        sar = specific_air_range(aircraft,altp,mass,mach,disa)
        return sar
    #---------------------------------------------------------------------------------------

    altp_ini = aircraft.design_driver.ref_cruise_altp

    d_altp = 250.

    fct = [fct_sar_max, mass,mach,disa,aircraft]

    (altp_sar_max,sar_max,rc) = maximize_1d(altp_ini,d_altp,fct)

    return sar_max,altp_sar_max


