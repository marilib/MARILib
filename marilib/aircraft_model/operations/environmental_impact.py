#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


from marilib.tools import units as unit

from marilib.aircraft_model.operations import mission as perfo, flight_mechanics as flight


#===========================================================================================================
def fuel_efficiency_metric(aircraft):
    """
    Fuel efficiency metric (CO2 metric)
    """

    cabin = aircraft.cabin
    design_driver = aircraft.design_driver
    weights = aircraft.weights

    rgf = cabin.projected_area      # Reference Geometric Factor (Pressurized floor area)

    high_weight = 0.92*weights.mtow
    low_weight = 0.45*weights.mtow + 0.63*weights.mtow**0.924
    medium_weight = 0.5*(high_weight+low_weight)

    disa = 0.

    mach = design_driver.cruise_mach    # take cruise mach instead of Maxi Range because SFC is constant

    # WARNING : Maximum SAR altitude or speed may be lowered by propulsion ceilings
    #-----------------------------------------------------------------------------------------------------------
    (sar_max_hw,altp_sar_max_hw) = perfo.sar_max(aircraft,high_weight,mach,disa)
    (altp_sar_max_hw,hw_ceiling) = check_ceiling(aircraft,high_weight,altp_sar_max_hw,mach,disa)
    if(hw_ceiling<0.):
        lower_mach = mach - 0.03
        (sar_max_hw,altp_sar_max_hw) = perfo.sar_max(aircraft,high_weight,lower_mach,disa)
        (altp_sar_max_hw,hw_ceiling) = check_ceiling(aircraft,high_weight,altp_sar_max_hw,lower_mach,disa)
        sar_max_hw = perfo.specific_air_range(aircraft,altp_sar_max_hw,high_weight,lower_mach,disa)
    else:
        sar_max_hw = perfo.specific_air_range(aircraft,altp_sar_max_hw,high_weight,mach,disa)

    (sar_max_mw,altp_sar_max_mw) = perfo.sar_max(aircraft,medium_weight,mach,disa)
    (altp_sar_max_mw,mw_ceiling) = check_ceiling(aircraft,medium_weight,altp_sar_max_mw,mach,disa)
    if(mw_ceiling<0.):
        lower_mach = mach - 0.03
        (sar_max_mw,altp_sar_max_mw) = perfo.sar_max(aircraft,medium_weight,lower_mach,disa)
        (altp_sar_max_mw,mw_ceiling) = check_ceiling(aircraft,medium_weight,altp_sar_max_mw,lower_mach,disa)
        sar_max_mw = perfo.specific_air_range(aircraft,altp_sar_max_mw,medium_weight,lower_mach,disa)
    else:
        sar_max_mw = perfo.specific_air_range(aircraft,altp_sar_max_mw,medium_weight,mach,disa)

    (sar_max_lw,altp_sar_max_lw) = perfo.sar_max(aircraft,low_weight,mach,disa)
    (altp_sar_max_lw,lw_ceiling) = check_ceiling(aircraft,low_weight,altp_sar_max_lw,mach,disa)
    if(lw_ceiling<0.):
        lower_mach = mach - 0.03
        (sar_max_lw,altp_sar_max_lw) = perfo.sar_max(aircraft,low_weight,lower_mach,disa)
        (altp_sar_max_lw,lw_ceiling) = check_ceiling(aircraft,low_weight,altp_sar_max_lw,lower_mach,disa)
        sar_max_lw = perfo.specific_air_range(aircraft,altp_sar_max_lw,low_weight,lower_mach,disa)
    else:
        sar_max_lw = perfo.specific_air_range(aircraft,altp_sar_max_lw,low_weight,mach,disa)

    CO2_metric = (1/rgf**0.24)*(1/sar_max_hw + 1/sar_max_mw + 1/sar_max_lw)/3        # kg/m/m2

    return CO2_metric,rgf


#===========================================================================================================
def check_ceiling(aircraft,mass,altp_ini,mach,disa):
    """
    Check reachable altitude
    """

    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    vz_req_mcl = unit.mps_ftpmin(300.)
    vz_req_mcr = unit.mps_ftpmin(0.)

    isomach = 2
    nei = 0

    altp = altp_ini
    ceiling = 0

    [slope,vz_clb] = flight.air_path(aircraft,nei,altp_ini,disa,isomach,mach,mass,MCL)

    if(vz_clb<vz_req_mcl):
        altp, rei = flight.propulsion_ceiling(aircraft,altp_ini,nei,vz_req_mcl,disa,isomach,mach,mass,MCL)
        if(rei==1):
            ceiling = 1
        else:
            ceiling = -1

    [slope,vz_crz] = flight.air_path(aircraft,nei,altp_ini,disa,isomach,mach,mass,MCR)

    if(vz_crz<vz_req_mcr):
        altp, rei = flight.propulsion_ceiling(aircraft,altp_ini,nei,vz_req_mcr,disa,isomach,mach,mass,MCR)

        if(rei==1):
            ceiling = 2
        else:
            ceiling = -1

    return altp,ceiling


