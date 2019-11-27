#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

from marilib import numpy

from marilib.tools.math import maximize_1d, lin_interp_1d

from marilib.earth import environment as earth

from marilib.airplane.propulsion import propulsion_models as propu


#===========================================================================================================
def high_lift(wing,hld_conf):
    """
# 0 =< hld_type =< 10
# 0 =< hld_conf =< 1
# Typically : hld_conf = 1 ==> cz_max_ld
#           : hld_conf = 0.1 to 0.5 ==> cz_max_to
    """

    # Maximum lift coefficients of different airfoils, DUBS 1987
    cz_max_ld = {
            0 : 1.45 ,  # Clean
            1 : 2.25 ,  # Flap only, Rotation without slot
            2 : 2.60 ,  # Flap only, Rotation single slot      (ATR)
            3 : 2.80 ,  # Flap only, Rotation double slot
            4 : 2.80 ,  # Fowler Flap
            5 : 2.00 ,  # Slat only
            6 : 2.45 ,  # Slat + Flap rotation without slot
            7 : 2.70 ,  # Slat + Flap rotation single slot
            8 : 2.90 ,  # Slat + Flap rotation double slot
            9 : 3.00 ,  # Slat + Fowler                      (A320)
            10 : 3.20,  # Slat + Fowler + Fowler double slot (A321)
            }.get(wing.hld_type, "Erreur - high_lift_, HLDtype out of range")    # 9 is default if x not found

    if(wing.hld_type<5):
        cz_max_base = 1.45      # Flap only
    else:
        if(hld_conf==0):
            cz_max_base = 1.45  # Clean
        else:
            cz_max_base = 2.00  # Slat + Flap

    cz_max = (1-hld_conf)*cz_max_base + hld_conf*cz_max_ld

    cz_0 = cz_max - cz_max_base  # Assumed the Lift vs AoA is just translated upward and Cz0 clean equal to zero

    return cz_max, cz_0


#===========================================================================================================
def drag(aircraft, pamb, tamb, mach, cz):
    """
    Total aircraft drag with the assumption that the wing takes all the lift
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing
    tanks = aircraft.tanks
    htp = aircraft.horizontal_tail
    vtp = aircraft.vertical_tail

    # Form and friction drag
    #-----------------------------------------------------------------------------------------------------------
    re = earth.reynolds_number(pamb,tamb,mach)

    nac_cxf,nac_nwa = propu.nacelle_drag(aircraft,re,mach)

    fac = ( 1 + 0.126*mach**2 )

    fuse_cxf = 1.05*((0.455/fac)*(numpy.log(10)/numpy.log(re*fuselage.length))**2.58 ) * fuselage.net_wetted_area / wing.area
    pod_cxf = 1.05*((0.455/fac)*(numpy.log(10)/numpy.log(re*tanks.pod_length))**2.58 ) * tanks.pod_net_wetted_area / wing.area
    wing_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*wing.mac))**2.58) * wing.net_wetted_area / wing.area
    htp_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*htp.mac))**2.58) * htp.net_wetted_area / wing.area
    vtp_cxf = 1.4*((0.455/fac)*(numpy.log(10)/numpy.log(re*vtp.mac))**2.58) * vtp.net_wetted_area / wing.area

    aircraft_cxf = fuse_cxf + wing_cxf + htp_cxf + vtp_cxf + nac_cxf + pod_cxf

    # Parasitic drag (seals, antennas, sensors, ...)
    #-----------------------------------------------------------------------------------------------------------
    aircraft_net_wetted_area =   fuselage.net_wetted_area + wing.net_wetted_area + htp.net_wetted_area + vtp.net_wetted_area \
                               + tanks.pod_net_wetted_area + nac_nwa

    aircraft_knwa = aircraft_net_wetted_area/1000.

    aircraft_kp = (0.0247*aircraft_knwa - 0.11)*aircraft_knwa + 0.166       # Parasitic drag factor

    aircraft_cx_parasitic = aircraft_cxf*aircraft_kp

    # Additional drag
    #-----------------------------------------------------------------------------------------------------------
    X = numpy.array([1.0, 1.5, 2.4, 3.3, 4.0, 5.0])
    Y = numpy.array([0.036, 0.020, 0.0075, 0.0025, 0., 0.])

    param = fuselage.tail_cone_length/fuselage.width

    cx_tap_base = lin_interp_1d(param,X,Y)     # Tapered fuselage drag (tail cone)

    fuse_cx_tapered = cx_tap_base*propu.tail_cone_drag_effect(aircraft)     # Effect of tail cone fan

    # Total zero lift drag
    #-----------------------------------------------------------------------------------------------------------
    aircraft_cx0 = aircraft_cxf + aircraft_cx_parasitic + fuse_cx_tapered

    # Induced drag
    #-----------------------------------------------------------------------------------------------------------
    ki = ((fuselage.width / wing.span)**2 + 1.05 )  / (numpy.pi * wing.aspect_ratio)
    aircraft_cxi = ki*cz**2  # Induced drag

    # Compressibility drag
    #-----------------------------------------------------------------------------------------------------------
    # Freely inspired from Korn equation
    cz_design = 0.5
    aircraft_mach_div = 0.03 + design_driver.cruise_mach + 0.1*(cz_design-cz)

    aircraft_cxc = 0.0025 * numpy.exp(40.*(mach - aircraft_mach_div) )

    # Sum up
    #-----------------------------------------------------------------------------------------------------------
    aircraft_cx = aircraft_cx0 + aircraft_cxi + aircraft_cxc

    aircraft_lod  = cz/aircraft_cx

    return aircraft_cx, aircraft_lod


#===========================================================================================================
def lod_max(aircraft,pamb,tamb,mach):
    """
    Maximum lift to drag ratio
    """

    #=======================================================================================
    def fct_lod_max(cz,aircraft,pamb,tamb,mach):
        [cx,lod] = drag(aircraft,pamb,tamb,mach,cz)
        return lod
    #---------------------------------------------------------------------------------------
    cz_ini = 0.5
    dcz = 0.05

    fct = [fct_lod_max, aircraft,pamb,tamb,mach]
    (lod_max_cz,lod_max,rc) = maximize_1d(cz_ini,dcz,fct)

    return lod_max,lod_max_cz


