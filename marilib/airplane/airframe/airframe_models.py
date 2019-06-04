#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy

from marilib.tools import units as unit


#===========================================================================================================
def  cza_wo_htp(mach,fuselage_width,aspect_ratio,span,sweep):
    """
    Polhamus formula
    """

    cza_wo_htp =  (numpy.pi*aspect_ratio*(1.07*(1+fuselage_width/span)**2)*(1-fuselage_width/span)) \
                 / (1+numpy.sqrt(1+0.25*aspect_ratio**2*(1+numpy.tan(sweep)**2-mach**2)))

    return cza_wo_htp


#===========================================================================================================
def  wing_aero_data(aircraft,mach,hld_conf):

    fuselage = aircraft.fuselage
    wing = aircraft.wing

     # Polhamus formula
    cza_wo_h = cza_wo_htp(mach,fuselage.width,wing.aspect_ratio,wing.span,wing.sweep)

    # Position of wing+fuselage center of lift
    xlc_wo_h = wing.x_mac + (0.25+0.10*hld_conf)*wing.mac

    ki_wing = wing_induced_factor(aircraft)

    return cza_wo_h,xlc_wo_h,ki_wing


#===========================================================================================================
def wing_downwash(aircraft,cz_wing):
    """
    Downwash angle generated behind the wing
    """

    dw_angle = cz_wing*wing_induced_factor(aircraft)

    return dw_angle


#===========================================================================================================
def wing_induced_factor(aircraft):
    """
    Induced drag coefficient of wing + fuselage
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing

    Ki = ((fuselage.width / wing.span)**2 + 1.05 )  / (numpy.pi * wing.aspect_ratio)

    return Ki


#===========================================================================================================
def htp_aero_data(aircraft):
    """
    WARNING : output values are in Wing reference area
    """

    wing = aircraft.wing
    htp = aircraft.horizontal_tail

    # Helmbold formula
    cza_htp =  (numpy.pi*htp.aspect_ratio)/(1+numpy.sqrt(1+(htp.aspect_ratio/2)**2))*(htp.area/wing.area)

    # Position of HTP center of lift
    xlc_htp = htp.x_mac + 0.25*htp.mac

    # Maximum angle of attack allowed for HTP
    aoa_max_htp = unit.rad_deg(9)

    # HTP induced drag coefficient
    ki_htp = 1.2/(numpy.pi*htp.aspect_ratio)

    return cza_htp,xlc_htp,aoa_max_htp,ki_htp


#===========================================================================================================
def vtp_aero_data(aircraft):
    """
    WARNING : output values are in Wing reference area
    """

    vtp = aircraft.vertical_tail

    # Helmbold formula
    cyb_vtp = vtp_lift_gradiant(aircraft) # (numpy.pi*vtp.aspect_ratio)/(1+numpy.sqrt(1+(vtp.aspect_ratio/2)**2))*(vtp.area/wing.area)

    # Position of VTP center of lift
    xlc_vtp = vtp.x_mac + 0.25*vtp.mac

    # Maximum angle of attack allowed for VTP
    aoa_max_vtp = unit.rad_deg(35)

    # VTP induced drag coefficient
    ki_vtp = 1.3/(numpy.pi*vtp.aspect_ratio)

    return cyb_vtp,xlc_vtp,aoa_max_vtp,ki_vtp


#===========================================================================================================
def vtp_lift_gradiant(aircraft):
    """
    Lift curve slope of the vertical tail versus angle of attack
    WARNING : output value is in Wing reference area
    """

    wing = aircraft.wing
    vtp = aircraft.vertical_tail

    # Helmbold formula
    cyb =  (numpy.pi*vtp.aspect_ratio)/(1+numpy.sqrt(1+(vtp.aspect_ratio/2)**2))*(vtp.area/wing.area)

    return cyb


