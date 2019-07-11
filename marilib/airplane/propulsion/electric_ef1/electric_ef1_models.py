#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

import numpy

from marilib.earth import environment as earth

from marilib.airplane.propulsion import jet_models as jet

from marilib.airplane.propulsion.turbofan.turbofan_models import turbofan_thrust


#===========================================================================================================
def electrofan_sec(aircraft,pamb,tamb,mach,rating,nei):
    """
    Specific Energy Consumption
    """

    fn,sec = electrofan_thrust(aircraft,pamb,tamb,mach,rating,nei)

    return sec,fn


#===========================================================================================================
def electrofan_thrust(aircraft,pamb,tamb,mach,rating,nei):

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    shaft_power = {"MTO":engine.mto_e_shaft_power,
                   "MCN":engine.mcn_e_shaft_power,
                   "MCL":engine.mcl_e_shaft_power,
                   "MCR":engine.mcr_e_shaft_power,
                   "FID":engine.fid_e_shaft_power}

    pw_shaft = shaft_power[rating]

    (fn_fan,q0) = jet.fan_thrust(nacelle,pamb,tamb,mach,pw_shaft)

    pw_elec = pw_shaft / (nacelle.motor_efficiency*nacelle.controller_efficiency)

    sec = pw_elec / fn_fan

    fn = fn_fan*(engine.n_engine - nei)

    return fn,sec


#===========================================================================================================
def electrofan_nacelle_drag(aircraft,nacelle,Re,Mach):
    """
    Turbofan nacelle drag
    """

    wing = aircraft.wing

    fac = (1. + 0.126*Mach**2)

    # All nacelle drag
    nac_nwa = nacelle.net_wetted_area
    nac_cxf =   1.15*((0.455/fac)*(numpy.log(10)/numpy.log(Re*nacelle.length))**2.58)*nac_nwa/wing.area

    return nac_cxf,nac_nwa


#===========================================================================================================
def electrofan_oei_drag(aircraft,nacelle,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    wing = aircraft.wing

    dCx = 0.12*nacelle.width**2 / wing.area

    return dCx


