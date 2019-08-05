#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

import numpy
from scipy.optimize import fsolve

from marilib.earth import environment as earth

from marilib.airplane.propulsion import jet_models as jet

from marilib.airplane.propulsion.turbofan.turbofan_models import turbofan_thrust


#===========================================================================================================
def ef1_sec(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    SFC for EF1 architecture
    """

    #===========================================================================================================
    def fct_sec(throttle,aircraft,pamb,tamb,mach,rating,thrust,nei):
        fn,sec,data = ef1_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)
        return (thrust - fn)

    x_ini = 0.8

    fct_arg = (aircraft,pamb,tamb,mach,rating,thrust,nei)

    output_dict = fsolve(fct_sec, x0=x_ini, args=fct_arg, full_output=True)

    throttle = output_dict[0][0]

    fn,sec,data = ef1_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    return sec


#===========================================================================================================
def ef1_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei):

    propulsion = aircraft.propulsion
    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    shaft_power = {"MTO":engine.mto_e_shaft_power,
                   "MCN":engine.mcn_e_shaft_power,
                   "MCL":engine.mcl_e_shaft_power,
                   "MCR":engine.mcr_e_shaft_power,
                   "FID":engine.fid_e_shaft_power}

    pw_shaft = throttle*shaft_power[rating]

    (fn_fan,q0) = jet.fan_thrust(nacelle,pamb,tamb,mach,pw_shaft)

    pw_elec = pw_shaft / (nacelle.motor_efficiency*nacelle.controller_efficiency)

    if (nacelle.rear_nacelle==1):

        r_engine = aircraft.rear_electric_engine
        r_nacelle = aircraft.rear_electric_nacelle

        r_shaft_power = {"MTO":r_engine.mto_r_shaft_power,
                         "MCN":r_engine.mcn_r_shaft_power,
                         "MCL":r_engine.mcl_r_shaft_power,
                         "MCR":r_engine.mcr_r_shaft_power,
                         "FID":r_engine.fid_r_shaft_power}

        r_pw_shaft = throttle*r_shaft_power[rating]

        if (propulsion.bli_effect>0):
            vsnd = earth.sound_speed(tamb)
            vair = vsnd*mach
            (r_fn_fan,r_q0,dVbli) = jet.fan_thrust_with_bli(r_nacelle,pamb,tamb,mach,r_pw_shaft)
            dvbli_o_v = dVbli/vair
        else:
            (r_fn_fan,r_q0) = jet.fan_thrust(r_nacelle,pamb,tamb,mach,r_pw_shaft)
            dvbli_o_v = 0.

        r_pw_elec = r_pw_shaft / (r_nacelle.motor_efficiency*r_nacelle.controller_efficiency)

    else:

        r_pw_shaft = 0.
        r_pw_elec = 0.
        r_fn_fan = 0.
        r_q0 = 0.
        dvbli_o_v = 0.

    pw = pw_elec*(nacelle.n_engine - nei) + r_pw_elec

    fn = fn_fan*(nacelle.n_engine - nei) + r_fn_fan

    sec = pw / fn

    data = (fn_fan,pw_elec,pw_shaft,q0,r_fn_fan,r_pw_elec,r_pw_shaft,r_q0,dvbli_o_v)

    return fn,sec,data


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


