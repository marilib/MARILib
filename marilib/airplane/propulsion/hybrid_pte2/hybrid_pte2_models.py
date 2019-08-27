#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

from marilib import numpy
from marilib import fsolve

from marilib.earth import environment as earth

from marilib.airplane.propulsion import jet_models as jet

from marilib.airplane.propulsion.turbofan.turbofan_models import turbofan_thrust


#===========================================================================================================
def pte2_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    SFC for PTE2 architecture
    """

    #===========================================================================================================
    def fct_sfc(throttle,aircraft,pamb,tamb,mach,rating,thrust,nei):
        fn,sfc,sec,data = pte2_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)
        return (thrust - fn)

    x_ini = 0.8

    fct_arg = (aircraft,pamb,tamb,mach,rating,thrust,nei)

    output_dict = fsolve(fct_sfc, x0=x_ini, args=fct_arg, full_output=True)

    throttle = output_dict[0][0]

    fn,sfc,sec,data = pte2_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    return sfc


#===========================================================================================================
def pte2_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei):

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    battery = aircraft.pte2_battery
    power_elec = aircraft.pte2_power_elec_chain
    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    rear_shaft_power = {"MTO":r_engine.mto_r_shaft_power,
                        "MCN":r_engine.mcn_r_shaft_power,
                        "MCL":r_engine.mcl_r_shaft_power,
                        "MCR":r_engine.mcr_r_shaft_power,
                        "FID":r_engine.fid_r_shaft_power}

    # Battery power feed is used in temporary phases only
    power_factor = battery.power_feed * r_nacelle.controller_efficiency * r_nacelle.motor_efficiency
    battery_power_feed = {"MTO":power_factor,
                          "MCN":0.,
                          "MCL":power_factor,
                          "MCR":0.,
                          "FID":0.}

    fn,sfc0,data = turbofan_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei)
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    if (nacelle.rear_nacelle==1):

        power_offtake = throttle * rear_shaft_power[rating]/(nacelle.n_engine-nei) / power_elec.overall_efficiency

        # Shaft power of rear fan
        shaft_power2 = power_offtake * (nacelle.n_engine - nei)     # Shaft power dedicated to electric generator

        # Effective eFan shaft power
        pw_shaft2 =   shaft_power2 * power_elec.overall_efficiency \
                    + r_nacelle.motor_efficiency * r_nacelle.controller_efficiency * battery_power_feed[rating]

        if (propulsion.bli_effect>0):
            (fn_fan2,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,Pamb,Tamb,Mach,pw_shaft2)
            dVbli_o_V2 = dVbli/Vair
        else:
            (fn_fan2,q0) = jet.fan_thrust(r_nacelle,Pamb,Tamb,Mach,pw_shaft2)
            dVbli_o_V2 = 0.

        sec = (pw_shaft2/(r_nacelle.motor_efficiency*r_nacelle.controller_efficiency))/fn_fan2

    else:
        power_offtake = 0.
        sec = 0.
        fn_fan2 = 0.
        dVbli_o_V2 = 0.
        shaft_power2 = 0.

    # Shaft power of main fans
    shaft_power1 = shaft_power0 - power_offtake     # Shaft power dedicated to the fan

    if (propulsion.bli_effect>0):
        (fn_fan1,q1,dVbli) = jet.fan_thrust_with_bli(nacelle,Pamb,Tamb,Mach,shaft_power1)
        dVbli_o_V1 = dVbli/Vair
    else:
        (fn_fan1,q0) = jet.fan_thrust(nacelle,Pamb,Tamb,Mach,shaft_power1)
        dVbli_o_V1 = 0.

    sfc = sfc0 * (fn0*(nacelle.n_engine - nei)) / ((fn_core + fn_fan1)*(nacelle.n_engine - nei) + fn_fan2)

    fn = (fn_core + fn_fan1)*(nacelle.n_engine - nei) + fn_fan2

    data = (fn_core,fn_fan1,dVbli_o_V1,shaft_power1,fn_fan2,dVbli_o_V2,shaft_power2,fn0,shaft_power0,sfc0)

    return fn,sfc,sec,data


#===========================================================================================================
def blimp_body_drag(aircraft,Re,Mach):
    """
    Turbofan nacelle drag
    """

    wing = aircraft.wing
    blimp_body = aircraft.pte2_blimp_body

    fac = (1. + 0.126*Mach**2)

    blp_nwa = blimp_body.net_wetted_area

    blp_cxf =   1.05*((0.455/fac)*(numpy.log(10)/numpy.log(Re*blimp_body.length))**2.58)*blp_nwa/wing.area

    return blp_cxf,blp_nwa


#===========================================================================================================
def pte2_oei_drag(aircraft,nacelle,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    wing = aircraft.wing

    dCx = 0.01*nacelle.width**2 / wing.area     # Very low because the nacelle is burried into the boundary layer of the blimp

    return dCx

