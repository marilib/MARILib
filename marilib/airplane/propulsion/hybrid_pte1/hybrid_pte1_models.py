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
def pte1_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    SFC for PTE1 architecture
    """

    #===========================================================================================================
    def fct_sfc(throttle,aircraft,pamb,tamb,mach,rating,thrust,nei):
        fn,sfc,sec,data = pte1_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)
        return (thrust - fn)

    x_ini = 0.8

    fct_arg = (aircraft,pamb,tamb,mach,rating,thrust,nei)

    output_dict = fsolve(fct_sfc, x0=x_ini, args=fct_arg, full_output=True)

    throttle = output_dict[0][0]

    fn,sfc,sec,data = pte1_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    return sfc


#===========================================================================================================
def pte1_sfc_old(aircraft,pamb,tamb,mach,rating,nei):
    """
    Bucket SFC for a turbofan
    """

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    power_elec = aircraft.pte1_power_elec_chain
    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    power_ratio = {"MTO":power_elec.mto_e_power_ratio,
                   "MCN":power_elec.mcn_e_power_ratio,
                   "MCL":power_elec.mcl_e_power_ratio,
                   "MCR":power_elec.mcr_e_power_ratio,
                   "FID":power_elec.fid_e_power_ratio}

    sfc0_ref = ( 0.4 + 1./engine.bpr**0.895 )/36000.

    sfc0 = sfc0_ref * earth.fuel_heat("Kerosene") / propulsion.fuel_heat

    if (propulsion.bli_effect>0):
        kBLIe = propulsion.bli_r_thrust_factor
    else:
        kBLIe = 1.

    kC = engine.core_thrust_ratio
    kW = power_ratio[rating]

    eff_prop = nacelle.efficiency_prop
    eff_e_prop = r_nacelle.efficiency_prop
    eff_chain = power_elec.overall_efficiency

    eff_h = kC + (1.-kC)*( kW*kBLIe*(eff_e_prop/eff_prop)*eff_chain + (1.-kW) )

    sfc = sfc0 / eff_h

    return sfc


#===========================================================================================================
def pte1_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei):

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    battery = aircraft.pte1_battery
    power_elec = aircraft.pte1_power_elec_chain
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

    power_offtake = throttle * rear_shaft_power[rating]/(nacelle.n_engine-nei) / power_elec.overall_efficiency

    shaft_power1 = shaft_power0 - power_offtake     # Shaft power dedicated to the fan

    fn_fan1 = nacelle.efficiency_prop*shaft_power1/Vair     # Effective fan thrust

    shaft_power2 = power_offtake * (nacelle.n_engine - nei)     # Shaft power dedicated to electric generator

    # Effective eFan shaft power
    pw_shaft2 =   shaft_power2 * power_elec.overall_efficiency \
                + r_nacelle.motor_efficiency * r_nacelle.controller_efficiency * battery_power_feed[rating]

    if (pw_shaft2 > 0.):

        if (propulsion.bli_effect>0):
            (fn_fan2,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,Pamb,Tamb,Mach,pw_shaft2)
            dVbli_o_V = dVbli/Vair
        else:
            (fn_fan2,q0) = jet.fan_thrust(r_nacelle,Pamb,Tamb,Mach,pw_shaft2)
            dVbli_o_V = 0.

        sec = (pw_shaft2/(r_nacelle.motor_efficiency*r_nacelle.controller_efficiency))/fn_fan2

    else:

        dVbli_o_V  = 0.
        fn_fan2 = 0.
        sec = 0.

    sfc = sfc0 * (fn0*(nacelle.n_engine - nei)) / ((fn_core + fn_fan1)*(nacelle.n_engine - nei) + fn_fan2)

    fn = (fn_core + fn_fan1)*(nacelle.n_engine - nei) + fn_fan2

    data = (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0,sfc0)

    return fn,sfc,sec,data


