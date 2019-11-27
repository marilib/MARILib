#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib import numpy
from marilib import fsolve

from marilib.earth import environment as earth


#===========================================================================================================
def turboprop_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    PSFC for a turboprop
    """

    #===========================================================================================================
    def fct_sfc(throttle,aircraft,pamb,tamb,mach,rating,thrust,nei):
        fn,sfc,data = turboprop_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)
        return (thrust - fn)

    x_ini = 0.8

    fct_arg = (aircraft,pamb,tamb,mach,rating,thrust,nei)

    output_dict = fsolve(fct_sfc, x0=x_ini, args=fct_arg, full_output=True)

    throttle = output_dict[0][0]

    fn,sfc,data = turboprop_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    return sfc


#===========================================================================================================
def turboprop_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei):
    """
    Calculation of thrust for pure turboprop airplane
    Warning : ALL engine thrust returned
    """

    propulsion = aircraft.propulsion
    engine = aircraft.turboprop_engine
    nacelle = aircraft.turboprop_nacelle

    factor = engine.rating_factor       # [MTO,MCN,MCL,MCR,FID]

    eta_prop = nacelle.efficiency_prop

    psfc_ref = 0.4*1.68969e-07   # 0.4 lb/shp/h

    psfc = psfc_ref * earth.fuel_heat("Kerosene") / propulsion.fuel_heat

    rho,sig = earth.air_density(Pamb,Tamb)
    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    shaft_power0 = throttle*factor[rating]*engine.reference_power*sig**0.5

    fn0 = eta_prop*shaft_power0/Vair

    fn = fn0*(nacelle.n_engine - nei)        # All turbofan thrust

    sfc = psfc*shaft_power0/fn0

    data = (fn0,shaft_power0)   # Data for ONE turbofan engine

    return fn,sfc,data


#===========================================================================================================
def turboprop_oei_drag(aircraft,nacelle,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    wing = aircraft.wing

    dCx = 0.12*nacelle.width**2 / wing.area

    return dCx
