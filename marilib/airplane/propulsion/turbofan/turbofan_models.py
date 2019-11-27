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


#===========================================================================================================
def turbofan_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    SFC for a turbofan
    """

    #===========================================================================================================
    def fct_sfc(throttle,aircraft,pamb,tamb,mach,rating,thrust,nei):
        fn,sfc,data = turbofan_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)
        return (thrust - fn)

    x_ini = 0.8

    fct_arg = (aircraft,pamb,tamb,mach,rating,thrust,nei)

    output_dict = fsolve(fct_sfc, x0=x_ini, args=fct_arg, full_output=True)

    throttle = output_dict[0][0]

    fn,sfc,data = turbofan_thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    return sfc


#===========================================================================================================
def turbofan_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei):
    """
    Calculation of thrust for pure turbofan airplane
    Warning : ALL engine thrust returned
    """

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    sfc_ref = ( 0.4 + 1./engine.bpr**0.895 )/36000.

    sfc = sfc_ref * earth.fuel_heat("Kerosene") / propulsion.fuel_heat

    factor = engine.rating_factor       # [MTO,MCN,MCL,MCR,FID]

    kth =  0.475*Mach**2 + 0.091*(engine.bpr/10.)**2 \
         - 0.283*Mach*engine.bpr/10. \
         - 0.633*Mach - 0.081*engine.bpr/10. + 1.192

    (rho,sig) = earth.air_density(Pamb,Tamb)

    fn0 = throttle*factor[rating]*kth*engine.reference_thrust*sig**0.75

    fn_core = fn0 * engine.core_thrust_ratio        # Core thrust

    fn_fan0 = fn0 * (1.-engine.core_thrust_ratio)    # Fan thrust

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    shaft_power0 = fn_fan0*Vair/nacelle.efficiency_prop   # Available total shaft power for one engine

    fn = fn0*(propulsion.n_engine - nei)        # All turbofan thrust

    data = (fn_core,fn_fan0,fn0,shaft_power0)   # Data for ONE turbofan engine

    return fn,sfc,data


#===========================================================================================================
def turbofan_oei_drag(aircraft,nacelle,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    wing = aircraft.wing

    dCx = 0.12*nacelle.width**2 / wing.area

    return dCx



