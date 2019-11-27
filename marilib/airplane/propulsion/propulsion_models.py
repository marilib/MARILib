#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

from marilib import numpy

from marilib import fsolve
from marilib.tools.math import lin_interp_1d

from marilib.earth import environment as earth

from marilib.airplane.propulsion.turbofan.turbofan_models \
    import turbofan_sfc, turbofan_thrust, turbofan_oei_drag

from marilib.airplane.propulsion.turboprop.turboprop_models \
    import turboprop_sfc, turboprop_thrust, turboprop_oei_drag

from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_models \
    import pte1_sfc, pte1_thrust

from marilib.airplane.propulsion.electric_ef1.electric_ef1_models \
    import ef1_sec, ef1_thrust, electrofan_oei_drag

from marilib.airplane.propulsion import jet_models as jet

#===========================================================================================================
def sfc(aircraft,pamb,tamb,mach,rating,thrust,nei):
    """
    Bucket SFC for a turbofan
    IMPORTANT REMARK : for EF1 architecture, SFC is in Energy unit per Thrust unit (J/N)
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        sfc = turbofan_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei)
    elif (propulsion.architecture=="TP"):
        sfc = turboprop_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei)
    elif (propulsion.architecture=="PTE1"):
        sfc = pte1_sfc(aircraft,pamb,tamb,mach,rating,thrust,nei)
    elif (propulsion.architecture=="EF1"):
        sfc = ef1_sec(aircraft,pamb,tamb,mach,rating,thrust,nei)
    else:
        raise Exception("propulsion.architecture index is out of range")

    return sfc


#===========================================================================================================
def thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei):
    """
    Calculation of the total thrust of the architecture
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        fn,sfc,data = turbofan_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei)
        sec = 0.
    elif (propulsion.architecture=="TP"):
        fn,sfc,data = turboprop_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei)
        sec = 0.
    elif (propulsion.architecture=="PTE1"):
        fn,sfc,sec,data = pte1_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei)
    elif (propulsion.architecture=="EF1"):
        sfc = 0.
        fn,sec,data = ef1_thrust(aircraft,Pamb,Tamb,Mach,rating,throttle,nei)
    else:
        raise Exception("propulsion.architecture index is out of range")

    return fn,sfc,sec,data


#===========================================================================================================
def nacelle_drag(aircraft,Re,Mach):
    """
    All nacelle drag coefficient
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        nacelle = aircraft.turbofan_nacelle
        nacelle_cxf,nacelle_nwa = jet.nacelle_generic_drag(aircraft,nacelle,Re,Mach)
    elif (propulsion.architecture=="TP"):
        nacelle = aircraft.turboprop_nacelle
        nacelle_cxf,nacelle_nwa = jet.nacelle_generic_drag(aircraft,nacelle,Re,Mach)
    elif (propulsion.architecture=="PTE1"):
        nacelle = aircraft.turbofan_nacelle
        t_nacelle_cxf,t_nacelle_nwa = jet.nacelle_generic_drag(aircraft,nacelle,Re,Mach)
        r_nacelle = aircraft.rear_electric_nacelle
        r_nacelle_cxf,r_nacelle_nwa = jet.nacelle_generic_drag(aircraft,r_nacelle,Re,Mach)
        nacelle_nwa = t_nacelle_nwa + r_nacelle_nwa
        nacelle_cxf = (t_nacelle_cxf*t_nacelle_nwa + r_nacelle_cxf*r_nacelle_nwa)/nacelle_nwa
    elif (propulsion.architecture=="EF1"):
        e_nacelle = aircraft.electrofan_nacelle
        e_nacelle_cxf,e_nacelle_nwa = jet.nacelle_generic_drag(aircraft,e_nacelle,Re,Mach)
        if (e_nacelle.rear_nacelle==1):
            r_nacelle = aircraft.rear_electric_nacelle
            r_nacelle_cxf,r_nacelle_nwa = jet.nacelle_generic_drag(aircraft,r_nacelle,Re,Mach)
        else:
            r_nacelle_cxf = 0.
            r_nacelle_nwa = 0.
        nacelle_nwa = e_nacelle_nwa + r_nacelle_nwa
        nacelle_cxf = (e_nacelle_cxf*e_nacelle_nwa + r_nacelle_cxf*r_nacelle_nwa)/nacelle_nwa
    else:
        raise Exception("propulsion.architecture index is out of range")

    return nacelle_cxf,nacelle_nwa


#===========================================================================================================
def oei_drag(aircraft,pamb,tamb):
    """
    Inoperative engine drag coefficient
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        nacelle = aircraft.turbofan_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)
    elif (propulsion.architecture=="TP"):
        nacelle = aircraft.turboprop_nacelle
        dcx = turboprop_oei_drag(aircraft,nacelle,pamb,tamb)
    elif (propulsion.architecture=="PTE1"):
        nacelle = aircraft.turbofan_nacelle
        dcx = turbofan_oei_drag(aircraft,nacelle,pamb,tamb)
    elif (propulsion.architecture=="EF1"):
        nacelle = aircraft.electrofan_nacelle
        dcx = electrofan_oei_drag(aircraft,nacelle,pamb,tamb)
    else:
        raise Exception("propulsion.architecture index is out of range")

    return dcx


#===========================================================================================================
def thrust_pitch_moment(aircraft,fn,pamb,mach,dcx_oei):

    propulsion = aircraft.propulsion
    wing = aircraft.wing

    gam = earth.heat_ratio()

    if (propulsion.architecture=="TF"):
        nacelle = aircraft.turbofan_nacelle
    elif (propulsion.architecture=="TP"):
        nacelle = aircraft.turboprop_nacelle
    elif (propulsion.architecture=="PTE1"):
        nacelle = aircraft.turbofan_nacelle
    elif (propulsion.architecture=="EF1"):
        nacelle = aircraft.electrofan_nacelle
    else:
        raise Exception("propulsion.architecture index is out of range")

    if (nacelle.n_engine==2):
        cm_prop = nacelle.z_ext*(dcx_oei - fn/(0.5*gam*pamb*mach**2*wing.area))
    elif (nacelle.n_engine==4):
        cm_prop =   nacelle.z_ext*(dcx_oei - (fn/3.)/(0.5*gam*pamb*mach**2*wing.area)) \
                  - nacelle.z_int*(2.*fn/3.)/(0.5*gam*pamb*mach**2*wing.area)
    else:
        raise Exception("thrust_pitch_moment, Number of engine is not supported")

    return cm_prop


#===========================================================================================================
def tail_cone_drag_effect(aircraft):
    """
    Effect of rear fan installation on fuselage tail cone drag
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        fac = 1.
    elif (propulsion.architecture=="TP"):
        fac = 1.
    elif (propulsion.architecture=="PTE1"):
        dfus = aircraft.fuselage.width
        dhub = aircraft.rear_electric_nacelle.hub_width
        fac = (1.-(dhub/dfus)**2)
    elif (propulsion.architecture=="EF1"):
        if (aircraft.electrofan_nacelle.rear_nacelle):
            dfus = aircraft.fuselage.width
            dhub = aircraft.rear_electric_nacelle.hub_width
            fac = (1.-(dhub/dfus)**2)
        else:
            fac = 1.
    else:
        raise Exception("propulsion.architecture index is out of range")

    return fac






