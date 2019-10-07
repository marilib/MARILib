#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

from marilib import numpy

from marilib.earth import environment as earth

from marilib.airplane.propulsion import propulsion_models as propu

from marilib.airplane.propulsion.turbofan.turbofan_design \
    import eval_turbofan_nacelle_design, eval_turbofan_engine_design, \
           eval_turbofan_pylon_mass, eval_turbofan_nacelle_mass

from marilib.airplane.propulsion.turboprop.turboprop_design \
    import eval_turboprop_nacelle_design, eval_turboprop_engine_design, \
           eval_turboprop_nacelle_mass

from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_design \
    import eval_pte1_nacelle_design, eval_pte1_engine_design, \
           eval_pte1_nacelle_mass, eval_pte1_battery_mass

from marilib.airplane.propulsion.electric_ef1.electric_ef1_design \
    import eval_ef1_nacelle_design, eval_ef1_engine_design, eval_battery_cg_range, \
           eval_ef1_pylon_mass, eval_ef1_nacelle_mass, eval_ef1_battery_mass

from marilib.airplane.airframe.airframe_design import eval_wing_tank_data, eval_fuel_cg_range


#===========================================================================================================
def eval_propulsion_design(aircraft):
    """
    Propulsion architecture design
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):
        eval_turbofan_engine_design(aircraft)
        eval_turbofan_nacelle_design(aircraft)
    elif (propulsion.architecture=="TP"):
        eval_turboprop_engine_design(aircraft)
        eval_turboprop_nacelle_design(aircraft)
    elif (propulsion.architecture=="PTE1"):
        eval_turbofan_engine_design(aircraft)
        eval_pte1_engine_design(aircraft)
        eval_pte1_nacelle_design(aircraft)
    elif (propulsion.architecture=="EF1"):
        eval_ef1_engine_design(aircraft)
        eval_ef1_nacelle_design(aircraft)
    else:
        raise Exception("propulsion.architecture index is out of range")


    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MTO]
    altp = propulsion.flight_data["altp"][MTO]
    mach = propulsion.flight_data["mach"][MTO]
    nei = propulsion.flight_data["nei"][MTO]

    throttle = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    Fn,SFC,SEC,Data = propu.thrust(aircraft,pamb,tamb,mach,MTO,throttle,nei)

    propulsion.mto_thrust_ref = Fn/propulsion.n_engine

    propulsion.reference_thrust_effective = (Fn/propulsion.n_engine)/0.80

    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MCN]
    altp = propulsion.flight_data["altp"][MCN]
    mach = propulsion.flight_data["mach"][MCN]
    nei = propulsion.flight_data["nei"][MCN]

    throttle = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,SFC,SEC,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCN,throttle,nei)

    propulsion.mcn_thrust_ref = Fn/(propulsion.n_engine-nei)

    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MCL]
    altp = propulsion.flight_data["altp"][MCL]
    mach = propulsion.flight_data["mach"][MCL]
    nei = propulsion.flight_data["nei"][MCL]

    throttle = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,SFC,SEC,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCL,throttle,nei)

    propulsion.mcl_thrust_ref = Fn/propulsion.n_engine

    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MCR]
    altp = propulsion.flight_data["altp"][MCR]
    mach = propulsion.flight_data["mach"][MCR]
    nei = propulsion.flight_data["nei"][MCR]

    throttle = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    # WARNING : SFC & SEC cruise reference corresponds to MCR thrust, actual values may be lower
    if (propulsion.architecture=="TF"):
        fn,sfc,data = propu.turbofan_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
        propulsion.sfc_cruise_ref = sfc
        propulsion.sec_cruise_ref = 0.
    elif (propulsion.architecture=="TP"):
        fn,sfc,data = propu.turboprop_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
        propulsion.sfc_cruise_ref = sfc
        propulsion.sec_cruise_ref = 0.
    elif (propulsion.architecture=="PTE1"):
        fn,sfc,sec,data = propu.pte1_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
        propulsion.sfc_cruise_ref = sfc
        propulsion.sec_cruise_ref = sec
    elif (propulsion.architecture=="EF1"):
        fn,sec,data = propu.ef1_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
        propulsion.sfc_cruise_ref = 0.
        propulsion.sec_cruise_ref = sec
    else:
        raise Exception("propulsion.architecture index is out of range")

    (Fn,SFC,SEC,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)

    propulsion.mcr_thrust_ref = Fn/propulsion.n_engine

    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][FID]
    altp = propulsion.flight_data["altp"][FID]
    mach = propulsion.flight_data["mach"][FID]
    nei = propulsion.flight_data["nei"][FID]

    throttle = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,SFC,SEC,Data) = propu.thrust(aircraft,pamb,tamb,mach,FID,throttle,nei)

    propulsion.fid_thrust_ref = Fn/propulsion.n_engine

    return


#===========================================================================================================
def eval_propulsion_mass(aircraft):
    """
    Propulsion mass & CG estimation
    """

    propulsion = aircraft.propulsion

    if (propulsion.architecture=="TF"):

        pylon = aircraft.turbofan_pylon
        nacelle = aircraft.turbofan_nacelle

        eval_turbofan_pylon_mass(aircraft)
        eval_turbofan_nacelle_mass(aircraft)

        propulsion.mass = pylon.mass + nacelle.mass
        propulsion.c_g = (pylon.c_g*pylon.mass + nacelle.c_g*nacelle.mass)/propulsion.mass

    elif (propulsion.architecture=="TP"):

        nacelle = aircraft.turboprop_nacelle

        eval_turboprop_nacelle_mass(aircraft)

        propulsion.mass = nacelle.mass
        propulsion.c_g = nacelle.c_g

    elif (propulsion.architecture=="PTE1"):

        pylon = aircraft.turbofan_pylon
        nacelle = aircraft.turbofan_nacelle
        r_nacelle = aircraft.rear_electric_nacelle
        power_elec = aircraft.pte1_power_elec_chain

        eval_turbofan_pylon_mass(aircraft)
        eval_pte1_nacelle_mass(aircraft)

        propulsion.mass = pylon.mass + nacelle.mass + r_nacelle.mass + power_elec.mass
        propulsion.c_g = (  pylon.c_g*pylon.mass + nacelle.c_g*nacelle.mass \
                          + r_nacelle.c_g*r_nacelle.mass + power_elec.c_g*power_elec.mass \
                          )/propulsion.mass

    elif (propulsion.architecture=="EF1"):

        pylon = aircraft.electrofan_pylon
        nacelle = aircraft.electrofan_nacelle
        r_nacelle = aircraft.rear_electric_nacelle
        power_elec = aircraft.ef1_power_elec_chain

        eval_ef1_pylon_mass(aircraft)
        eval_ef1_nacelle_mass(aircraft)

        propulsion.mass = pylon.mass + nacelle.mass + r_nacelle.mass + power_elec.mass
        propulsion.c_g = (  pylon.c_g*pylon.mass + nacelle.c_g*nacelle.mass \
                          + r_nacelle.c_g*r_nacelle.mass + power_elec.c_g*power_elec.mass \
                          )/propulsion.mass

    else:
        raise Exception("propulsion.architecture index is out of range")

    return


#===========================================================================================================
def eval_tank_mass(aircraft):
    """
    Tank predesign
    """

    propulsion = aircraft.propulsion

    eval_wing_tank_data(aircraft)

    if (propulsion.fuel_type=="Kerosene"):
        eval_fuel_cg_range(aircraft)
    elif (propulsion.fuel_type=="Battery"):
        eval_battery_cg_range(aircraft)
    else:
        raise Exception("propulsion.fuel_type is not allowed")


#===========================================================================================================
def eval_battery_mass(aircraft):
    """
    Battery mass and CG estimation
    """

    if (aircraft.propulsion.architecture=="TF"):
       aircraft.propulsion.battery_energy_density = 0.
       aircraft.center_of_gravity.battery = 0.
       aircraft.weights.battery_in_owe = 0.
    elif (aircraft.propulsion.architecture=="TP"):
       aircraft.propulsion.battery_energy_density = 0.
       aircraft.center_of_gravity.battery = 0.
       aircraft.weights.battery_in_owe = 0.
    elif (aircraft.propulsion.architecture=="PTE1"):
        eval_pte1_battery_mass(aircraft)
    elif (aircraft.propulsion.architecture=="EF1"):
        eval_ef1_battery_mass(aircraft)

    return



