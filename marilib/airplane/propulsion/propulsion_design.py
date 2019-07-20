#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy

from marilib.earth import environment as earth

from marilib.airplane.propulsion import propulsion_models as propu

from marilib.airplane.propulsion.turbofan.turbofan_design \
    import eval_turbofan_nacelle_design, eval_turbofan_engine_design, \
           eval_turbofan_pylon_mass, eval_turbofan_nacelle_mass

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

        engine = aircraft.turbofan_engine

        eval_turbofan_engine_design(aircraft)
        eval_turbofan_nacelle_design(aircraft)

    elif (propulsion.architecture=="PTE1"):

        engine = aircraft.turbofan_engine

        eval_turbofan_engine_design(aircraft)
        eval_pte1_engine_design(aircraft)
        eval_pte1_nacelle_design(aircraft)

    elif (propulsion.architecture=="EF1"):

        engine = aircraft.electrofan_engine

        eval_ef1_engine_design(aircraft)
        eval_ef1_nacelle_design(aircraft)

    else:
        raise Exception("propulsion.architecture index is out of range")


    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    disa = 15.
    altp = 0.
    mach = 0.25
    nei = 0.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,Data) = propu.thrust(aircraft,pamb,tamb,mach,MTO,nei)

    propulsion.reference_thrust_effective = (Fn/engine.n_engine)/0.80
    propulsion.mto_thrust_ref = Fn/engine.n_engine


    disa = aircraft.low_speed.disa_oei
    altp = aircraft.low_speed.req_oei_altp
    mach = 0.5*aircraft.design_driver.cruise_mach
    nei = 1.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCN,nei)

    propulsion.mcn_thrust_ref = Fn/(engine.n_engine-nei)


    disa = 0.
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    propulsion.sfc_cruise_ref = propu.sfc(aircraft,pamb,tamb,mach,MCR,nei)

    if (propulsion.architecture=="TF"):

        sec = 0.

    elif (propulsion.architecture=="PTE1"):

        fn,sec,data = propu.pte1_thrust(aircraft,pamb,tamb,mach,MCR,nei)

    elif (propulsion.architecture=="EF1"):

        fn,sec = propu.ef1_thrust(aircraft,pamb,tamb,mach,MCR,nei)

    else:
        raise Exception("propulsion.architecture index is out of range")

    propulsion.sec_cruise_ref = sec

    (Fn,Data) = propu.thrust(aircraft,pamb,tamb,mach,FID,nei)

    propulsion.fid_thrust_ref = Fn/engine.n_engine


    disa = 0.
    altp = aircraft.design_driver.top_of_climb_altp
    mach = aircraft.design_driver.cruise_mach
    nei = 0.

    (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

    (Fn,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCL,nei)

    propulsion.mcl_thrust_ref = Fn/engine.n_engine

    (Fn,Data) = propu.thrust(aircraft,pamb,tamb,mach,MCR,nei)

    propulsion.mcr_thrust_ref = Fn/engine.n_engine

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

    elif (propulsion.architecture=="PTE1"):

        pylon = aircraft.turbofan_pylon
        nacelle = aircraft.turbofan_nacelle

        e_nacelle = aircraft.rear_electric_nacelle
        power_elec = aircraft.pte1_power_elec_chain

        eval_turbofan_pylon_mass(aircraft)
        eval_pte1_nacelle_mass(aircraft)

        propulsion.mass = pylon.mass + nacelle.mass + e_nacelle.mass + power_elec.mass
        propulsion.c_g = (  pylon.c_g*pylon.mass + nacelle.c_g*nacelle.mass \
                          + e_nacelle.c_g*e_nacelle.mass + power_elec.c_g*power_elec.mass \
                          )/propulsion.mass

    elif (propulsion.architecture=="EF1"):

        pylon = aircraft.electrofan_pylon
        nacelle = aircraft.electrofan_nacelle

        eval_ef1_pylon_mass(aircraft)
        eval_ef1_nacelle_mass(aircraft)

        propulsion.mass = pylon.mass + nacelle.mass
        propulsion.c_g = (pylon.c_g*pylon.mass + nacelle.c_g*nacelle.mass)/propulsion.mass

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
       aircraft.weights.battery = 0.
    elif (aircraft.propulsion.architecture=="PTE1"):
       eval_pte1_battery_mass(aircraft)
    elif (aircraft.propulsion.architecture=="EF1"):
       eval_ef1_battery_mass(aircraft)

    return



