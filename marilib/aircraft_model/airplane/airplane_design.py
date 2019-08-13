#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

from marilib.earth import environment as earth

from marilib.aircraft_model.airplane import aerodynamics as airplane_aero


#===========================================================================================================
def eval_payload_mass(aircraft):
    """
    Payload mass & CG estimation
    """

    cabin = aircraft.cabin
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    payload = aircraft.payload

    payload.m_container_pallet = 4.36*fuselage.width*fuselage.length        # Containers and pallets

    if (aircraft.propulsion.fuel_type=="Battery"):
        if (aircraft.ef1_battery.stacking=="Variable"):
            payload.maximum = cabin.n_pax_ref*payload.m_pax_max
        else:
            payload.maximum = cabin.n_pax_ref*payload.m_pax_nominal # Because in this case, MZFW = MTOW
    else:
        payload.maximum = cabin.n_pax_ref*payload.m_pax_max

    payload.nominal = cabin.n_pax_ref*payload.m_pax_nominal

    payload.max_fwd_req_cg = cabin.fwd_limit + 0.35*cabin.length        # Payload max forward CG
    payload.max_fwd_mass = 0.60*cabin.n_pax_ref*payload.m_pax_max       # Payload mass for max forward CG

    payload.max_bwd_req_cg = cabin.fwd_limit + 0.70*cabin.length        # Payload max backward CG
    payload.max_bwd_mass = 0.70*cabin.n_pax_ref*payload.m_pax_max       # Payload mass for max backward CG

    payload.cg_container_pallet = 0.25*(cabin.fwd_limit + wing.x_root) + 0.25*(wing.x_root  +wing.c_root+cabin.fwd_limit+cabin.length)

    return


#===========================================================================================================
def eval_system_mass(aircraft):
    """
    Systems mass & CG estimation
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing
    ldg = aircraft.landing_gears
    vtp = aircraft.vertical_tail
    htp = aircraft.horizontal_tail
    weights = aircraft.weights

    systems = aircraft.systems

    propulsion = aircraft.propulsion

    systems.mass = 0.545*weights.mtow**0.8    # global mass of all systems

    systems.c_g = 0.50*fuselage.c_g + 0.20*wing.c_g + 0.10*ldg.c_g + 0.1*propulsion.c_g + 0.05*htp.c_g + 0.05*vtp.c_g

    return


#===========================================================================================================
def eval_aircraft_weights(aircraft):
    """
    Weights estimation
    """

    cabin = aircraft.cabin
    fuselage = aircraft.fuselage
    payload = aircraft.payload
    wing = aircraft.wing
    htp = aircraft.horizontal_tail
    vtp = aircraft.vertical_tail
    ldg = aircraft.landing_gears
    systems = aircraft.systems
    propulsion = aircraft.propulsion
    tanks = aircraft.tanks

    weights = aircraft.weights

    weights.mwe =  cabin.m_furnishing + fuselage.mass + wing.mass + htp.mass + vtp.mass \
                 + ldg.mass + systems.mass + propulsion.mass

    weights.owe = weights.mwe + cabin.m_op_item + payload.m_container_pallet + weights.battery_in_owe

    if (propulsion.fuel_type=="Battery"):
        mzfw = weights.mtow
    else:
        mzfw = weights.owe + payload.maximum

    weights.mass_constraint_1 = weights.mzfw - mzfw

    if (propulsion.fuel_type=="Battery"):
        mlw = weights.mtow
    else:
        if (cabin.n_pax_ref>100):
            mlw = min(weights.mtow , (1.07*weights.mzfw))
        else:
            mlw = weights.mtow

    weights.mass_constraint_2 = weights.mlw - mlw

    # WARNING : for EF1 architecture, MFW corresponds to max battery weight
    weights.mfw = min(tanks.mfw_volume_limited, weights.mtow - weights.owe)

    return


#===========================================================================================================
def eval_mass_coupling(aircraft):
    """
    Weights estimation internal coupling
    This relation is put apart from aircraft_weights because GEMS does not manage functions that compute their own input
    """

    cabin = aircraft.cabin
    payload = aircraft.payload
    weights = aircraft.weights

    if (aircraft.propulsion.fuel_type!="Battery"):

        weights.mzfw = weights.owe + payload.maximum

        if (cabin.n_pax_ref>100):
            weights.mlw = min(weights.mtow , (1.07*weights.mzfw))
        else:
            weights.mlw = weights.mtow

    return


#===========================================================================================================
def eval_aircraft_cg(aircraft):
    """
    Center of gravity estimation
    """

    cabin = aircraft.cabin
    fuselage = aircraft.fuselage
    payload = aircraft.payload
    wing = aircraft.wing
    htp = aircraft.horizontal_tail
    vtp = aircraft.vertical_tail
    ldg = aircraft.landing_gears
    systems = aircraft.systems
    propulsion = aircraft.propulsion
    tanks = aircraft.tanks
    weights = aircraft.weights

    c_g = aircraft.center_of_gravity

    c_g.mwe = (  fuselage.c_g*fuselage.mass + cabin.cg_furnishing*cabin.m_furnishing + wing.c_g*wing.mass \
               + ldg.c_g*ldg.mass + propulsion.c_g*propulsion.mass + htp.c_g*htp.mass \
               + vtp.c_g*vtp.mass + systems.c_g*systems.mass \
               )/weights.mwe

    c_g.owe = (  c_g.mwe*weights.mwe + cabin.cg_op_item*cabin.m_op_item + c_g.battery*weights.battery_in_owe \
               + payload.cg_container_pallet*payload.m_container_pallet \
               ) / weights.owe

    c_g.max_fwd_mass = weights.owe + tanks.fuel_max_fwd_mass + payload.max_fwd_mass
    c_g.max_fwd_req_cg =  (  c_g.owe*weights.owe + tanks.fuel_max_fwd_cg*tanks.fuel_max_fwd_mass \
                           + payload.max_fwd_req_cg*payload.max_fwd_mass \
                           )/c_g.max_fwd_mass

    c_g.max_bwd_mass = weights.owe + tanks.fuel_max_bwd_mass + payload.max_bwd_mass
    c_g.max_bwd_req_cg = (  c_g.owe*weights.owe + tanks.fuel_max_bwd_cg*tanks.fuel_max_bwd_mass \
                          + payload.max_bwd_req_cg*payload.max_bwd_mass \
                          )/c_g.max_bwd_mass

    return


#===========================================================================================================
def eval_aerodynamics_design(aircraft):
    """
    Defines high lift movable deflection settings
    HLDconf varies from 0 (clean) to 1 (full deflected)
    Typically : HLDconf = 1 ==> CzmaxLD
              : HLDconf = 0.1 to 0.5 ==> CzmaxTO
    """

    design_driver = aircraft.design_driver
    wing = aircraft.wing

    aerodynamics = aircraft.aerodynamics

    mach = design_driver.cruise_mach
    altp = design_driver.ref_cruise_altp
    disa = 0.

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)

    aerodynamics.cruise_lod_max, aerodynamics.cz_cruise_lod_max = airplane_aero.lod_max(aircraft, pamb, tamb, mach)

    aerodynamics.hld_conf_clean = 0.0   # By definition (0=<hld_conf=<1)
    aerodynamics.cz_max_clean,Cz0 = airplane_aero.high_lift(wing, aerodynamics.hld_conf_clean)

    aerodynamics.hld_conf_to = 0.3      # Take off (empirical setting)
    aerodynamics.cz_max_to,Cz0 = airplane_aero.high_lift(wing, aerodynamics.hld_conf_to)

    aerodynamics.hld_conf_ld = 1.0      # By definition (0=<hld_conf=<1), 1 is full landing
    aerodynamics.cz_max_ld,Cz0 = airplane_aero.high_lift(wing, aerodynamics.hld_conf_ld)

    return


