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


#===========================================================================================================
def eval_ef1_pylon_mass(aircraft):
    """
    Electrofan pylon mass & CG estimation
    """

    nacelle = aircraft.electrofan_nacelle
    engine = aircraft.electrofan_engine

    pylon = aircraft.electrofan_pylon

    pylon.mass = 0.0031*engine.reference_thrust*engine.n_engine

    if (engine.n_engine==2):
        pylon.c_g = nacelle.x_ext + 0.75*nacelle.length
    elif (engine.n_engine==4):
        pylon.c_g = 0.5*(nacelle.x_int + nacelle.x_ext) + 0.75*nacelle.length
    else:
        raise Exception("Number of engine is not allowed")

    return


#===========================================================================================================
def eval_ef1_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    low_speed = aircraft.low_speed
    propulsion = aircraft.propulsion

    power_elec = aircraft.ef1_power_elec_chain
    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    engine.rating_factor = {"MTO":1.00, "MCN":0.80, "MCL":0.80, "MCR":0.80, "FID":0.05}

    # Propulsion architecture design, definition reference conditions for engine performances
    #-----------------------------------------------------------------------------------------------------------

    # Initialisation
    crm = design_driver.cruise_mach
    toc = design_driver.top_of_climb_altp
    rca = design_driver.ref_cruise_altp
    roa = low_speed.req_oei_altp

    #                      MTO   MCN    MCL  MCR  FID
    fd_disa = {"MTO":15. , "MCN":0.   , "MCL":0. , "MCR":0. , "FID":0. }
    fd_altp = {"MTO":0.  , "MCN":roa  , "MCL":toc, "MCR":rca, "FID":rca}
    fd_mach = {"MTO":0.25, "MCN":crm/2, "MCL":crm, "MCR":crm, "FID":crm}
    fd_nei  = {"MTO":0.  , "MCN":1.   , "MCL":0. , "MCR":0. , "FID":0. }

    propulsion.flight_data = {"disa":fd_disa, "altp":fd_altp, "mach":fd_mach, "nei":fd_nei}

    # Propulsion architecture design, definition of e-fan power in each flight phase
    #-----------------------------------------------------------------------------------------------------------
    # Main engines
    #-----------------------------------------------------------------------------------------------------------
    disa = fd_disa[MTO]
    altp = fd_altp[MTO]
    mach = fd_mach[MTO]

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    Vsnd = earth.sound_speed(tamb)
    Vair = mach*Vsnd

    # reference_thrust is driving the design through ref_shat_power but actual thrust is computed later
    ref_shaft_power = 0.8 * engine.reference_thrust * Vair / nacelle.efficiency_prop

    engine.reference_power = ref_shaft_power

    engine.mto_e_shaft_power = ref_shaft_power * engine.rating_factor[MTO]
    engine.mcn_e_shaft_power = ref_shaft_power * engine.rating_factor[MCN]
    engine.mcl_e_shaft_power = ref_shaft_power * engine.rating_factor[MCL]
    engine.mcr_e_shaft_power = ref_shaft_power * engine.rating_factor[MCR]
    engine.fid_e_shaft_power = ref_shaft_power * engine.rating_factor[FID]

    # Max rear fan shaft power, if any
    #-----------------------------------------------------------------------------------------------------------
    r_shaft_power = numpy.array([engine.mto_r_shaft_power,
                                 engine.mcn_r_shaft_power,
                                 engine.mcl_r_shaft_power,
                                 engine.mcr_r_shaft_power,
                                 engine.fid_r_shaft_power])

    power_elec.max_power = max(r_shaft_power)
    power_elec.max_power_rating = numpy.argmax(r_shaft_power)

    return


#===========================================================================================================
def eval_ef1_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing
    vtp = aircraft.vertical_tail

    propulsion = aircraft.propulsion

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    disa = 0.
    altp = design_driver.ref_cruise_altp
    mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)

    shaft_power = engine.mcr_e_shaft_power

    hub_width = 0.2

    jet.efan_nacelle_design(nacelle,Pamb,Tamb,mach,shaft_power,hub_width)

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1):

        if (engine.n_engine==2):

            nacelle.y_ext = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = - 0.5 * fuselage.height \
                            + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width

        elif (engine.n_engine==4):

            nacelle.y_int = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_int = wing.x_root + (nacelle.y_int-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_int = - 0.5 * fuselage.height \
                            + (nacelle.y_int - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width

            nacelle.y_ext = 2.0 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = - 0.5 * fuselage.height \
                            + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width
        else:
            raise Exception("engine.n_engine, number of engine not supported")

    elif (nacelle.attachment == 2):

        if (engine.n_engine==2):

            nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

            nacelle.x_ext = vtp.x_root - 0.5*nacelle.length

            nacelle.z_ext = 0.5 * fuselage.height

        else:
            raise Exception("engine.n_engine, number of engine not supported")

    else:
        raise Exception("nacelle.attachment, index is out of range")

    return


#===========================================================================================================
def eval_ef1_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    power_elec = aircraft.ef1_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    shaftPowerMax = power_elec.max_power

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                      ) * shaftPowerMax

    nacelle.mass = (  1./nacelle.controller_pw_density + 1./nacelle.motor_pw_density \
                    + 1./nacelle.nacelle_pw_density \
                   ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    nacelle.c_g = fuselage.length + 0.5*nacelle.length

    return


#===========================================================================================================
def eval_battery_cg_range(aircraft):
    """
    Wing battery predesign using tank data structure
    """

    tanks = aircraft.tanks

    # Need to take account of any possible battery loading
    # WARNING : in early design steps, it may occur that the resulting weight of the airplane would be higher than MTOW
    tanks.fuel_max_fwd_cg = tanks.fuel_central_cg    # Battery max forward CG, central volume is forward only within backward swept wing
    tanks.fuel_max_fwd_mass = tanks.central_volume*tanks.fuel_density

    tanks.fuel_max_bwd_cg = tanks.fuel_cantilever_cg    # Battery max Backward CG
    tanks.fuel_max_bwd_mass = tanks.cantilever_volume*tanks.fuel_density

    return


#===========================================================================================================
def eval_ef1_battery_mass(aircraft):
    """
    Nevertheless, for simplicity reason, battery CG is supposed constant
    """

    tanks = aircraft.tanks
    propulsion = aircraft.propulsion
    battery = aircraft.ef1_battery

    propulsion.battery_energy_density = battery.energy_density

    battery.mass_max = tanks.mfw_volume_limited
    battery.c_g = tanks.fuel_total_cg

    aircraft.weights.battery = 0.
    aircraft.center_of_gravity.battery = battery.c_g

    return


