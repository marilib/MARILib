#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib import numpy

from marilib.earth import environment as earth

#===========================================================================================================
def eval_turboprop_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    low_speed = aircraft.low_speed

    engine = aircraft.turboprop_engine
    nacelle = aircraft.turboprop_nacelle
    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    engine.rating_factor = {"MTO":1.00, "MCN":0.90, "MCL":0.90, "MCR":0.80, "FID":0.10}

    engine.reference_thrust = propulsion.reference_thrust

    # Propulsion architecture design, definition reference conditions for engine performances
    #-----------------------------------------------------------------------------------------------------------

    # Initialisation
    crm = design_driver.cruise_mach
    toc = design_driver.top_of_climb_altp
    rca = design_driver.ref_cruise_altp
    roa = low_speed.req_oei_altp

    #                      MTO   MCN    MCL  MCR  FID
    fd_disa = {"MTO":15. , "MCN":15.  , "MCL":0. , "MCR":0. , "FID":0. }
    fd_altp = {"MTO":0.  , "MCN":roa  , "MCL":toc, "MCR":rca, "FID":rca}
    fd_mach = {"MTO":0.15, "MCN":crm/2, "MCL":crm, "MCR":crm, "FID":crm}
    fd_nei  = {"MTO":0.  , "MCN":1.   , "MCL":0. , "MCR":0. , "FID":0. }

    propulsion.flight_data = {"disa":fd_disa, "altp":fd_altp, "mach":fd_mach, "nei":fd_nei}

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

    return


#===========================================================================================================
def eval_turboprop_nacelle_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing
    vtp = aircraft.vertical_tail
    engine = aircraft.turboprop_engine
    propulsion = aircraft.propulsion

    nacelle = aircraft.turboprop_nacelle

    nacelle.n_engine = propulsion.n_engine

    nacelle.width = 0.25*(engine.reference_power/1.e3)**0.2        # statistical regression

    nacelle.length = 0.84*(engine.reference_power/1.e3)**0.2       # statistical regression

    nacelle.net_wetted_area = (2.3*(engine.reference_power/1.e3)**0.2)*nacelle.n_engine     # statistical regression

    nacelle.hub_width = 0.2

    nacelle.propeller_width = numpy.sqrt((4./numpy.pi)*(engine.reference_thrust/3000.))      # Assuming 3000 N/m2

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1):

        if (nacelle.n_engine==2):

            nacelle.y_ext = 0.5 * fuselage.width + 0.8 * nacelle.propeller_width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = wing.z_root + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

        elif (nacelle.n_engine==4):

            nacelle.y_int = 0.5 * fuselage.width + 0.8 * nacelle.propeller_width      # statistical regression

            nacelle.x_int = wing.x_root + (nacelle.y_int-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_int = wing.z_root + (nacelle.y_int - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

            nacelle.y_ext = 0.5 * fuselage.width + 2.1 * nacelle.propeller_width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = wing.z_root + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

        else:
            raise Exception("nacelle.n_engine, number of engine not supported")

    else:
        raise Exception("nacelle.attachment, index is out of range")

    propulsion.y_ext_nacelle = nacelle.y_ext

    nacelle.rear_nacelle = 0    # No rear nacelle by default in this architecture

    return


#===========================================================================================================
def eval_turboprop_nacelle_mass(aircraft):
    """
    Thermal propulsive nacelle mass estimation
    """

    engine = aircraft.turboprop_engine
    nacelle = aircraft.turboprop_nacelle

    nacelle.mass = (1.266*(engine.reference_power/1.e3)**0.9)*nacelle.n_engine       # statistical regression

    nacelle.c_g = nacelle.x_ext + 0.7 * nacelle.length      # statistical regression

    return


