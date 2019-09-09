#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

from marilib import numpy


#===========================================================================================================
def eval_turbofan_pylon_mass(aircraft):
    """
    Turbofan pylon mass & CG estimation
    """

    nacelle = aircraft.turbofan_nacelle
    engine = aircraft.turbofan_engine

    pylon = aircraft.turbofan_pylon

    pylon.mass = 0.0031*engine.reference_thrust*nacelle.n_engine

    if (nacelle.n_engine==2):
        pylon.c_g = nacelle.x_ext + 0.75*nacelle.length
    elif (nacelle.n_engine==4):
        pylon.c_g = 0.5*(nacelle.x_int + nacelle.x_ext) + 0.75*nacelle.length
    else:
        raise Exception("Number of engine is not allowed")

    return


#===========================================================================================================
def eval_turbofan_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    low_speed = aircraft.low_speed

    engine = aircraft.turbofan_engine
    propulsion = aircraft.propulsion

#    engine.rating_factor = {"MTO":0.800, "MCN":0.688, "MCL":0.624, "MCR":0.560, "FID":0.100}

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

    return


#===========================================================================================================
def eval_turbofan_nacelle_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    fuselage = aircraft.fuselage
    wing = aircraft.wing
    vtp = aircraft.vertical_tail
    engine = aircraft.turbofan_engine
    propulsion = aircraft.propulsion

    nacelle = aircraft.turbofan_nacelle

    nacelle.width = 0.5 * engine.bpr ** 0.7 + 5.E-6 * engine.reference_thrust

    nacelle.length = 0.86 * nacelle.width + engine.bpr ** 0.37      # statistical regression

    Knac = numpy.pi * nacelle.width * nacelle.length

    nacelle.net_wetted_area = Knac*(1.48 - 0.0076*Knac)*nacelle.n_engine        # statistical regression

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1):

        if (nacelle.n_engine==2):

            nacelle.y_ext = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

        elif (nacelle.n_engine==4):

            nacelle.y_int = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_int = wing.x_root + (nacelle.y_int-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_int = (nacelle.y_int - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

            nacelle.y_ext = 2.0 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) - 0.5*nacelle.width

        else:
            raise Exception("nacelle.n_engine, number of engine not supported")

    elif (nacelle.attachment == 2):

        if (nacelle.n_engine==2):

            nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

            nacelle.x_ext = vtp.x_root - 0.5*nacelle.length

            nacelle.z_ext = fuselage.height

        else:
            raise Exception("nacelle.n_engine, number of engine not supported")

    else:
        raise Exception("nacelle.attachment, index is out of range")

    nacelle.rear_nacelle = 0    # No rear nacelle in this architecture

    return


#===========================================================================================================
def eval_turbofan_nacelle_mass(aircraft):
    """
    Thermal propulsive nacelle mass estimation
    """

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    nacelle.mass = (1250. + 0.021*engine.reference_thrust)*nacelle.n_engine       # statistical regression

    nacelle.c_g = nacelle.x_ext + 0.7 * nacelle.length      # statistical regression

    return


