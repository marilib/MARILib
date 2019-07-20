#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy


#===========================================================================================================
def eval_turbofan_pylon_mass(aircraft):
    """
    Turbofan pylon mass & CG estimation
    """

    nacelle = aircraft.turbofan_nacelle
    engine = aircraft.turbofan_engine

    pylon = aircraft.turbofan_pylon

    pylon.mass = 0.0031*engine.reference_thrust*engine.n_engine

    if (engine.n_engine==2):
        pylon.c_g = nacelle.x_ext + 0.75*nacelle.length
    elif (engine.n_engine==4):
        pylon.c_g = 0.5*(nacelle.x_int + nacelle.x_ext) + 0.75*nacelle.length
    else:
        raise Exception("Number of engine is not allowed")

    return


#===========================================================================================================
def eval_turbofan_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    engine = aircraft.turbofan_engine

    engine.rating_factor = {"MTO":0.800, "MCN":0.688, "MCL":0.624, "MCR":0.560, "FID":0.100}

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

    nacelle = aircraft.turbofan_nacelle

    nacelle.width = 0.5 * engine.bpr ** 0.7 + 5.E-6 * engine.reference_thrust

    nacelle.length = 0.86 * nacelle.width + engine.bpr ** 0.37      # statistical regression

    Knac = numpy.pi * nacelle.width * nacelle.length

    nacelle.net_wetted_area = Knac*(1.48 - 0.0076*Knac)*engine.n_engine        # statistical regression

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
def eval_turbofan_nacelle_mass(aircraft):
    """
    Thermal propulsive nacelle mass estimation
    """

    engine = aircraft.turbofan_engine

    nacelle = aircraft.turbofan_nacelle

    nacelle.mass = (1250. + 0.021*engine.reference_thrust)*engine.n_engine       # statistical regression

    nacelle.c_g = nacelle.x_ext + 0.7 * nacelle.length      # statistical regression

    return


