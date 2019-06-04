#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import initialization as init

from marilib.airplane.airframe.airframe_design \
    import eval_cabin_design, eval_fuselage_design, eval_vtp_design, eval_htp_design, eval_wing_design

from marilib.airplane.propulsion.propulsion_design \
    import eval_propulsion_design

from marilib.aircraft_model.airplane.airplane_design \
    import eval_aerodynamics_design

from marilib.processes.assembly \
    import aircraft_initialize, eval_mass_breakdown, eval_climb_performances, \
           eval_payload_range_analysis, eval_handling_quality_analysis, eval_hq0, eval_mda0

from marilib.processes.component \
    import eval_nominal_mission, eval_take_off_performances, eval_landing_performances, \
           eval_co2_metric, eval_cost_mission, eval_economics


#-----------------------------------------------------------------------------------------------------------
def aircraft_initialization(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine):
    aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

#-----------------------------------------------------------------------------------------------------------
def fuselage_design(aircraft):
    eval_cabin_design(aircraft)
    eval_fuselage_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def predesign_initialization(aircraft):
    # Variables :
    # aircraft.turbofan_nacelle.width
    # aircraft.turbofan_nacelle.y_ext
    # Must be initialized before running lifting_plane_design
    #---------------------------------------------------------------------------
    bpr = aircraft.turbofan_engine.bpr
    reference_thrust = aircraft.turbofan_engine.reference_thrust
    aircraft.turbofan_nacelle.width = init.turbofan_nacelle_width(bpr,reference_thrust)

    nacelle_attachment = aircraft.turbofan_nacelle.attachment
    fuselage_width = aircraft.fuselage.width
    nacelle_width = aircraft.turbofan_nacelle.width
    aircraft.turbofan_nacelle.y_ext = init.turbofan_nacelle_y_ext(nacelle_attachment,fuselage_width,nacelle_width)
    return

#-----------------------------------------------------------------------------------------------------------
def lifting_plane_design(aircraft):
    eval_vtp_design(aircraft)
    eval_wing_design(aircraft)
    eval_htp_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def propulsion(aircraft):
    eval_propulsion_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def aircraft_aerodynamics(aircraft):
    eval_aerodynamics_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def aircraft_mass(aircraft):
    eval_mass_breakdown(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def handling_quality_analysis(aircraft):
#    eval_handling_quality_analysis(aircraft)    # Computes CG constraints only
    eval_hq0(aircraft)                          # Compute Wing X position, HTP & VTP areas without solving mass constraints

#-----------------------------------------------------------------------------------------------------------
def nominal_mission(aircraft):
    eval_nominal_mission(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def performance_analysis(aircraft):
    eval_take_off_performances(aircraft)
    eval_climb_performances(aircraft)
    eval_landing_performances(aircraft)
    eval_payload_range_analysis(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def criteria(aircraft):
    eval_co2_metric(aircraft)
    eval_cost_mission(aircraft)
    eval_economics(aircraft)
    return



# Initialize aircraft data structure
#---------------------------------------------------------------------------
aircraft = Aircraft()

n_pax_ref = 150                     # Reference number of passengers
design_range = unit.m_NM(3000)      # Design range
cruise_mach = 0.78                  # Nominal cruise mach number

propu_config = 1    # 1: turbofan, 2: partial turbo electric
n_engine = 2        # Number of engine

aircraft_initialization(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

fuselage_design(aircraft)

predesign_initialization(aircraft)

lifting_plane_design(aircraft)

propulsion(aircraft)

aircraft_aerodynamics(aircraft)

aircraft_mass(aircraft)

nominal_mission(aircraft)

handling_quality_analysis(aircraft)

performance_analysis(aircraft)

criteria(aircraft)


