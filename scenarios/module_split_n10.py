#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

------------------------------------------------------------------------------------------------------------
This scenario allows to play with GEMS a full design process with statistical empennage sizing

All processes must be managed at MDO level

Global design parameter n°1 : aircraft.turbofan_engine.reference_thrust, bounds = (50000,150000)
Global design parameter n°2 : aircraft.wing.area, bounds = (50,200)

Circular dependencies on : aircraft.turbofan_nacelle.width
                         : aircraft.turbofan_nacelle.y_ext
                         : aircraft.horizontal_tail.area
                         : aircraft.vertical_tail.area

Mass design parameter n°1 : aircraft.weights.mtow
Mass design parameter n°2 : aircraft.weights.mlw
Mass design parameter n°3 : aircraft.weights.mzfw
Mass constraint n°1 : aircraft.weights.mass_constraint_1 ==> 0
Mass constraint n°2 : aircraft.weights.mass_constraint_2 ==> 0
Mass constraint n°3 : aircraft.weights.mass_constraint_3 ==> 0

Perfo constraint n°1 : aircraft.high_speed.perfo_constraint_1 > 0
Perfo constraint n°2 : aircraft.high_speed.perfo_constraint_2 > 0
Perfo constraint n°3 : aircraft.low_speed.perfo_constraint_3 > 0
Perfo constraint n°4 : aircraft.high_speed.perfo_constraint_3 > 0
Perfo constraint n°5 : aircraft.low_speed.perfo_constraint_1 > 0
Perfo constraint n°6 : aircraft.low_speed.perfo_constraint_2 > 0

Possible criteria : aircraft.weights.mtow
                  : aircraft.cost_mission.block_fuel
                  : aircraft.environmental_impact.CO2_metric
                  : aircraft.economics.cash_operating_cost
                  : aircraft.economics.direct_operating_cost
------------------------------------------------------------------------------------------------------------
"""

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.airplane.airframe.airframe_design \
    import eval_cabin_design, eval_fuselage_design, eval_vtp_design, eval_vtp_statistical_sizing, \
           eval_htp_design, eval_htp_statistical_sizing, eval_wing_design

from marilib.airplane.propulsion.propulsion_design \
    import eval_propulsion_design

from marilib.aircraft_model.airplane.airplane_design \
    import eval_aerodynamics_design, eval_mass_coupling

from marilib.aircraft_model.operations.mission \
    import eval_nominal_mission, eval_cost_mission

from marilib.processes.component \
    import eval_take_off_performances, eval_landing_performances, eval_climb_performances, \
           eval_oei_performances, eval_co2_metric, eval_economics

from marilib.processes.assembly \
    import aircraft_initialize, eval_mass_breakdown, eval_payload_range_analysis


#-----------------------------------------------------------------------------------------------------------
def aircraft_initialization(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine):
    aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

#-----------------------------------------------------------------------------------------------------------
def fuselage_design(aircraft):
    eval_cabin_design(aircraft)
    eval_fuselage_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def lifting_plane_design(aircraft):
    eval_wing_design(aircraft)
    eval_vtp_design(aircraft)
    eval_htp_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def geometry_coupling(aircraft):
    eval_vtp_statistical_sizing(aircraft)
    eval_htp_statistical_sizing(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def propulsion(aircraft):
    eval_propulsion_design(aircraft)

#-----------------------------------------------------------------------------------------------------------
def aircraft_aerodynamics(aircraft):
    eval_aerodynamics_design(aircraft)

#-----------------------------------------------------------------------------------------------------------
def aircraft_mass(aircraft):
    eval_mass_breakdown(aircraft)

#-----------------------------------------------------------------------------------------------------------
def mass_coupling(aircraft):
    eval_mass_coupling(aircraft)

#-----------------------------------------------------------------------------------------------------------
def nominal_mission(aircraft):
    eval_nominal_mission(aircraft)

#-----------------------------------------------------------------------------------------------------------
def climb_performances(aircraft):
    eval_climb_performances(aircraft)
    eval_oei_performances(aircraft)

#-----------------------------------------------------------------------------------------------------------
def low_speed_performances(aircraft):
    eval_take_off_performances(aircraft)
    eval_landing_performances(aircraft)

#-----------------------------------------------------------------------------------------------------------
def co2_metric(aircraft):
    eval_co2_metric(aircraft)

#-----------------------------------------------------------------------------------------------------------
def economics(aircraft):
    eval_cost_mission(aircraft)
    eval_economics(aircraft)

#-----------------------------------------------------------------------------------------------------------
def payload_range_analysis(aircraft):
    eval_payload_range_analysis(aircraft)


# Initialize aircraft data structure
#---------------------------------------------------------------------------
aircraft = Aircraft()

n_pax_ref = 150                     # Reference number of passengers
design_range = unit.m_NM(3000)      # Design range
cruise_mach = 0.78                  # Nominal cruise mach number

propu_config = "TF"    # "TF": turbofan, "PTE1": partial turbo electric
n_engine = 2           # Number of engine


aircraft_initialization(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)

fuselage_design(aircraft)

lifting_plane_design(aircraft)

propulsion(aircraft)

geometry_coupling(aircraft)

aircraft_aerodynamics(aircraft)

aircraft_mass(aircraft)

mass_coupling(aircraft)

nominal_mission(aircraft)

climb_performances(aircraft)

low_speed_performances(aircraft)

co2_metric(aircraft)

economics(aircraft)

payload_range_analysis(aircraft)


