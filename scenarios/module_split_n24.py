#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

------------------------------------------------------------------------------------------------------------
This scenario allows to play with GEMS a full design process with Mass - Mission adaptation AND HQ based empennage
sizing treated as an optimization problem

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
Mass constraint n°1 : aircraft.weights.mass_constraint_1 >= 0
Mass constraint n°2 : aircraft.weights.mass_constraint_2 >= 0
Mass constraint n°3 : aircraft.weights.mass_constraint_3 >= 0

Criterion to be used is : aircraft.weights.mtow


HQ design parameter n°1 : aircraft.wing.x_root
HQ design parameter n°2 : aircraft.horizontal_tail.area
HQ design parameter n°3 : aircraft.vertical_tail.area
HQ constraint n°1 : aircraft.center_of_gravity.cg_constraint_1 >= 0
HQ constraint n°2 : aircraft.center_of_gravity.cg_constraint_2 >= 0
HQ constraint n°3 : aircraft.center_of_gravity.cg_constraint_3 >= 0

Possible criteria : aircraft.horizontal_tail.area & aircraft.vertical_tail.area
                  : aircraft.cost_mission.mtow

REMARK :
HQ optimization can be treated as two coupled optimization of HTP area and VTP area (to be minimized)
                               or one single optimization of the MTOW


Perfo constraint n°1 : aircraft.high_speed.perfo_constraint_1 >= 0
Perfo constraint n°2 : aircraft.high_speed.perfo_constraint_2 >= 0
Perfo constraint n°3 : aircraft.low_speed.perfo_constraint_3 >= 0
Perfo constraint n°4 : aircraft.high_speed.perfo_constraint_3 >= 0
Perfo constraint n°5 : aircraft.low_speed.perfo_constraint_1 >= 0
Perfo constraint n°6 : aircraft.low_speed.perfo_constraint_2 >= 0

Possible criteria : aircraft.weights.mtow
                  : aircraft.cost_mission.block_fuel
                  : aircraft.environmental_impact.CO2_metric
                  : aircraft.economics.cash_operating_cost
                  : aircraft.economics.direct_operating_cost

Additionally :
It is possible to experiment bi-level optimization by managing two disciplinary optimizations :

1- Find best cruise altitude for nominal mission (variable : aircraft.nominal_mission.nominal_cruise_altp)
   which minimizes the mission block fuel (variable : aircraft.nominal_mission.block_fuel)
   under the following bound constraints : m_ft(25000) <= aircraft.nominal_mission.nominal_cruise_altp <= m_ft(45000)
   under ceiling constraints which must be kept positive :  (aircraft.nominal_mission.vz_climb_margin >= 0)
                                                            (aircraft.nominal_mission.vz_cruise_margin >= 0)

2- Find best flying mach number when one engine is inoperative (variable : aircraft.low_speed.oei_best_speed)
   which maximizes the fly path (variable : aircraft.low_speed.eff_oei_path)
   under the following bound constraints : 0.25 <= mach <= aircraft.nominal_mission.nominal_cruise_mach

------------------------------------------------------------------------------------------------------------
"""
import marilib

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.airplane.airframe.airframe_design \
    import eval_cabin_design, eval_fuselage_design, eval_vtp_design, \
           eval_htp_design, eval_wing_design

from marilib.airplane.propulsion.propulsion_design \
    import eval_propulsion_design

from marilib.aircraft_model.airplane.airplane_design \
    import eval_aerodynamics_design

from marilib.aircraft_model.operations.mission \
    import eval_nominal_mission, eval_nominal_climb_constraints, eval_cost_mission

from marilib.processes.component \
    import eval_take_off_performances, eval_landing_performances, eval_climb_performances, \
           eval_oei_path, eval_co2_metric, eval_economics

from marilib.processes.assembly \
    import aircraft_initialize, eval_mass_breakdown, eval_payload_range_analysis, \
           eval_handling_quality_analysis


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
def propulsion(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_propulsion_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def aircraft_aerodynamics(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_aerodynamics_design(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def aircraft_mass(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_mass_breakdown(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def handling_quality_analysis(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_handling_quality_analysis(aircraft)

#-----------------------------------------------------------------------------------------------------------
def nominal_mission(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_nominal_mission(aircraft)
    eval_nominal_climb_constraints(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def performance_analysis(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_take_off_performances(aircraft)
    eval_climb_performances(aircraft)
    eval_landing_performances(aircraft)
    eval_payload_range_analysis(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def oei_performance_analysis(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_oei_path(aircraft)
    return

#-----------------------------------------------------------------------------------------------------------
def criteria(aircraft):
    """
    @constants : [rating_factor,rating_code,architecture,fuel_type,flight_data]
    """
    eval_co2_metric(aircraft)
    eval_cost_mission(aircraft)
    eval_economics(aircraft)
    return


if __name__ == "__main__":
    # Initialize aircraft data structure
    #---------------------------------------------------------------------------
    aircraft = Aircraft()
    
    n_pax_ref = 150                     # Reference number of passengers
    design_range = unit.m_NM(3000)      # Design range
    cruise_mach = 0.78                  # Nominal cruise mach number
    
    propu_config = "TF"    # "TF": turbofan, "PTE1": partial turbo electric
    n_engine = 2           # Number of engine
    
    aircraft_initialization(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)
    
    #---------------------------------------------------------------------------
    aircraft.nominal_mission.nominal_cruise_altp = unit.m_ft(35000)
    
    aircraft.low_speed.oei_best_speed = 0.77        # Mach number
    
    
    #---------------------------------------------------------------------------
    # Setting HQ optimization mode
    aircraft.center_of_gravity.cg_range_optimization = 1
    #---------------------------------------------------------------------------
    
    fuselage_design(aircraft)
    
    lifting_plane_design(aircraft)
    
    propulsion(aircraft)
        
    aircraft_aerodynamics(aircraft)
    
    aircraft_mass(aircraft)
        
    handling_quality_analysis(aircraft)
    
    nominal_mission(aircraft)
    
    performance_analysis(aircraft)
    
    oei_performance_analysis(aircraft)
    
    criteria(aircraft)
    
    
    print(aircraft.nominal_mission.block_fuel)
    print(aircraft.nominal_mission.vz_climb_margin)
    print(aircraft.nominal_mission.vz_cruise_margin)
    
    print(aircraft.low_speed.eff_oei_path)
