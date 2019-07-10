#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

from marilib.tools import units as unit

from marilib.aircraft_model.airplane import viewer as show

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run

#======================================================================================================
# Initialization
#======================================================================================================
propulsive_architecture = "TF" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 150
design_range = unit.m_NM(3000)
cruise_mach = 0.78

#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

print("-------------------------------------------")
print("Initialization : done")


#======================================================================================================
# Model sequence
#======================================================================================================
# This sequence evaluates all models without any constraint satifaction solving

# Estimate geometric without any constraint satisfaction
#------------------------------------------------------------------------------------------------------
run.eval_aircraft_geom_analysis(aircraft)

# Estimate all mass and CGs without any constraint satisfaction
#------------------------------------------------------------------------------------------------------
run.eval_mass_breakdown(aircraft)

# Handling quality analysis
#------------------------------------------------------------------------------------------------------
run.eval_handling_quality_analysis(aircraft)

# Calculate all airplane performances
#------------------------------------------------------------------------------------------------------
run.eval_performance_analysis(aircraft)

# Calculate Payload-Range diagram
#------------------------------------------------------------------------------------------------------
run.eval_payload_range_analysis(aircraft)

print("-------------------------------------------")
print("Sequence : done")


#======================================================================================================
# Print aircraft to file
#======================================================================================================
aircraft.export_to_file(file = "aircraft_data.txt")

print("-------------------------------------------")
print("aircraft printed to file")

#======================================================================================================
# Print some results
#======================================================================================================
print("-------------------------------------------")
print("turbofan_nacelle_width = ","%.3f"%aircraft.turbofan_nacelle.width," m")
print("turbofan_nacelle_y_ext = ","%.3f"%aircraft.turbofan_nacelle.y_ext," m")
print("htp_area_statistical_design = ","%.3f"%aircraft.horizontal_tail.area," m2")
print("vtp_area_statistical_design = ","%.3f"%aircraft.vertical_tail.area," m2")
print("-------------------------------------------")
print("mass_constraint_1 = ","%.3f"%aircraft.weights.mass_constraint_1," kg")
print("mass_constraint_2 = ","%.3f"%aircraft.weights.mass_constraint_2," kg")
print("mass_constraint_3 = ","%.3f"%aircraft.weights.mass_constraint_3," kg")
print("-------------------------------------------")
print("cg_constraint_1 = ","%.3f"%aircraft.center_of_gravity.cg_constraint_1," m")
print("cg_constraint_2 = ","%.3f"%aircraft.center_of_gravity.cg_constraint_2," m")
print("cg_constraint_3 = ","%.3f"%aircraft.center_of_gravity.cg_constraint_3," m")
print("-------------------------------------------")
print("low_speed_perfo_constraint_1 = ","%.3f"%aircraft.low_speed.perfo_constraint_1," no_dim")
print("low_speed_perfo_constraint_2 = ","%.3f"%aircraft.low_speed.perfo_constraint_2," no_dim")
print("low_speed_perfo_constraint_3 = ","%.3f"%aircraft.low_speed.perfo_constraint_3," no_dim")
print("-------------------------------------------")
print("high_speed_perfo_constraint_1 = ","%.3f"%aircraft.high_speed.perfo_constraint_1," m/s")
print("high_speed_perfo_constraint_2 = ","%.3f"%aircraft.high_speed.perfo_constraint_2," m/s")
print("high_speed_perfo_constraint_3 = ","%.3f"%aircraft.high_speed.perfo_constraint_3," no_dim")

#======================================================================================================
# airplane 3D view
#======================================================================================================
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"Design example","This plane")


