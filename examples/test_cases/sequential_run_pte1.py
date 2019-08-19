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
propulsive_architecture = "PTE1" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 150
design_range = unit.m_NM(3000)
cruise_mach = 0.78

# initialize aircraft
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

print("-------------------------------------------")
print("Initialization : done")

#======================================================================================================
# Sequential process
#======================================================================================================
# This sequence goes through all models without doing any solving

# Solve the geometric coupling between airframe and engines
#------------------------------------------------------------------------------------------------------
run.eval_aircraft_geom_analysis(aircraft)

# Estimate all mass and CGs with or without Mass-Mission adaptation
#------------------------------------------------------------------------------------------------------
run.eval_mass_breakdown(aircraft)

# Calculate Payload-Range diagram
#------------------------------------------------------------------------------------------------------
run.eval_payload_range_analysis(aircraft)

# Calculate all airplane performances
#------------------------------------------------------------------------------------------------------
run.eval_performance_analysis(aircraft)

# Handling quality analysis
#------------------------------------------------------------------------------------------------------
run.eval_handling_quality_analysis(aircraft)

print("-------------------------------------------")
print("Sequence : done")

# Print output aircraft
#------------------------------------------------------------------------------------------------------
aircraft.export_to_file(filename="aircraft_data.txt", write_detail=True)

# Draw 3D view
#------------------------------------------------------------------------------------------------------
show.draw_3d_view(aircraft,"Design example","This plane")
