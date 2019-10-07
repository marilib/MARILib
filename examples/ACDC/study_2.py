#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study illustrates the mass mission adaptation process
"""

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.earth import environment as earth

from marilib.processes import assembly as run

from marilib.aircraft_model.operations import mission as miss

from marilib.aircraft_model.airplane import viewer as show

# Initialize aircraft data structure
#------------------------------------------------------------------------------------------------------
aircraft = Aircraft()

design_range = unit.m_NM(3000)
cruise_mach = 0.78
n_pax_ref = 150

propu_config = "TF" # "TF": turbofan, "PTE1": partial turboelectric
n_engine = 2

# Initialize all input data
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine)


# Possibility to modify initial values
#------------------------------------------------------------------------------------------------------
mtow_i = 72000.  # aircraft.weights.mtow

# Reloading updated values
#------------------------------------------------------------------------------------------------------
aircraft.weights.mtow = mtow_i


# Solve the geometric coupling between airframe and engines
#------------------------------------------------------------------------------------------------------
run.eval_aircraft_pre_design(aircraft)

# Estimate all mass and CGs
#------------------------------------------------------------------------------------------------------
run.eval_mass_estimation(aircraft)

# Evaluate required MTOW to satisfy nominal mission range
#------------------------------------------------------------------------------------------------------
miss.eval_nominal_mission(aircraft)

mtow_req = aircraft.weights.owe + aircraft.nominal_mission.payload + aircraft.nominal_mission.total_fuel


# Print relevant output
#------------------------------------------------------------------------------------------------------
print("")
print("MTOW input = ","%.2f"%mtow_i," kg")
print("OWE structure = ","%.2f"%aircraft.weights.owe," kg")

print("")
print("Total mission fuel = ","%.2f"%aircraft.nominal_mission.total_fuel," kg")
print("Payload = ","%.3f"%aircraft.nominal_mission.payload," kg")
print("MTOW required = ","%.2f"%mtow_req," kg")
print("OWE mission = ","%.2f"%(mtow_i-aircraft.nominal_mission.payload-aircraft.nominal_mission.total_fuel)," kg")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n2","This plane")

