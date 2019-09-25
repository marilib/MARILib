#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study illustrates the impact of wing deformation mode on
wing mass and aerodynamics and thus airplane mission performance
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
morphing_i = 1 #aircraft.wing.morphing     # (1: AR mode, 2: Span mode)
wing_span_i = aircraft.wing.span
wing_aspect_ratio_i = aircraft.wing.aspect_ratio
wing_area_i = aircraft.wing.area

mtow_i = aircraft.weights.mtow

print("")
print("Morphing = ","%i"%morphing_i," (1: AR mode, 2: Span mode)")
print("Initial wing span = ","%.2f"%wing_span_i," m")
print("Initial wing aspect ratio = ","%.2f"%wing_aspect_ratio_i," no_dim")
print("Wing area = ","%.2f"%wing_area_i," m2")
print("MTOW input = ","%.0f"%mtow_i," kg")

# Reloading updated values
#------------------------------------------------------------------------------------------------------
aircraft.wing.morphing = morphing_i
aircraft.wing.area = wing_area_i
aircraft.wing.span = wing_span_i
aircraft.wing.aspect_ratio = wing_aspect_ratio_i

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


# Print relevant data
#------------------------------------------------------------------------------------------------------
delta_mtow = mtow_i - mtow_req     # Should be positive or zero

print("")
print("MTOW required = ","%.0f"%mtow_req," kg")
print("Delta MTOW = ","%.2f"%delta_mtow," kg (Should be positive or zero)")
disa = 0
altp = aircraft.design_driver.ref_cruise_altp
mach = aircraft.design_driver.cruise_mach
vtas = earth.vtas_from_mach(altp,disa,mach)

print("")
print("Wing mass = ","%.2f"%aircraft.wing.mass," kg")
print("Wing span = ","%.2f"%aircraft.wing.span," m")
print("Wing aspect ratio = ","%.2f"%aircraft.wing.aspect_ratio," no_dim")

print("")
print("True air speed = ","%.2f"%unit.kt_mps(vtas)," kt")
print("Fuel mission = ","%.2f"%aircraft.nominal_mission.block_fuel," kg")
print("LoD cruise (LoD max) = ","%.2f"%aircraft.aerodynamics.cruise_lod_max," no_dim")
print("SFC cruise = ","%.3f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
show.draw_3d_view(aircraft,"study_n1","This plane")

