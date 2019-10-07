#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry

This study allows to compare span driven and aspect ratio driven optimum design
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
aircraft.wing.morphing = 2     # (1: AR is fixed, 2: Span is fixed)

aircraft.wing.aspect_ratio = 9
aircraft.wing.span = 32.


# Automatic optimization of the airplane
#------------------------------------------------------------------------------------------------------
thrust_bnd = (110000,150000)
area_bnd = (100,200)
search_domain = (thrust_bnd,area_bnd)

criterion = "MTOW"       # Criterion to be chosen among  "MTOW", "block_fuel", "CO2_metric", "COC", "DOC"
mda_type = "MDA2"

run.mdf_process(aircraft,search_domain,criterion,mda_type)


# Print relevant output
#------------------------------------------------------------------------------------------------------
print("")
print("Engine thrust = ","%.1f"%(aircraft.propulsion.reference_thrust_effective/10)," daN")
print("Engine mass = ","%.1f"%(aircraft.propulsion.mass)," kg")
print("SFC cruise = ","%.3f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("")
print("Wing area = ","%.1f"%aircraft.wing.area," m2")
print("Wing span = ","%.1f"%aircraft.wing.span," m")
print("Wing aspect ratio = ","%.1f"%aircraft.wing.aspect_ratio," no_dim")
print("Wing mass = ","%.1f"%aircraft.wing.mass," kg")
print("LoD cruise (LoD max) = ","%.2f"%aircraft.aerodynamics.cruise_lod_max," no_dim")
print("")
print("MTOW = ","%.0f"%aircraft.weights.mtow," kg")
print("OWE = ","%.0f"%aircraft.weights.owe," kg")

print("")
print("Cash Operating Cost = ","%.1f"%aircraft.economics.cash_operating_cost," $/trip")
print("Cost mission block fuel = ","%.1f"%(aircraft.cost_mission.block_fuel)," kg/trip")
print("Carbon dioxide emission = ","%.1f"%(aircraft.cost_mission.block_CO2)," kg/trip")

print("")
print("Take off field length required = "+"%.1f"%aircraft.low_speed.req_tofl+" m")
print("Take off field length effective = "+"%.1f"%aircraft.low_speed.eff_tofl+" m")
print("")
print("Approach speed required = "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
print("Approach speed effective = "+"%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)+" kt")
print("")
print("Vertical speed required = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print("Vertical speed effective = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Time to climb required = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print("Time to climb effective = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
#show.draw_3d_view(aircraft,"study_n4","This plane")

