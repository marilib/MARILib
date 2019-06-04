#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""


from marilib.tools import units as unit

from marilib.aircraft_model.airplane import viewer as show

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run, initialization as init

#======================================================================================================
# Initialization
#======================================================================================================
propulsive_architecture = 1 # 1:turbofan, 2:partial turboelectric

aircraft = Aircraft(propulsive_architecture)

n_pax_ref = 8
design_range = unit.m_NM(3300)
cruise_mach = 0.79

number_of_engine = init.n_engine()

#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

print("-------------------------------------------")
print("Initialization : done")

#======================================================================================================
# Design process
#======================================================================================================

#------------------------------------------------------------------------------------------------------
thrust_bnd = (20000,35000)
area_bnd = (20,50)
search_domain = (thrust_bnd,area_bnd)

# Perform MDF optimization
#------------------------------------------------------------------------------------------------------
criterion = "MTOW"

run.mdf_process(aircraft,search_domain,criterion)

print("-------------------------------------------")
print("Optimization : done")

#======================================================================================================
# Print some results
#======================================================================================================
print("-------------------------------------------")
print("Number of passengers = ","%.0f"%aircraft.cabin.n_pax_ref," int")
print("Design range = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print("Cruise Mach number = ","%.2f"%aircraft.design_driver.cruise_mach," Mach")
print("-------------------------------------------")
print("Reference thrust = ","%.0f"%aircraft.turbofan_engine.reference_thrust," N")
print("Bypass ratio = ","%.1f"%aircraft.turbofan_engine.bpr," no_dim")
print("-------------------------------------------")
print("Wing area = ","%.2f"%aircraft.wing.area," m2")
print("Wing span = ","%.2f"%aircraft.wing.span," m")
print("-------------------------------------------")
print("Fuselage length = ","%.2f"%aircraft.fuselage.length," m")
print("Fuselage width = ","%.2f"%aircraft.fuselage.width," m")
print("-------------------------------------------")
print("MTOW = ","%.2f"%aircraft.weights.mtow," kg")
print("MLW = ","%.2f"%aircraft.weights.mlw," kg")
print("OWE = ","%.2f"%aircraft.weights.owe," kg")
print("MWE = ","%.2f"%aircraft.weights.mwe," kg")
print("-------------------------------------------")
print("Design range = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print("Effective nominal range = "+"%.0f"%unit.NM_m(aircraft.nominal_mission.range)+" NM")
print("")
print("Take off field length required = "+"%.0f"%aircraft.low_speed.req_tofl+" m")
print("Take off field length effective = "+"%.0f"%aircraft.low_speed.eff_tofl+" m")
print("")
print("Approach speed required = "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
print("Approach speed effective = "+"%.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)+" kt")
print("")
print("One engine ceiling path required = "+"%.1f"%(aircraft.low_speed.req_oei_path*100)+" %")
print("One engine ceiling path effective = "+"%.1f"%(aircraft.low_speed.eff_oei_path*100)+" %")
print("")
print("Climb speed required in MCL rating = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_climb)+" ft/min")
print("Climb speed effective in MCL rating = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)+" ft/min")
print("")
print("Climb speed required in MCR rating = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.req_vz_cruise)+" ft/min")
print("Climb speed effective in MCR rating = "+"%.1f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_cruise)+" ft/min")
print("")
print("Time to climb required = "+"%.1f"%unit.min_s(aircraft.high_speed.req_ttc)+" min")
print("Time to climb effective = "+"%.1f"%unit.min_s(aircraft.high_speed.eff_ttc)+" min")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"Design example","This plane")
