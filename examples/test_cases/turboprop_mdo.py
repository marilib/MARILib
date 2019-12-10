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
propulsive_architecture = "TP" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 70
design_range = unit.m_NM(600)
cruise_mach = 0.44

# initialize aircraft
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

print("-------------------------------------------")
print("Initialization : done")

aircraft.propulsion.reference_thrust = 45000.
aircraft.wing.area = 65.

#------------------------------------------------------------------------------------------------------
thrust_bnd = (30000,100000)
area_bnd = (40,100)
search_domain = (thrust_bnd,area_bnd)

# Perform MDF optimization
#------------------------------------------------------------------------------------------------------
criterion = "MTOW"
mda_type = "MDA2"

#run.eval_mda2(aircraft)

run.mdf_process(aircraft,search_domain,criterion,mda_type)

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
print("Reference thrust = ","%.0f"%aircraft.propulsion.reference_thrust," N")
#print("Reference shaft power turboprop = ","%.1f"%(aircraft.turboprop_engine.reference_power/1000.)," kW")
print("Reference thrust effective = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
print("Cruise SFC = ","%.4f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("Cruise SEC = ","%.4f"%(aircraft.propulsion.sec_cruise_ref/100)," kW/daN")
print("Cruise LoD = ","%.4f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
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
print("Approach speed required= "+"%.1f"%unit.kt_mps(aircraft.low_speed.req_app_speed)+" kt")
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
print("-------------------------------------------")
print("Evaluation mission range = ","%.1f"%unit.NM_m(aircraft.cost_mission.range)," NM")
print("Evaluation mission battery mass = ","%.3f"%aircraft.cost_mission.req_battery_mass," kg")
print("Evaluation mission cash op cost = ","%.0f"%aircraft.economics.cash_operating_cost," $")
print("CO2 metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1000)," kg/km/m0.48")



# Print output aircraft
#------------------------------------------------------------------------------------------------------
aircraft.export_to_file(filename="aircraft_data.txt", write_detail=True)

# Draw 3D view
#------------------------------------------------------------------------------------------------------
show.draw_3d_view(aircraft,"Design example","This plane")
