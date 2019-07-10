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

from marilib.aircraft_model.operations import handling_qualities as h_q

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
# Modify initial values here
#======================================================================================================

aircraft.turbofan_engine.reference_thrust = 119000.
aircraft.wing.area = 151.9

#======================================================================================================
# Design process
#======================================================================================================
# Global design parameter n°1 : aircraft.turbofan_engine.reference_thrust, bounds = (50000,150000)
# Global design parameter n°2 : aircraft.wing.area, bounds = (50,200)
#
# Geometrical coupling on : aircraft.horizontal_tail.area     (statistical sizing of HTP area)
# Geometrical coupling on : aircraft.vertical_tail.area       (statistical sizing of VTP area)
#
# Geometrical coupling on : aircraft.turbofan_nacelle.width
# Geometrical coupling on : aircraft.turbofan_nacelle.y_ext
#
# Mass design parameter n°1 : aircraft.weights.mtow
# Mass design parameter n°2 : aircraft.weights.mlw
# Mass design parameter n°3 : aircraft.weights.mzfw
# Mass constraint n°1 : aircraft.weights.mass_constraint_1 ==> 0
# Mass constraint n°2 : aircraft.weights.mass_constraint_2 ==> 0
# Mass constraint n°3 : aircraft.weights.mass_constraint_3 ==> 0
#
# Perfo constraint n°1 : aircraft.high_speed.perfo_constraint_1 >= 0
# Perfo constraint n°2 : aircraft.high_speed.perfo_constraint_2 >= 0
# Perfo constraint n°3 : aircraft.low_speed.perfo_constraint_3 >= 0
# Perfo constraint n°4 : aircraft.high_speed.perfo_constraint_3 >= 0
# Perfo constraint n°5 : aircraft.low_speed.perfo_constraint_1 >= 0
# Perfo constraint n°6 : aircraft.low_speed.perfo_constraint_2 >= 0
#
# Possible criterion : aircraft.weights.mtow
# Possible criterion : aircraft.cost_mission.block_fuel
# Possible criterion : aircraft.environmental_impact.CO2_metric
# Possible criterion : aircraft.economics.cash_operating_cost
# Possible criterion : aircraft.economics.direct_operating_cost

#------------------------------------------------------------------------------------------------------
thrust_bnd = (50000,150000)
area_bnd = (50,200)
search_domain = (thrust_bnd,area_bnd)

# Perform MDF optimization
#------------------------------------------------------------------------------------------------------
criterion = "block_fuel"
mda_type = "MDA2"

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
print("Reference thrust turbofan = ","%.0f"%aircraft.turbofan_engine.reference_thrust," N")
print("Reference thrust effective = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
print("Turbofan mass = ","%.0f"%aircraft.turbofan_nacelle.mass," kg")
print("Cruise SFC = ","%.4f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("Cruise LoD = ","%.4f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
print("-------------------------------------------")
print("Wing area = ","%.2f"%aircraft.wing.area," m2")
print("Wing span = ","%.2f"%aircraft.wing.span," m")
print("-------------------------------------------")
print("Wing position = ","%.2f"%aircraft.wing.x_root," m")
print("HTP area = ","%.2f"%aircraft.horizontal_tail.area," m2")
print("VTP area = ","%.2f"%aircraft.vertical_tail.area," m2")
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
print("Evaluation mission range = ","%.0f"%unit.NM_m(aircraft.cost_mission.range)," NM")
print("Evaluation mission block fuel = ","%.0f"%aircraft.cost_mission.block_fuel," kg")
print("Evaluation mission cash op cost = ","%.0f"%aircraft.economics.cash_operating_cost," $")
print("CO2 metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1000)," kg/km/m0.48")



# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"Design example","This plane")
