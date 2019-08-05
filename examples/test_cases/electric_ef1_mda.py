#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""


from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run, initialization as init

from marilib.aircraft_model.operations import handling_qualities as h_q

from marilib.earth import environment as earth

from marilib.airplane.propulsion import jet_models as jet

from marilib.aircraft_model.airplane import viewer as show

#======================================================================================================
# Initialization
#======================================================================================================
propulsive_architecture = "EF1" # TF:turbofan, PTE1:partial turboelectric 1, EF1:electric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 19
design_range = unit.m_NM(100)
cruise_mach = 0.50

#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

#======================================================================================================
# Modify initial values here
#======================================================================================================

aircraft.electrofan_engine.reference_thrust = 50000.
aircraft.wing.area = 75.

aircraft.electrofan_nacelle.rear_nacelle = 1

r_power = 0.5e6       # Watts, electric motor power

aircraft.rear_electric_engine.mto_r_shaft_power = r_power
aircraft.rear_electric_engine.mcn_r_shaft_power = r_power
aircraft.rear_electric_engine.mcl_r_shaft_power = r_power
aircraft.rear_electric_engine.mcr_r_shaft_power = r_power

aircraft.propulsion.bli_effect = 1                      # boundary layer ingestion on rear fan

#======================================================================================================
# Design process
#======================================================================================================

# Solve the geometric coupling between airframe and engines
#------------------------------------------------------------------------------------------------------
run.eval_aircraft_statistical_pre_design(aircraft)

# Estimate all mass and CGs
#------------------------------------------------------------------------------------------------------
#run.eval_mass_estimation(aircraft)
run.eval_mass_mission_adaptation(aircraft)

# Calculate Payload-Range diagram
#------------------------------------------------------------------------------------------------------
run.eval_payload_range_analysis(aircraft)

# Calculate all airplane performances
#------------------------------------------------------------------------------------------------------
run.eval_performance_analysis(aircraft)

# Handling quality analysis
#------------------------------------------------------------------------------------------------------
run.eval_handling_quality_analysis(aircraft)

#------------------------------------------------------------------------------------------------------
# aircraft.export_to_file(filename = "aircraft_data.txt", write_detail=True)

#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("Number of passengers = ","%.0f"%aircraft.cabin.n_pax_ref," int")
print("Design range = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print("Cruise Mach number = ","%.2f"%aircraft.design_driver.cruise_mach," Mach")
print("-------------------------------------------")
print("Reference thrust electrofan = ","%.0f"%aircraft.electrofan_engine.reference_thrust," N")
print("Reference shaft power electrofan = ","%.1f"%(aircraft.electrofan_engine.reference_power/1000.)," kW")
print("Reference thrust effective = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
print("Electrofan mass = ","%.0f"%aircraft.electrofan_nacelle.mass," kg")
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
print("Nominal mission total energy = ","%.3f"%unit.MWh_J(aircraft.nominal_mission.total_enrg)," MWh")
print("Nominal mission battery mass = ","%.3f"%aircraft.nominal_mission.battery_mass," kg")
print("Battery = ","%.2f"%aircraft.weights.battery," kg")
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
print("Evaluation mission block energy = ","%.3f"%unit.MWh_J(aircraft.cost_mission.block_enrg)," MWh")
print("Evaluation mission total energy = ","%.3f"%unit.MWh_J(aircraft.cost_mission.total_enrg)," MWh")
print("Evaluation mission battery mass = ","%.3f"%aircraft.cost_mission.battery_mass," kg")
print("Evaluation mission cash op cost = ","%.0f"%aircraft.economics.cash_operating_cost," $")
print("CO2 metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1000)," kg/km/m0.48")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"Design example","This plane")
