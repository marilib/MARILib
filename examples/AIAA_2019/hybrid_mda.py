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
from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_models import pte1_thrust

from marilib.aircraft_model.airplane import viewer as show

#======================================================================================================
# Initialization
#======================================================================================================
propulsive_architecture = "PTE1" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 150
design_range = unit.m_NM(3000)
cruise_mach = 0.78

#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

#======================================================================================================
# Modify initial values here
#======================================================================================================

aircraft.propulsion.reference_thrust = 118482.
aircraft.wing.area = 152.9

e_power = 1.0e6       # Watts, electric motor power

aircraft.rear_electric_engine.mto_r_shaft_power = e_power
aircraft.rear_electric_engine.mcn_r_shaft_power = e_power
aircraft.rear_electric_engine.mcl_r_shaft_power = e_power
aircraft.rear_electric_engine.mcr_r_shaft_power = e_power

aircraft.propulsion.bli_effect = 1                      #init.boundary_layer_effect()
aircraft.pte1_power_elec_chain.overall_efficiency = 0.90     # 0.90 from init.e_chain_efficiency()

pw_density_factor = 1.

aircraft.pte1_power_elec_chain.generator_pw_density = init.generator_power_density() * pw_density_factor
aircraft.pte1_power_elec_chain.rectifier_pw_density = init.rectifier_pw_density() * pw_density_factor
aircraft.pte1_power_elec_chain.wiring_pw_density = init.wiring_pw_density() * pw_density_factor
aircraft.pte1_power_elec_chain.cooling_pw_density = init.cooling_pw_density() * pw_density_factor

aircraft.rear_electric_nacelle.controller_pw_density = init.controller_pw_density() * pw_density_factor
aircraft.rear_electric_nacelle.motor_pw_density = init.e_motor_pw_density() * pw_density_factor
aircraft.rear_electric_nacelle.nacelle_pw_density = init.e_nacelle_pw_density() * pw_density_factor

#======================================================================================================
# Design process
#======================================================================================================

# Solve the geometric coupling between airframe and engines
#------------------------------------------------------------------------------------------------------
run.eval_aircraft_pre_design(aircraft)

# Estimate all mass and CGs
#------------------------------------------------------------------------------------------------------
#run.eval_mass_estimation(aircraft)
run.eval_mass_mission_adaptation(aircraft)

# Calculate all airplane performances
#------------------------------------------------------------------------------------------------------
run.eval_performance_analysis(aircraft)

# Handling quality analysis
#------------------------------------------------------------------------------------------------------
run.eval_handling_quality_analysis(aircraft)


# Some performances about electric chain
#------------------------------------------------------------------------------------------------------
gen_pwd = aircraft.pte1_power_elec_chain.generator_pw_density
rec_pwd = aircraft.pte1_power_elec_chain.rectifier_pw_density
wire_pwd = aircraft.pte1_power_elec_chain.wiring_pw_density
cool_pwd = aircraft.pte1_power_elec_chain.cooling_pw_density

cont_pwd = aircraft.rear_electric_nacelle.controller_pw_density
mot_pwd = aircraft.rear_electric_nacelle.motor_pw_density
nac_pwd = aircraft.rear_electric_nacelle.nacelle_pw_density

global_e_mass = (1/gen_pwd + 1/rec_pwd + 1/wire_pwd + 1/cool_pwd + 1/cont_pwd + 1/mot_pwd + 1/nac_pwd)*e_power

#------------------------------------------------------------------------------------------------------
if (propulsive_architecture=="PTE1"):

    shaft_power = aircraft.rear_electric_engine.mcr_r_shaft_power

    disa = 0.
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    (pamb,tamb,tstd,dt_o_dz) = earth.atmosphere(altp,disa)

    (e_fan_thrust,q_air,dv_bli) = jet.fan_thrust_with_bli(aircraft.rear_electric_nacelle,pamb,tamb,mach,shaft_power)

    vair = mach*earth.sound_speed(tamb)

    kVbli = dv_bli / vair

    #------------------------------------------------------------------------------------------------------
    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code

    throttle = 1.
    nei = 0

    fn,sfc,sec,data = pte1_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)

    (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0,sfc0) = data

    sfc_factor = sfc/sfc0   # factor on cruise SFC due to rear fuselage electric nacelle with bli

    # Print some results
    #------------------------------------------------------------------------------------------------------
    print("-------------------------------------------")
    print("Global mass of electric chain = ","%.0f"%global_e_mass," kg")
    print("Electric fan length = ","%.2f"%aircraft.rear_electric_nacelle.length," m")
    print("Electric fan diameter = ","%.2f"%aircraft.rear_electric_nacelle.fan_width," m")
    print("Electric nozzle diameter = ","%.2f"%aircraft.rear_electric_nacelle.nozzle_width," m")
    print("relative decrease of e-fan inlet velocity in cruise = ","%.3f"%kVbli," no_dim")
    print("Factor on e-fan thrust due to BLI in cruise = ","%.3f"%aircraft.propulsion.bli_r_thrust_factor," no_dim")
    print("Global factor on SFC in cruise = ","%.4f"%sfc_factor," no_dim")

print("-------------------------------------------")
print("Number of passengers = ","%.0f"%aircraft.cabin.n_pax_ref," int")
print("Design range = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print("Cruise Mach number = ","%.2f"%aircraft.design_driver.cruise_mach," Mach")
print("-------------------------------------------")
print("Reference thrust turbofan = ","%.0f"%aircraft.propulsion.reference_thrust," N")
print("Reference thrust effective = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
print("Turbofan mass = ","%.0f"%aircraft.turbofan_nacelle.mass," kg")
print("Cruise SFC = ","%.4f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
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
print("Evaluation mission range = ","%.0f"%unit.NM_m(aircraft.cost_mission.range)," NM")
print("Evaluation mission block fuel = ","%.0f"%aircraft.cost_mission.block_fuel," kg")
print("Evaluation mission cash op cost = ","%.0f"%aircraft.economics.cash_operating_cost," $")
print("CO2 metric = ","%.4f"%(aircraft.environmental_impact.CO2_metric*1000)," kg/km/m0.48")

# airplane 3D view
#------------------------------------------------------------------------------------------------------
print("-------------------------------------------")
print("3 view drawing : launched")

show.draw_3d_view(aircraft,"Design example","This plane")
