#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

from marilib import numpy

from marilib.tools import units as unit

from marilib.aircraft_model.airplane import viewer as show

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run, initialization as init

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

aircraft.propulsion.reference_thrust = 120000.
aircraft.wing.area = 149

e_power = 1.0e6       # Watts, electric motor power

result = numpy.array([["e-fan power              (kW)"],
                      ["e-chain mass             (kg)"],
                      ["TF ref thrust           (daN)"],
                      ["Effective ref thrust    (daN)"],
                      ["TF nacelle mass          (kg)"],
                      ["Wing area                (m2)"],
                      ["Wing span                 (m)"],
                      ["Fuselage length           (m)"],
                      ["Fuselage width            (m)"],
                      ["MTOW                     (kg)"],
                      ["MLW                      (kg)"],
                      ["OWE                      (kg)"],
                      ["MWE                      (kg)"],
                      ["Cruise SFC         (kg/daN/h)"],
                      ["Cruise L/D           (no_dim)"],
                      ["Take Off Field Length     (m)"],
                      ["Approach speed           (kt)"],
                      ["One engine inop path      (%)"],
                      ["Vz TOC MCL rating    (ft/min)"],
                      ["Vz TOC MCR rating    (ft/min)"],
                      ["Time to climb           (min)"],
                      ["Block fuel               (kg)"],
                      ["Cash Op Cost         ($/trip)"],
                      ["CO2 metric (10e-4kg/km/m0.48)"]])

for e_power in (0.05e6, 0.15e6, 0.25e6, 0.5e6, 1.0e6, 1.5e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6, 4.0e6):

    aircraft.pte1_power_elec_chain.mto = e_power
    aircraft.pte1_power_elec_chain.mcn = e_power
    aircraft.pte1_power_elec_chain.mcl = e_power
    aircraft.pte1_power_elec_chain.mcr = e_power

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

    #------------------------------------------------------------------------------------------------------
    thrust_bnd = (50000,150000)
    area_bnd = (50,200)
    search_domain = (thrust_bnd,area_bnd)

    # Perform MDF optimization
    #------------------------------------------------------------------------------------------------------
    criterion = "Block_fuel"
    mda_type = "MDA2"

    run.mdf_process(aircraft,search_domain,criterion,mda_type)

    print("-------------------------------------------")
    print("Optimization : done")

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

    res = numpy.array([["%8.0f"%(e_power/1000.)],
                       ["%8.0f"%global_e_mass],
                       ["%8.0f"%(aircraft.propulsion.reference_thrust/10)],
                       ["%8.0f"%(aircraft.propulsion.reference_thrust_effective/10)],
                       ["%8.0f"%aircraft.turbofan_nacelle.mass],
                       ["%8.1f"%aircraft.wing.area],
                       ["%8.1f"%aircraft.wing.span],
                       ["%8.1f"%aircraft.fuselage.length],
                       ["%8.1f"%aircraft.fuselage.width],
                       ["%8.0f"%aircraft.weights.mtow],
                       ["%8.0f"%aircraft.weights.mlw],
                       ["%8.0f"%aircraft.weights.owe],
                       ["%8.0f"%aircraft.weights.mwe],
                       ["%8.4f"%(aircraft.propulsion.sfc_cruise_ref*36000)],
                       ["%8.4f"%(aircraft.aerodynamics.cruise_lod_max)],
                       ["%8.0f"%aircraft.low_speed.eff_tofl],
                       ["%8.1f"%unit.kt_mps(aircraft.low_speed.eff_app_speed)],
                       ["%8.1f"%(aircraft.low_speed.eff_oei_path*100)],
                       ["%8.0f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_climb)],
                       ["%8.0f"%unit.ftpmin_mps(aircraft.high_speed.eff_vz_cruise)],
                       ["%8.1f"%unit.min_s(aircraft.high_speed.eff_ttc)],
                       ["%8.0f"%aircraft.cost_mission.block_fuel],
                       ["%8.0f"%aircraft.economics.cash_operating_cost],
                       ["%8.0f"%(aircraft.environmental_impact.CO2_metric*10000000)]])

    result = numpy.hstack([result,res])

numpy.savetxt("scan_result.txt",result,delimiter=" ;",fmt='%s')

#======================================================================================================
# Print some results
#======================================================================================================
# print("-------------------------------------------")
print("Number of passengers = ","%.0f"%aircraft.cabin.n_pax_ref," int")
print("Design range = ","%.0f"%unit.NM_m(aircraft.design_driver.design_range)," NM")
print("Cruise Mach number = ","%.2f"%aircraft.design_driver.cruise_mach," Mach")
print("-------------------------------------------")
print("Reference thrust turbofan = ","%.0f"%aircraft.propulsion.reference_thrust," N")
print("Reference thrust effective = ","%.0f"%aircraft.propulsion.reference_thrust_effective," N")
print("Cruise SFC = ","%.4f"%(aircraft.propulsion.sfc_cruise_ref*36000)," kg/daN/h")
print("Cruise LoD = ","%.4f"%(aircraft.aerodynamics.cruise_lod_max)," no_dim")
print("-------------------------------------------")
print("Turbofan mass = ","%.0f"%aircraft.turbofan_nacelle.mass," kg")
print("Global mass of electric chain = ","%.0f"%global_e_mass," kg")
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
