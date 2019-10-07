#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run

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

run.eval_geometrical_analysis(aircraft)

run.eval_mass_breakdown(aircraft)

run.eval_handling_quality_analysis(aircraft)

run.eval_performance_analysis(aircraft)

run.eval_payload_range_analysis(aircraft)


#------------------------------------------------------------------------------------------------------
print(unit.smart_format(unit.convert_to(aircraft.design_driver.INFO["design_range"]["unit"],
                                        aircraft.design_driver.design_range)),
      aircraft.design_driver.INFO["design_range"]["unit"],
      ":",aircraft.design_driver.INFO["design_range"]["txt"])

print(unit.smart_format(unit.convert_to(aircraft.wing.INFO["area"]["unit"],
                                        aircraft.wing.area)),
      aircraft.wing.INFO["area"]["unit"],
      ":",aircraft.wing.INFO["area"]["txt"])

print(unit.smart_format(unit.convert_to(aircraft.propulsion.INFO["reference_thrust_effective"]["unit"],
                                        aircraft.propulsion.reference_thrust_effective)),
      aircraft.propulsion.INFO["reference_thrust_effective"]["unit"],
      ":",aircraft.propulsion.INFO["reference_thrust_effective"]["txt"])


#------------------------------------------------------------------------------------------------------
aircraft.export_to_file(filename = "aircraft_data.txt", write_detail = True)

