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
propulsive_architecture = "TF" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 150
design_range = unit.m_NM(3000)
cruise_mach = 0.78

txt = numpy.array([["Cruise Mach number     (mach)"],
                   ["Cruise altitude          (ft)"],
                   ["Effective ref thrust    (daN)"],
                   ["Wing area                (m2)"],
                   ["Wing span                 (m)"],
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
                   ["CO2 metric (10e-3kg/km/m0.48)"]])

#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)

# Initialization of design variables
aircraft.turbofan_engine.reference_thrust = 150000.
aircraft.wing.area = 150

#======================================================================================================
# Modify initial values here
#======================================================================================================

for mach in (0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86):

    for altp in (20000., 24000., 28000., 32000., 36000., 40000.):

        aircraft.design_driver.cruise_mach = mach
        aircraft.design_driver.ref_cruise_altp = unit.m_ft(altp)

        print("-------------------------------------------")
        print("Doing optimization for : mach = ",mach,"    altp = ",altp, " ft")
        print("")

       # Perform MDF optimization
        #------------------------------------------------------------------------------------------------------
        thrust_bnd = (50000,150000)
        area_bnd = (50,200)
        search_domain = (thrust_bnd,area_bnd)

        criterion = "block_fuel"
        mda_type = "MDA2"

        run.mdf_process(aircraft,search_domain,criterion,mda_type)

        print("Done")

        # Store results
        #------------------------------------------------------------------------------------------------------
        res = numpy.array([["%8.2f"%(aircraft.design_driver.cruise_mach)],
                           ["%8.0f"%unit.ft_m(aircraft.design_driver.ref_cruise_altp)],
                           ["%8.0f"%(aircraft.propulsion.reference_thrust_effective/10)],
                           ["%8.1f"%aircraft.wing.area],
                           ["%8.1f"%aircraft.wing.span],
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
                           ["%8.1f"%(aircraft.environmental_impact.CO2_metric*1000000)]])

        txt = numpy.hstack([txt,res])

#------------------------------------------------------------------------------------------------------
numpy.savetxt("scan_result_1.txt",txt,delimiter=" ;",fmt='%s')

# USE "altp_mach_graphic.py" TO DRAW SOME RESULTS
