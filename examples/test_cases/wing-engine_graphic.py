#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

import numpy as np

from marilib.tools import units as unit

from marilib.aircraft_data.aircraft_description import Aircraft

from marilib.processes import assembly as run

from marilib.aircraft_model.airplane.viewer import draw_design_space

from marilib.processes.assembly import explore_design_space


#======================================================================================================
# Initialization
#======================================================================================================
propulsive_architecture = "TF" # TF:turbofan, PTE1:partial turboelectric 1
number_of_engine = 2

aircraft = Aircraft()

n_pax_ref = 150
design_range = unit.m_NM(3000)
cruise_mach = 0.78

# Build initial point
#------------------------------------------------------------------------------------------------------
run.aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propulsive_architecture, number_of_engine)


# Perform MDF optimization
#------------------------------------------------------------------------------------------------------
thrust_bnd = (50000,150000)
area_bnd = (50,200)
search_domain = (thrust_bnd,area_bnd)

criterion = "Block_fuel"
mda_type = "MDA2"

run.mdf_process(aircraft,search_domain,criterion,mda_type)

print("-------------------------------------------")
print("Optimization : done")
print("")
print("Engine thrust = ","%.1f"%(unit.daN_N(aircraft.propulsion.reference_thrust))," daN")
print("Wing area = ","%.1f"%aircraft.wing.area," m2")


# Compute grid points
#------------------------------------------------------------------------------------------------------
res = [aircraft.propulsion.reference_thrust,
       aircraft.wing.area]
step = [0.05, 0.05]    # Relative grid step
file = "explore_design.txt"

explore_design_space(aircraft, res, step, mda_type, file)


# Draw graphic
#------------------------------------------------------------------------------------------------------
field = criterion
const = ['TOFL', 'App_speed', 'OEI_path', 'Vz_MCL', 'Vz_MCR', 'TTC']
color = ['red', 'blue', 'violet', 'orange', 'brown', 'yellow']
limit = [aircraft.low_speed.req_tofl,
         unit.kt_mps(aircraft.low_speed.req_app_speed),
         unit.pc_no_dim(aircraft.low_speed.req_oei_path),
         unit.ftpmin_mps(aircraft.high_speed.req_vz_climb),
         unit.ftpmin_mps(aircraft.high_speed.req_vz_cruise),
         unit.min_s(aircraft.high_speed.req_ttc)]       # Limit values
bound = np.array(["ub", "ub", "lb", "lb", "lb", "ub"])                 # ub: upper bound, lb: lower bound

draw_design_space(file, res, field, const, color, limit, bound)


