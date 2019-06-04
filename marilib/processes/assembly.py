#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

Changed name from "design.py" to "assembly.py" on 21:05:2019
"""

import numpy as np
from scipy.optimize import fsolve,minimize,SR1, NonlinearConstraint,BFGS
from marilib.tools import units as unit

from marilib.earth import environment as earth
from marilib.aircraft_model.airplane import airplane_design as airplane, regulation as regul

from marilib.airplane.airframe import airframe_design as airframe

from marilib.airplane.propulsion import propulsion_design as propulsion

from marilib.aircraft_model.operations import handling_qualities as h_q, \
                                              mission as perfo

from marilib.processes import component as sub_proc, initialization as init


#===========================================================================================================
def aircraft_initialize(aircraft, n_pax_ref, design_range, cruise_mach, propu_config, n_engine):
    """
    Initialize a generic aircraft
    """

    aircraft.propulsion.architecture = propu_config

    aircraft.propulsion.fuel_type = init.fuel_type()

    aircraft.name = "my_test_airplane"
    aircraft.design_driver.design_range = design_range        # TLR
    aircraft.design_driver.cruise_mach = cruise_mach          # TLR
    aircraft.cabin.n_pax_ref = n_pax_ref                      # TLR

    aircraft.design_driver.ref_cruise_altp = init.ref_cruise_altp(propu_config)        # TLR
    aircraft.design_driver.top_of_climb_altp = init.top_of_climb_altp(propu_config)    # TLR

    aircraft.aerodynamics.hld_conf_clean = init.hld_conf_clean()
    aircraft.aerodynamics.hld_conf_ld = init.hld_conf_ld()

    aircraft.low_speed.altp_tofl = init.altp_tofl()
    aircraft.low_speed.disa_tofl = init.disa_tofl()
    aircraft.low_speed.kvs1g_tofl = regul.kvs1g_min_take_off()       # Regulation
    aircraft.low_speed.req_tofl = init.req_tofl(design_range)        # TLR

    aircraft.low_speed.altp_app_speed = init.altp_app_speed()
    aircraft.low_speed.disa_app_speed = init.disa_app_speed()
    aircraft.low_speed.kvs1g_app_speed = regul.kvs1g_min_landing()   # Regulation
    aircraft.low_speed.req_app_speed = init.req_app_speed(n_pax_ref)          # TLR

    aircraft.low_speed.disa_oei = init.disa_oei()
    aircraft.low_speed.req_oei_path = regul.ceil_oei_min_path(n_engine)     # Regulation
    aircraft.low_speed.req_oei_altp = init.req_oei_altp(propu_config)       # TLR

    aircraft.high_speed.disa_climb = init.disa_climb()
    aircraft.high_speed.req_vz_climb = init.req_vz_climb()           # TLR
    aircraft.high_speed.req_vz_cruise = init.req_vz_cruise()         # TLR
    aircraft.high_speed.req_toc_altp = init.top_of_climb_altp(propu_config)
    aircraft.high_speed.cas1_ttc = init.cas1_ttc(cruise_mach)
    aircraft.high_speed.cas2_ttc = init.cas2_ttc(cruise_mach)
    aircraft.high_speed.req_ttc = init.req_ttc()                     # TLR

    aircraft.cost_mission.disa = init.cost_mission_disa()
    aircraft.cost_mission.range = init.cost_mission_range(design_range)

    aircraft.economics.fuel_price = init.fuel_price()
    aircraft.economics.elec_price = init.elec_price()
    aircraft.economics.battery_price = init.battery_price()
    aircraft.economics.labor_cost = init.labor_cost()
    aircraft.economics.irp = init.irp()
    aircraft.economics.period = init.period()
    aircraft.economics.interest_rate = init.interest_rate()
    aircraft.economics.utilisation = init.utilisation(design_range)

    aircraft.environmental_impact.CO2_index = earth.emission_index("CO2")
    aircraft.environmental_impact.H2O_index = earth.emission_index("H2O")
    aircraft.environmental_impact.SO2_index = earth.emission_index("SO2")
    aircraft.environmental_impact.NOx_index = earth.emission_index("NOx")
    aircraft.environmental_impact.CO_index = earth.emission_index("CO")
    aircraft.environmental_impact.HC_index = earth.emission_index("HC")
    aircraft.environmental_impact.sulfuric_acid_index = earth.emission_index("sulfuric_acid")
    aircraft.environmental_impact.nitrous_acid_index = earth.emission_index("nitrous_acid")
    aircraft.environmental_impact.nitric_acid_index = earth.emission_index("nitric_acid")
    aircraft.environmental_impact.soot_index = earth.emission_index("soot")

    aircraft.weights.mzfw = init.mzfw(n_pax_ref,design_range)
    aircraft.weights.mtow = init.mtow(n_pax_ref,design_range)
    aircraft.weights.mlw = init.mlw(n_pax_ref,aircraft.weights.mtow,aircraft.weights.mzfw)

    aircraft.cabin.n_pax_front = init.n_pax_front(n_pax_ref)
    aircraft.cabin.n_aisle = init.n_aisle(aircraft.cabin.n_pax_front)

    aircraft.payload.m_pax_nominal = init.m_pax_nominal(design_range)      # TLR
    aircraft.payload.m_pax_max = init.m_pax_max(design_range)              # TLR

    aircraft.center_of_gravity.cg_range_optimization = init.cg_range_optimization()

    aircraft.wing.attachment = init.wing_attachment()
    aircraft.wing.morphing = init.wing_morphing()
    aircraft.wing.hld_type = init.hld_type(n_pax_ref)

    aircraft.wing.area = init.wing_area(n_pax_ref,design_range)                                              # Main design variable
    aircraft.wing.aspect_ratio = init.wing_aspect_ratio()
    aircraft.wing.span = init.wing_span(aircraft.wing.area,aircraft.wing.aspect_ratio)

    if (aircraft.propulsion.architecture<4):
        aircraft.turbofan_engine.n_engine = n_engine
        aircraft.turbofan_engine.bpr = init.bpr(n_pax_ref)
        aircraft.turbofan_engine.reference_thrust = init.reference_thrust(n_pax_ref,design_range,n_engine)                                            # Main design variable

        aircraft.turbofan_engine.core_thrust_ratio = init.core_thrust_ratio()
        aircraft.turbofan_engine.core_width_ratio = init.core_width_ratio()
        aircraft.turbofan_engine.core_weight_ratio = init.core_weight_ratio()

        aircraft.turbofan_nacelle.attachment = init.nacelle_attachment(n_pax_ref)
        aircraft.turbofan_nacelle.efficiency_fan = init.efficiency_fan()
        aircraft.turbofan_nacelle.efficiency_prop = init.efficiency_prop()

    if (aircraft.propulsion.architecture==3):
        aircraft.body_nacelle.length = init.nacelle_body_length()
        aircraft.body_nacelle.width = init.nacelle_body_width()
        aircraft.body_nacelle.hub_width = init.nacelle_body_hub_width()

    if (aircraft.propulsion.architecture==4):
        aircraft.turboprop_engine.n_engine = n_engine
        aircraft.turboprop_engine.reference_thrust = init.reference_thrust(n_pax_ref,design_range,n_engine)                                            # Main design variable
        aircraft.turboprop_engine.propeller_efficiency = init.propeller_efficiency()

    aircraft.horizontal_tail.attachment = init.htp_attachment(aircraft.turbofan_nacelle.attachment)

    aircraft.power_elec_chain.mto = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcn = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcl = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.mcr = init.electric_shaft_power()       # Watts, electric motor power
    aircraft.power_elec_chain.fid = 0.01

    aircraft.battery.strategy = init.battery_strategy()
    aircraft.battery.power_feed = init.battery_power_feed()
    aircraft.battery.time_feed = init.battery_time_feed()
    aircraft.battery.energy_cruise = init.battery_energy_cruise()
    aircraft.battery.energy_density = init.battery_energy_density()
    aircraft.battery.power_density = init.battery_power_density()

    aircraft.power_elec_chain.overall_efficiency = init.e_chain_efficiency()
    aircraft.power_elec_chain.generator_pw_density = init.generator_power_density()
    aircraft.power_elec_chain.rectifier_pw_density = init.rectifier_pw_density()
    aircraft.power_elec_chain.wiring_pw_density = init.wiring_pw_density()
    aircraft.power_elec_chain.cooling_pw_density = init.cooling_pw_density()

    aircraft.electric_nacelle.efficiency_fan = init.efficiency_fan()
    aircraft.electric_nacelle.efficiency_prop = init.efficiency_prop()
    aircraft.electric_nacelle.motor_efficiency = init.e_motor_efficiency()
    aircraft.electric_nacelle.controller_efficiency = init.controller_efficiency()
    aircraft.electric_nacelle.controller_pw_density = init.controller_pw_density()
    aircraft.electric_nacelle.motor_pw_density = init.e_motor_pw_density()
    aircraft.electric_nacelle.nacelle_pw_density = init.e_nacelle_pw_density()

    aircraft.propulsion.bli_effect = init.boundary_layer_effect()

    return


#===========================================================================================================
def eval_aircraft_pre_design(aircraft):
    """
    Perform geometrical pre design
    Solves the coupling carried by nacelle geometry
    """

    airframe.eval_cabin_design(aircraft)
    airframe.eval_fuselage_design(aircraft)

    #===========================================================================================================
    def fct_aircraft_pre_design(x_in,aircraft):

        ac = aircraft

        ac.turbofan_nacelle.width = x_in[0]                         # Coupling variable
        ac.turbofan_nacelle.y_ext = x_in[1]                         # Coupling variable

        airframe.eval_vtp_design(ac)
        airframe.eval_wing_design(ac)
        airframe.eval_htp_design(ac)

        propulsion.eval_propulsion_design(ac)

        y_out = np.array([x_in[0] - ac.turbofan_nacelle.width,
                          x_in[1] - ac.turbofan_nacelle.y_ext])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    bpr = aircraft.turbofan_engine.bpr
    reference_thrust = aircraft.turbofan_engine.reference_thrust
    nacelle_attachment = aircraft.turbofan_nacelle.attachment
    fuselage_width = aircraft.fuselage.width

    nacelle_width_i = init.turbofan_nacelle_width(bpr,reference_thrust)
    nacelle_y_ext_i = init.turbofan_nacelle_y_ext(nacelle_attachment,fuselage_width,nacelle_width_i)

    x_ini = np.array([nacelle_width_i,nacelle_y_ext_i])

    fct_arg = aircraft

    output_dict = fsolve(fct_aircraft_pre_design, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.turbofan_nacelle.width = output_dict[0][0]                         # Coupling variable
    aircraft.turbofan_nacelle.y_ext = output_dict[0][1]                         # Coupling variable

    airframe.eval_vtp_design(aircraft)
    airframe.eval_wing_design(aircraft)
    airframe.eval_htp_design(aircraft)

    propulsion.eval_propulsion_design(aircraft)

    airplane.eval_aerodynamics_design(aircraft)

    return


#===========================================================================================================
def eval_mass_breakdown(aircraft):
    """
    Estimate mass and CGs of the airplane
    Takes MTOW,k MZFW & MLW as input
    """

    airframe.eval_cabin_mass(aircraft)
    airframe.eval_fuselage_mass(aircraft)
    airframe.eval_vtp_mass(aircraft)
    airframe.eval_wing_mass(aircraft)
    airframe.eval_htp_mass(aircraft)
    airframe.eval_landing_gear_mass(aircraft)

    propulsion.eval_propulsion_mass(aircraft)
    propulsion.eval_battery_mass(aircraft)
    propulsion.eval_tank_data(aircraft)

    airplane.eval_system_mass(aircraft)
    airplane.eval_payload_mass(aircraft)
    airplane.eval_aircraft_weights(aircraft)

    airplane.eval_aircraft_cg(aircraft)

    return


#===========================================================================================================
def eval_mass_estimation(aircraft):
    """
    Estimate mass and CGs of the airplane
    Takes MTOW as input but solves the coupling carried by MZFW and MLW
    """

    #===========================================================================================================
    def fct_mass(x_in,aircraft):

        aircraft.weights.mlw = x_in[0]      # Coupling variable
        aircraft.weights.mzfw = x_in[1]     # Coupling variable

        eval_mass_breakdown(aircraft)

        y_out = np.array([aircraft.weights.mass_constraint_1,
                          aircraft.weights.mass_constraint_2])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    x_ini = np.array([aircraft.weights.mlw,
                      aircraft.weights.mzfw])

    fct_arg = aircraft

    output_dict = fsolve(fct_mass, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.weights.mlw = output_dict[0][0]                          # Coupling variable
    aircraft.weights.mzfw = output_dict[0][1]                         # Coupling variable

    # Update mass
    #------------------------------------------------------------------------------------------------------
    eval_mass_breakdown(aircraft)

    return


#===========================================================================================================
def eval_mass_mission_adaptation(aircraft):
    """
    Perform mass - mission adaptation and update mass and CGs
    """

    #===========================================================================================================
    def fct_mass_mission(x_in,aircraft):

        aircraft.weights.mtow = x_in[0]             # Coupling variable
        aircraft.weights.mlw = x_in[1]              # Coupling variable
        aircraft.weights.mzfw = x_in[2]             # Coupling variable

        # Mass
        #------------------------------------------------------------------------------------------------------
        eval_mass_breakdown(aircraft)

        # Mission
        #------------------------------------------------------------------------------------------------------
        sub_proc.eval_nominal_mission(aircraft)

        y_out = np.array([aircraft.weights.mass_constraint_1,
                          aircraft.weights.mass_constraint_2,
                          aircraft.weights.mass_constraint_3])

        return y_out
    #-----------------------------------------------------------------------------------------------------------

    x_ini = np.array([aircraft.weights.mtow,
                      aircraft.weights.mlw,
                      aircraft.weights.mzfw])

    fct_arg = aircraft

    output_dict = fsolve(fct_mass_mission, x0=x_ini, args=fct_arg, full_output=True)

    aircraft.weights.mtow = output_dict[0][0]                         # Coupling variable
    aircraft.weights.mlw = output_dict[0][1]                          # Coupling variable
    aircraft.weights.mzfw = output_dict[0][2]                         # Coupling variable

    # Update mass data
    #------------------------------------------------------------------------------------------------------
    eval_mass_breakdown(aircraft)

    # Update mission data
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_nominal_mission(aircraft)

    return


#===========================================================================================================
def eval_payload_range_analysis(aircraft):
    """
    Compute Payload - Range diagram corner points
    """
    disa = 0
    altp = aircraft.design_driver.ref_cruise_altp
    mach = aircraft.design_driver.cruise_mach

    # Max payload mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    payload = aircraft.payload.maximum

    aircraft.max_payload_mission.tow = tow
    aircraft.max_payload_mission.payload = payload

    range,block_fuel,block_time,total_fuel = sub_proc.mission_range(aircraft,tow,payload,altp,mach,disa)

    aircraft.max_payload_mission.range = range
    aircraft.max_payload_mission.block_fuel = block_fuel
    aircraft.max_payload_mission.block_time = block_time
    aircraft.max_payload_mission.total_fuel = total_fuel

    # Max fuel mission
    #------------------------------------------------------------------------------------------------------
    tow = aircraft.weights.mtow
    total_fuel = aircraft.weights.mfw

    aircraft.max_fuel_mission.tow = tow
    aircraft.max_fuel_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.max_fuel_mission.payload = payload
    aircraft.max_fuel_mission.range = range
    aircraft.max_fuel_mission.block_fuel = block_fuel
    aircraft.max_fuel_mission.block_time = block_time

    # zero fuel mission
    #------------------------------------------------------------------------------------------------------
    total_fuel = aircraft.weights.mfw
    tow = aircraft.weights.owe + total_fuel

    aircraft.zero_payload_mission.tow = tow
    aircraft.zero_payload_mission.total_fuel = total_fuel

    range,payload,block_fuel,block_time = sub_proc.mission_fuel_limited(aircraft,tow,total_fuel,altp,mach,disa)

    aircraft.zero_payload_mission.range = range
    aircraft.zero_payload_mission.block_fuel = block_fuel
    aircraft.zero_payload_mission.block_time = block_time

    return


#===========================================================================================================
def eval_climb_performances(aircraft):
    """
    Compute climb performances
    """

    # Ceilings
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.design_driver.top_of_climb_altp
    oei_ceil_req = aircraft.low_speed.req_oei_altp

    vz_clb,vz_crz,oei_path,oei_mach = perfo.ceilings(aircraft,toc,oei_ceil_req)

    aircraft.low_speed.eff_oei_path = oei_path
    aircraft.high_speed.eff_vz_climb = vz_clb
    aircraft.high_speed.eff_vz_cruise = vz_crz

    aircraft.low_speed.perfo_constraint_3 = (oei_path - aircraft.low_speed.req_oei_path) / aircraft.low_speed.req_oei_path

    aircraft.high_speed.perfo_constraint_1 = vz_clb - aircraft.high_speed.req_vz_climb
    aircraft.high_speed.perfo_constraint_2 = vz_crz - aircraft.high_speed.req_vz_cruise

    # Time to climb to requested altitude
    #------------------------------------------------------------------------------------------------------
    toc = aircraft.high_speed.req_toc_altp
    disa = 0
    mass = aircraft.weights.mtow
    vcas1 = aircraft.high_speed.cas1_ttc
    vcas2 = aircraft.high_speed.cas2_ttc
    mach = aircraft.design_driver.cruise_mach

    ttc = perfo.time_to_climb(aircraft,toc,disa,mass,vcas1,vcas2,mach)

    aircraft.high_speed.eff_ttc = ttc

    aircraft.high_speed.perfo_constraint_3 = (aircraft.high_speed.req_ttc - ttc) / aircraft.high_speed.req_ttc

    return


#===========================================================================================================
def eval_performance_analysis(aircraft):
    """
    Compute operational performances
    """

    # Nominal mission
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_nominal_mission(aircraft)

    # Take off field length
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_take_off_performances(aircraft)

    # Approach speed
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_landing_performances(aircraft)

    # Climb performances
    #------------------------------------------------------------------------------------------------------
    eval_climb_performances(aircraft)

    # Environment
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_co2_metric(aircraft)

    # Cost mission
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    sub_proc.eval_cost_mission(aircraft)

    # Economics
    #------------------------------------------------------------------------------------------------------
    sub_proc.eval_economics(aircraft)

    return


#===========================================================================================================
def eval_handling_quality_analysis(aircraft):
    """
    Compute CG limits from handling qualities point of view
    """

    c_g = aircraft.center_of_gravity

    # Forward limit : trim landing
    #------------------------------------------------------------------------------------------------------
    altp = unit.m_ft(0)
    disa = 0
    nei = 0
    speed_mode = 1
    hld_conf = aircraft.aerodynamics.hld_conf_ld
    mass = c_g.max_fwd_mass

    cg_max_fwd_stall,speed,fn,aoa,ih,c_z,cx_trimmed = h_q.forward_cg_stall(aircraft,altp,disa,nei,hld_conf,speed_mode,mass)

    c_g.max_fwd_trim_cg = cg_max_fwd_stall         # Forward cg limit

    c_g.cg_constraint_1 = c_g.max_fwd_trim_cg - c_g.max_fwd_req_cg

    # Backward limit : static stability
    #------------------------------------------------------------------------------------------------------
    stability_margin = regul.static_stability_margin()

    cg_max_bwd_stab = h_q.backward_cg_stab(aircraft,stability_margin)

    c_g.max_bwd_stab_cg = cg_max_bwd_stab          # Backward cg limit

    c_g.cg_constraint_2 = c_g.max_bwd_req_cg - c_g.max_bwd_stab_cg

    # Vertical tail sizing
    #------------------------------------------------------------------------------------------------------

    h_q.vertical_tail_sizing(aircraft)

    c_g.cg_constraint_3 = c_g.max_bwd_oei_req_cg - c_g.max_bwd_oei_cg

    return


#===========================================================================================================
def eval_hq0(aircraft):
    """
    Perform hq based empennage sizing without updating characteristic masses MTOW, MLW & MZFW
    """

    aircraft.center_of_gravity.cg_range_optimization = 1    # Start HQ optimization mode

    #===========================================================================================================
    def fct_hq_optim(x_in,aircraft):

        c_g = aircraft.center_of_gravity

        aircraft.vertical_tail.lever_arm = x_in[0]
        aircraft.horizontal_tail.area = x_in[1]
        aircraft.vertical_tail.area = x_in[2]

        eval_mda0(aircraft)

        eval_handling_quality_analysis(aircraft)

        y_out = np.array([c_g.cg_constraint_1,
                          c_g.cg_constraint_2,
                          c_g.cg_constraint_3])
        return y_out
    #-----------------------------------------------------------------------------------------------------------

    x_ini = np.array([aircraft.vertical_tail.lever_arm,
                      aircraft.horizontal_tail.area * 1.05,     # Down scaling initial HTP area reduces the number of convergence problems ...
                      aircraft.vertical_tail.area])

    fct_arg = aircraft

    output_dict = fsolve(fct_hq_optim, x0=x_ini, args=fct_arg, full_output=True)

    if (output_dict[2]!=1):
        print(output_dict[3])
        raise Exception("Convergence problem in HQ optimization")

    aircraft.vertical_tail.lever_arm = output_dict[0][0]
    aircraft.horizontal_tail.area = output_dict[0][1]
    aircraft.vertical_tail.area = output_dict[0][2]

    eval_mda0(aircraft)

    eval_handling_quality_analysis(aircraft)

    return


#===========================================================================================================
def eval_mda0(aircraft):
    """
    Run the design sequence with statistical empennage sizing
    """

    eval_aircraft_pre_design(aircraft)  # Solves geometrical coupling

    eval_mass_breakdown(aircraft)       # Just mass analysis without any solving

    eval_performance_analysis(aircraft)

    return


#===========================================================================================================
def eval_mda1(aircraft):
    """
    Run the design sequence with statistical empennage sizing and local mass constraints satisfaction
    """

    eval_aircraft_pre_design(aircraft)  # Solves geometrical coupling

    eval_mass_estimation(aircraft)      # Solves mass coupling on MZFW and MLW

    eval_performance_analysis(aircraft)

    return


#===========================================================================================================
def eval_mda2(aircraft):
    """
    Run the design sequence with mass-mission adaptation and statistical empennage sizing
    """

    eval_aircraft_pre_design(aircraft)      # Solves geometrical coupling

    eval_mass_mission_adaptation(aircraft)  # Solves mass coupling on MZFW, MLW and MTOW

    eval_performance_analysis(aircraft)

    return


#===========================================================================================================
def eval_mda3(aircraft):
    """
    Run full MDA design with mass-mission adaptation and hq based empennage sizing
    """

    aircraft.center_of_gravity.cg_range_optimization = 1    # Start HQ optimization mode

    #===========================================================================================================
    def fct_hq_optim(x_in,aircraft):

        c_g = aircraft.center_of_gravity

        aircraft.vertical_tail.lever_arm = x_in[0]
        aircraft.horizontal_tail.area = x_in[1]
        aircraft.vertical_tail.area = x_in[2]

        eval_mda2(aircraft)

        eval_handling_quality_analysis(aircraft)

        y_out = np.array([c_g.cg_constraint_1,
                          c_g.cg_constraint_2,
                          c_g.cg_constraint_3])
        return y_out
    #-----------------------------------------------------------------------------------------------------------

    x_ini = np.array([aircraft.vertical_tail.lever_arm,
                      aircraft.horizontal_tail.area * 1.05,     # Down scaling initial HTP area reduces the number of convergence problems ...
                      aircraft.vertical_tail.area])

    fct_arg = aircraft

    output_dict = fsolve(fct_hq_optim, x0=x_ini, args=fct_arg, full_output=True)

    if (output_dict[2]!=1):
        raise Exception("Convergence problem in HQ optimization")

    aircraft.vertical_tail.lever_arm = output_dict[0][0]
    aircraft.horizontal_tail.area = output_dict[0][1]
    aircraft.vertical_tail.area = output_dict[0][2]

    eval_mda1(aircraft)
    eval_handling_quality_analysis(aircraft)

    return


#===========================================================================================================
def eval_optim_data(x_in,ac,crit_index,crit_ref):
    """
    Compute criterion and constraints
    """

    if (ac.propulsion.architecture==1):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==2):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==3):

        ac.turbofan_engine.reference_thrust = x_in[0]

    elif (ac.propulsion.architecture==4):

        ac.turboprop_engine.reference_thrust = x_in[0]

    else:
        raise Exception("propulsion.architecture index is out of range")

    ac.wing.area = x_in[1]

    # Run MDA
    #------------------------------------------------------------------------------------------------------

    eval_mda2(ac)

    # Constraints are violated if negative
    #------------------------------------------------------------------------------------------------------
    cst = np.zeros(6)

    cst[0] = ac.high_speed.perfo_constraint_1
    cst[1] = ac.high_speed.perfo_constraint_2
    cst[2] = ac.low_speed.perfo_constraint_3

    cst[3] = ac.high_speed.perfo_constraint_3
    cst[4] = ac.low_speed.perfo_constraint_1
    cst[5] = ac.low_speed.perfo_constraint_2

#    omag = [ 1. , 1., 1.e-5,  1.e2, 1.e2, 1.e-2 ]

#    cst /= omag

    # All criteria have to be minimized
    #------------------------------------------------------------------------------------------------------
    crt = np.zeros(5)

    crt[0] = ac.weights.mtow
    crt[1] = ac.cost_mission.block_fuel
    crt[2] = ac.environmental_impact.CO2_metric
    crt[3] = ac.economics.cash_operating_cost
    crt[4] = ac.economics.direct_operating_cost

    crit = crt[crit_index] * (20./crit_ref)

    return crit,cst


#===========================================================================================================
def eval_optim_cst(x_in,aircraft,crit_index,crit_ref):
    """
    Retrieve constraints
    """

    crit,cst = eval_optim_data(x_in,aircraft,crit_index,crit_ref)

    print("cst :",cst)

    return cst


#===========================================================================================================
def eval_optim_crt(x_in,aircraft,crit_index,crit_ref):
    """
    Retreve criteria
    """

    crit,cst = eval_optim_data(x_in,aircraft,crit_index,crit_ref)

    print("Design :",x_in)
    print("Crit :",crit)

    return crit


#===========================================================================================================
def mdf_process(aircraft,search_domain,criterion):
    """
    Compute criterion and constraints
    """

    if (aircraft.propulsion.architecture==1):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==2):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==3):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==4):
        start_value = (aircraft.turboprop_engine.reference_thrust,aircraft.wing.area)
    else:
        raise Exception("propulsion.architecture index is out of range")


    if (criterion=="MTOW"):
        crit_index = 0
    elif (criterion=="block_fuel"):
        crit_index = 1
    elif (criterion=="CO2_metric"):
        crit_index = 2
    elif (criterion=="COC"):
        crit_index = 3
    elif (criterion=="DOC"):
        crit_index = 4
    else:
        raise Exception("Criterion name is unknown")

    eval_mda0(aircraft)     # Initialization (compulsory only with mda3)

    crit_ref,cst_ref = eval_optim_data(start_value,aircraft,crit_index,1.)

    res = minimize(eval_optim_crt, start_value, args=(aircraft,crit_index,crit_ref,), method="trust-constr",
                   jac="3-point", hess=SR1(), hessp=None, bounds=search_domain, tol=1e-5,
                   constraints=NonlinearConstraint(fun=lambda x:eval_optim_cst(x,aircraft,crit_index,crit_ref),
                                                   lb=0., ub=np.inf, jac='3-point'),
                   options={'maxiter':500,'gtol': 1e-13})
    #              tol=None, callback=None,
    #              options={'grad': None, 'xtol': 1e-08, 'gtol': 1e-08, 'barrier_tol': 1e-08,
    #                       'sparse_jacobian': None, 'maxiter': 1000, 'verbose': 0,
    #                       'finite_diff_rel_step': None, 'initial_constr_penalty': 1.0,
    #                       'initial_tr_radius': 1.0, 'initial_barrier_parameter': 0.1,
    #                       'initial_barrier_tolerance': 0.1, 'factorization_method': None, 'disp': False})

    # res = minimize(eval_optim_crt, start_value, args=(aircraft,crit_index,crit_ref,), method="SLSQP", bounds=search_domain,
    #                constraints={"type":"ineq","fun":eval_optim_cst,"args":(aircraft,crit_index,crit_ref,)},
    #                jac="2-point",options={"maxiter":30,"ftol":0.0001,"eps":0.01},tol=1e-14)

    #res = minimize(eval_optim_crt, x_in, args=(aircraft,crit_index,crit_ref,), method="COBYLA", bounds=((110000,140000),(120,160)),
    #               constraints={"type":"ineq","fun":eval_optim_cst,"args":(aircraft,crit_index,crit_ref,)},
    #               options={"maxiter":100,"tol":0.1,"catol":0.0002,'rhobeg': 1.0})
    print(res)

    return res


#===========================================================================================================
def plot_mdf_process(aircraft,search_domain,criterion):
    """
    Compute criterion and constraints
    """

    if (aircraft.propulsion.architecture==1):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==2):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==3):
        start_value = (aircraft.turbofan_engine.reference_thrust,aircraft.wing.area)
    elif (aircraft.propulsion.architecture==4):
        start_value = (aircraft.turboprop_engine.reference_thrust,aircraft.wing.area)
    else:
        raise Exception("propulsion.architecture index is out of range")


    if (criterion=="MTOW"):
        crit_index = 0
    elif (criterion=="block_fuel"):
        crit_index = 1
    elif (criterion=="CO2_metric"):
        crit_index = 2
    elif (criterion=="COC"):
        crit_index = 3
    elif (criterion=="DOC"):
        crit_index = 4
    else:
        raise Exception("Criterion name is unknown")


    # res = minimize(eval_optim_crt, start_value, args=(aircraft,crit_index,), method="SLSQP", bounds=search_domain,
    #                constraints={"type":"ineq","fun":eval_optim_cst,"args":(aircraft,crit_index,)},
    #                jac="2-point",options={"maxiter":30,"ftol":0.0001,"eps":0.01})

    def obj_catch(xval):
        try:
            out=eval_optim_cst(xval,aircraft,crit_index)
        except:
            out=[0,0,0,0,0,0]
        return out

    pmot=np.linspace(search_domain[0][0],search_domain[0][1],21)
    sref=np.linspace(search_domain[1][0],search_domain[1][1],21)

    func_vals1=[obj_catch([pmoti,sref[10]]) for pmoti in pmot]
    func_vals2=[obj_catch([pmot[10],srefi]) for srefi in sref]

    from matplotlib import  pylab
    pylab.plt.plot(pmot,func_vals1,label="pmot obj")
    pylab.figure()
    pylab.plt.plot(sref,func_vals2,label="sref obj")
    pylab.legend()
    pylab.show()

