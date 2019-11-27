#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

from marilib import numpy

from marilib.tools import units as unit

from marilib.earth import environment as earth

#===========================================================================================================
def ref_cruise_altp(propulsive_architecture):
    if (propulsive_architecture=="TF"):
        ref_cruise_altp_i = unit.m_ft(35000.)
    elif (propulsive_architecture=="TP"):
        ref_cruise_altp_i = unit.m_ft(20000.)
    elif (propulsive_architecture=="PTE1"):
        ref_cruise_altp_i = unit.m_ft(35000.)
    elif (propulsive_architecture=="EF1"):
        ref_cruise_altp_i = unit.m_ft(25000.)
    else:
        raise Exception("propulsion.architecture index is out of range")
    return ref_cruise_altp_i

#===========================================================================================================
def top_of_climb_altp(propulsive_architecture):
    if (propulsive_architecture=="TF"):
        top_of_climb_altp_i = unit.m_ft(31000.)
    elif (propulsive_architecture=="TP"):
        top_of_climb_altp_i = unit.m_ft(15000.)
    elif (propulsive_architecture=="PTE1"):
        top_of_climb_altp_i = unit.m_ft(31000.)
    elif (propulsive_architecture=="EF1"):
        top_of_climb_altp_i = unit.m_ft(25000.)
    else:
        raise Exception("propulsion.architecture index is out of range")
    return top_of_climb_altp_i


#===========================================================================================================
def n_pax_front(n_pax_ref):
    if  (n_pax_ref<=8):   n_pax_front_i = 2
    elif(n_pax_ref<=16):  n_pax_front_i = 3
    elif(n_pax_ref<=70):  n_pax_front_i = 4
    elif(n_pax_ref<=120): n_pax_front_i = 5
    elif(n_pax_ref<=225): n_pax_front_i = 6
    elif(n_pax_ref<=300): n_pax_front_i = 8
    elif(n_pax_ref<=375): n_pax_front_i = 9
    else:                 n_pax_front_i = 10
    return n_pax_front_i

#===========================================================================================================
def n_aisle(n_pax_front):
    if(n_pax_front <= 6): n_aisle_i = 1
    else:                 n_aisle_i = 2
    return n_aisle_i

#===========================================================================================================
def fuselage_width(n_pax_front,n_aisle):
    fuselage_width_i = 0.38*n_pax_front + 1.05*n_aisle + 0.55
    return fuselage_width_i

#===========================================================================================================
def m_pax_nominal(design_range):
    if(design_range <= unit.m_NM(1500.)):
        m_pax_nominal_i = 95.
    elif(design_range <= unit.m_NM(3500.)):
        m_pax_nominal_i = 100.
    elif(design_range <= unit.m_NM(5500.)):
        m_pax_nominal_i = 105.
    else:
        m_pax_nominal_i = 110.
    return m_pax_nominal_i

#===========================================================================================================
def m_pax_max(design_range):
    if(design_range <= unit.m_NM(3500.)):
        m_pax_max_i = 120.
    elif(design_range <= unit.m_NM(5500.)):
        m_pax_max_i = 135.
    else:
        m_pax_max_i = 150.
    return m_pax_max_i

#===========================================================================================================
def cg_range_optimization():
    cg_range_optimization_i = 0       # Wing position, HTP area and VTP area optimized according to HQ criteria, 0: no, 1:yes
    return cg_range_optimization_i

#===========================================================================================================
def wing_attachment(propulsive_architecture):
    """
    wing_attachment = 1 : low wing
    wing_attachment = 2 : high wing
    """
    if (propulsive_architecture=="TP"):
        wing_attachment_i = 2
    else:
        wing_attachment_i = 1
    return wing_attachment_i

#===========================================================================================================
def wing_morphing():
    """
    wing_morphing = 1 : aspect ratio is the driver
    wing_morphing = 2 : wing span is the driver
    """
    wing_morphing_i = 1
    return wing_morphing_i

#===========================================================================================================
def hld_type(propulsive_architecture,n_pax_ref):
    """
    hld_type = 0  : Clean
    hld_type = 1  : Flap only, Rotation without slot
    hld_type = 2  : Flap only, Rotation with slot      (ATR)
    hld_type = 3  : Flap only, Rotation double slot
    hld_type = 4  : Flap only, Fowler
    hld_type = 5  : Slat only
    hld_type = 6  : Slat + Flap rotation double slot
    hld_type = 7  : Slat + Flap rotation with slot
    hld_type = 8  : Slat + Flap rotation double slot
    hld_type = 9  : Slat + Fowler                      (A320)
    hld_type = 10 : Slat + Slotted Fowler (A321)
    """
    if (n_pax_ref>100):
        hld_type_i = 9
    else:
        hld_type_i = 7
    return hld_type_i

#===========================================================================================================
def wing_area(n_pax_ref,design_range):
    wing_area_i = 60. + 88.*n_pax_ref*design_range*1.e-9
    return wing_area_i

#===========================================================================================================
def blimp_body_length(wing_span):
    blimp_body_length_i = wing_span
    return blimp_body_length_i

#===========================================================================================================
def blimp_body_width(wing_span):
    blimp_body_width_i = 0.18*wing_span
    return blimp_body_width_i

#===========================================================================================================
def blimp_body_gas():
    blimp_body_gas_i = "helium"
    return blimp_body_gas_i

#===========================================================================================================
def htp_area(wing_area):
    htp_area_i = 0.33*wing_area
    return htp_area_i

#===========================================================================================================
def vtp_area(wing_area):
    vtp_area_i = 0.2*wing_area
    return vtp_area_i

#===========================================================================================================
def hld_conf_clean():
    hld_conf_clean_i = 0        # By definition (0=<hld_conf=<1)
    return hld_conf_clean_i

#===========================================================================================================
def hld_conf_to():
    hld_conf_to_i = 0.3       # by definition (0=<hld_conf=<1)
    return hld_conf_to_i

#===========================================================================================================
def hld_conf_ld():
    hld_conf_ld_i = 1       # by definition (0=<hld_conf=<1)
    return hld_conf_ld_i

#===========================================================================================================
def wing_aspect_ratio(propulsive_architecture):
    wing_aspect_ratio_i = 9
    return wing_aspect_ratio_i

#===========================================================================================================
def wing_span(wing_area,wing_aspect_ratio):
    wing_span_i = numpy.sqrt(wing_area*wing_aspect_ratio)
    return wing_span_i

#===========================================================================================================
def wing_x_root(wing_aspect_ratio,wing_span):
    wing_x_root_i = 0.3*wing_span*numpy.sqrt(9/wing_aspect_ratio)
    return wing_x_root_i

#===========================================================================================================
def wing_sweep(cruise_mach):
    wing_sweep_i = 1.6*max(0,(cruise_mach-0.5))     # Empirical law
    return wing_sweep_i

#===========================================================================================================
def wing_mac(wing_area_i,wing_aspect_ratio_i):
    wing_mac_i = 1.2*numpy.sqrt(wing_area_i / wing_aspect_ratio_i)
    return wing_mac_i


#===========================================================================================================
def fuel_type(propulsive_architecture):
    if (propulsive_architecture=="EF1"):
        fuel_type_i = "Battery"     # Kerosene, Hydrogene, Methane, Battery
    else:
        fuel_type_i = "Kerosene"     # Kerosene, Hydrogene, Methane, Battery
    return fuel_type_i

#===========================================================================================================
def nacelle_body_length():
    nacelle_body_length_i = 4.
    return nacelle_body_length_i

#===========================================================================================================
def nacelle_body_width():
    nacelle_body_width_i = 1.5
    return nacelle_body_width_i

#===========================================================================================================
def nacelle_body_hub_width():
    nacelle_body_hub_width_i = 0.5
    return nacelle_body_hub_width_i

#===========================================================================================================
def ef1_rear_nacelle():
    efi_rear_nacelle_i = 0
    return efi_rear_nacelle_i

#===========================================================================================================
def nacelle_attachment(propulsive_architecture,n_pax_ref):
    if (propulsive_architecture=="TP"):
        nacelle_attachment_i = 1     # 1: underwing
    else:
        if (80<n_pax_ref):
            nacelle_attachment_i = 1     # 1: underwing
        else:
            nacelle_attachment_i = 2     # 2: on rear fuselage
    return nacelle_attachment_i

#===========================================================================================================
def rating_code():
    rating_code_i = ("MTO","MCN","MCL","MCR","FID")     # Engine rating codes
    return rating_code_i

#===========================================================================================================
def rating_factor(propulsive_architecture):
    if (propulsive_architecture=="TF"):
        rating_factor_i = {"MTO":0.800, "MCN":0.688, "MCL":0.624, "MCR":0.560, "FID":0.100}
    elif (propulsive_architecture=="EF1"):
        rating_factor_i = {"MTO":1.00, "MCN":0.80, "MCL":0.80, "MCR":0.80, "FID":0.05}
    else:
        rating_factor_i = None
    return rating_factor_i

#===========================================================================================================
def efficiency_fan():
    efficiency_fan_i = 0.95     # efficiency to convert shaft power into kinetic energy
    return efficiency_fan_i

#===========================================================================================================
def efficiency_prop():
    efficiency_prop_i = 0.82     # efficiency of a fan to convert shaft power into propulsive power
    return efficiency_prop_i

	
#===========================================================================================================
def propeller_efficiency():
    propeller_efficiency_i = 0.85     # efficiency of a propeller to convert shaft power into propulsive power
    return propeller_efficiency_i


#===========================================================================================================
def n_engine():
    n_engine_i = 2
    return n_engine_i

#===========================================================================================================
def bpr(n_pax_ref,propulsive_architecture):
    if (propulsive_architecture=="TF" or \
        propulsive_architecture=="PTE1" or \
        propulsive_architecture=="EF1"):
        if (80<n_pax_ref):
            bpr_i = 9.
        else:
            bpr_i = 5.
    elif (propulsive_architecture=="TP"):
        bpr_i = 30.
    else:
        raise Exception("propulsive architecture is not supported")
    return bpr_i

#===========================================================================================================
def reference_thrust(n_pax_ref,design_range,n_engine):
    reference_thrust_i = (1.e5 + 177.*n_pax_ref*design_range*1.e-6)/n_engine
    return reference_thrust_i

#===========================================================================================================
def nacelle_width(bpr,reference_thrust):
    turbofan_nacelle_width_i = (9.e-6*bpr*reference_thrust)**0.4
    return turbofan_nacelle_width_i

#===========================================================================================================
def propeller_width(reference_thrust):
    propeller_width_i = numpy.sqrt((4./numpy.pi)*(reference_thrust/3000.))      # Assuming 3000 N/m2
    return propeller_width_i

#===========================================================================================================
def nacelle_y_ext(propulsive_architecture,n_engine,attachment,fuselage_width_i,width_i):
    if (propulsive_architecture=="TP"):
        if (n_engine==2):
            nacelle_y_ext_i = 0.5 * fuselage_width_i + 0.7 * width_i
        elif (n_engine==4):
            nacelle_y_ext_i = 0.5 * fuselage_width_i + 1.9 * width_i
    else:
        if (attachment==1):
            if (n_engine==2):
                nacelle_y_ext_i = 0.7 * fuselage_width_i + 1.5 * width_i
            elif (n_engine==4):
                nacelle_y_ext_i = 1.9 * fuselage_width_i + 1.5 * width_i
        else:
            nacelle_y_ext_i = 0.5 * fuselage_width_i + 0.5 * width_i
    return nacelle_y_ext_i

#===========================================================================================================
def core_thrust_ratio():
    core_thrust_ratio_i = 0.13    # Mean contribution of the core to the total thrust
    return core_thrust_ratio_i

#===========================================================================================================
def core_width_ratio():
    core_width_ratio_i = 0.70    # Ratio of core mean diameter over fan diameter
    return core_width_ratio_i

#===========================================================================================================
def core_weight_ratio():
    core_weight_ratio_i = 0.13    # Mean contribution of the core to the total engine weight
    return core_weight_ratio_i

	
#===========================================================================================================
def prop_architecture():
    prop_architecture_i = "TF"    # prop_architecture, TF: turbofan, PTE1: partial turbo electric 1
    return prop_architecture_i

#===========================================================================================================
def tank_architecture():
    tank_architecture_i = 0    # tank_architecture, 0: wing box, 1:wing pods, 2:piggyback pod
    return tank_architecture_i

#===========================================================================================================
def tank_pod_length(tanks_architecture):
    if (tanks_architecture==0):
        tank_pod_length_i = 1.e-6
    elif (tanks_architecture==1):
        tank_pod_length_i = 6.
    elif (tanks_architecture==2):
        tank_pod_length_i = 10.
    return tank_pod_length_i

#===========================================================================================================
def tank_pod_width(tanks_architecture):
    if (tanks_architecture==0):
        tank_pod_width_i = 1.e-6
    elif (tanks_architecture==1):
        tank_pod_width_i = 2.
    elif (tanks_architecture==2):
        tank_pod_width_i = 2.
    return tank_pod_width_i

#===========================================================================================================
def pod_surface_mass():
    pod_surface_mass_i = 50.    # kg/m2
    return pod_surface_mass_i

#===========================================================================================================
def electric_shaft_power():
    electric_shaft_power_i = 1.e6    # Watts, electric motor power
    return electric_shaft_power_i

#===========================================================================================================
def idle_electric_shaft_power():
    idle_electric_shaft_power_i = 1.e4    # Watts, electric motor power
    return idle_electric_shaft_power_i

#===========================================================================================================
def battery_strategy():
    battery_strategy_i = 0    # Battery sizing strategy, 0: no battery, 1: power_feed & energy_cruise driven, 2: battery mass driven
    return battery_strategy_i

#===========================================================================================================
def battery_power_feed():
    battery_power_feed_i = 0.    # Power delivered to e-fan(s) at take off and(or) climb during a total of time_feed
    return battery_power_feed_i

#===========================================================================================================
def battery_time_feed():
    battery_time_feed_i = unit.s_min(15.)    # Maximum duration of the power_feed delivered to e-fan(s)
    return battery_time_feed_i

#===========================================================================================================
def battery_energy_cruise():
    battery_energy_cruise_i = 0.    # Total battery energy dedicated to cruise
    return battery_energy_cruise_i

#===========================================================================================================
def battery_energy_density(propulsive_architecture):
    if (propulsive_architecture in ("PTE1", "EF1")):
        battery_energy_density_i = unit.J_kWh(0.5)    # Battery energy density (kWh/kg)
    else:
        battery_energy_density_i = 0. # TF propulsion has no battery
    return battery_energy_density_i

#===========================================================================================================
def battery_power_density():
    battery_power_density_i = 1.e3    # Battery power density (capability to release power per mass unit (W/kg)
    return battery_power_density_i

#===========================================================================================================
def battery_density():
    battery_density_i = earth.fuel_density("Battery")    # Battery density (kg/m3)
    return battery_density_i

#===========================================================================================================
def e_chain_efficiency():
    e_chain_efficiency_i = 0.90    # Overall efficiency of the electric chain, from generator to motor
    return e_chain_efficiency_i

#===========================================================================================================
def controller_efficiency():
    controller_efficiency_i = 0.98    # no_dim
    return controller_efficiency_i

#===========================================================================================================
def e_motor_efficiency():
    motor_efficiency_i = 0.98    # no_dim
    return motor_efficiency_i

#===========================================================================================================
def generator_power_density():
    generator_power_density_i = 10.e3    # W/kg, Electric generator
    return generator_power_density_i

#===========================================================================================================
def rectifier_pw_density():
    rectifier_pw_density_i = 20.e3    # W/kg, Rectifier
    return rectifier_pw_density_i

#===========================================================================================================
def wiring_pw_density():
    wiring_pw_density_i = 20.e3    # W/kg, Wiring
    return wiring_pw_density_i

#===========================================================================================================
def cooling_pw_density():
    cooling_pw_density_i = 15.e3    # W/kg, Cooling
    return cooling_pw_density_i

#===========================================================================================================
def controller_pw_density():
    controller_pw_density_i = 20.e3    # W/kg, Electric motor
    return controller_pw_density_i

#===========================================================================================================
def e_motor_pw_density():
    e_motor_pw_density_i = 10.e3    # W/kg, Electric motor
    return e_motor_pw_density_i

#===========================================================================================================
def e_nacelle_pw_density():
    e_nacelle_pw_density_i = 5.e3    # W/kg, Electric nacelle
    return e_nacelle_pw_density_i

#===========================================================================================================
def boundary_layer_effect():
    bli_effect_i = 1    # 0: without, 1: with
    return bli_effect_i


#===========================================================================================================
def htp_attachment(propulsive_architecture,nacelle_attachment_i):
    if (propulsive_architecture=="TP"):
        htp_attachment_i = 2        # 2: T-tail (on top of fin)
    else:
        if (nacelle_attachment_i == 1):
            htp_attachment_i = 1        # 1: Classical (on fuselage tail cone)
        else:
            htp_attachment_i = 2        # 2: T-tail (on top of fin)
    return htp_attachment_i


#===========================================================================================================
def vtp_attachment():
    vtp_attachment_i = 1        # 1: Classical (on fuselage tail cone)
    return vtp_attachment_i


#===========================================================================================================
def mtow(n_pax_ref,design_range):
    mtow_i =  20500. + 67.e-6*n_pax_ref*design_range
    return mtow_i

#===========================================================================================================
def mzfw(n_pax_ref,design_range):
    mzfw_i = 25000. + 41.e-6*n_pax_ref*design_range
    return mzfw_i

#===========================================================================================================
def mlw(n_pax_ref,mtow_i,mzfw_i):
    if (n_pax_ref>100):
        mlw_i = min(mtow_i , (1.07*mzfw_i))
    else:
        mlw_i = mtow_i
    return mlw_i

#===========================================================================================================
def battery_mass():
    battery_mass_i = 0.
    return battery_mass_i

#===========================================================================================================
def battery_stacking():
    battery_stacking_i = "Variable"
    return battery_stacking_i


#===========================================================================================================
def disa_oei():
    init_disa_oei_i = 15.
    return init_disa_oei_i

#===========================================================================================================
def req_oei_altp(propulsive_architecture):
    req_oei_altp_i = unit.m_ft(11000.)
    return req_oei_altp_i

#===========================================================================================================
def oei_best_speed(cruise_mach):
    oei_best_speed_i = cruise_mach - 0.01
    return oei_best_speed_i


#===========================================================================================================
def altp_tofl():
    altp_tofl_i = 0.
    return altp_tofl_i

#===========================================================================================================
def disa_tofl():
    disa_tofl_i = 15.
    return disa_tofl_i

#===========================================================================================================
def req_tofl(design_range):
    if(design_range <= unit.m_NM(1500.)):
        req_tofl_i = 1500.
    elif(design_range <= unit.m_NM(3500.)):
        req_tofl_i = 2000.
    elif(design_range <= unit.m_NM(5500.)):
        req_tofl_i = 2500.
    else:
        req_tofl_i = 3000.
    return req_tofl_i


#===========================================================================================================
def altp_app_speed():
    altp_app_speed_i = unit.m_ft(0.)
    return altp_app_speed_i

#===========================================================================================================
def disa_app_speed():
    disa_app_speed_i = 0.
    return disa_app_speed_i

#===========================================================================================================
def req_app_speed(n_pax_ref):
    if (n_pax_ref<=100):
        req_app_speed_i = unit.mps_kt(120.)
    elif (n_pax_ref<=200):
        req_app_speed_i = unit.mps_kt(137.)
    else:
        req_app_speed_i = unit.mps_kt(140.)
    return req_app_speed_i


#===========================================================================================================
def cas1_ttc(cruise_mach):
    if (cruise_mach>=0.6):
        cas1_ttc_i = unit.mps_kt(250.)
    elif (cruise_mach>=0.4):
        cas1_ttc_i = unit.mps_kt(180.)
    else:
        cas1_ttc_i = unit.mps_kt(70.)
    return cas1_ttc_i

#===========================================================================================================
def cas2_ttc(cruise_mach):
    if (cruise_mach>=0.6):
        cas2_ttc_i = unit.mps_kt(300.)
    elif (cruise_mach>=0.4):
        cas2_ttc_i = unit.mps_kt(200.)
    else:
        cas2_ttc_i = unit.mps_kt(70.)
    return cas2_ttc_i

#===========================================================================================================
def req_ttc():
    req_ttc_i = unit.s_min(25.)
    return req_ttc_i


#===========================================================================================================
def disa_climb():
    disa_climb_i = 15.
    return disa_climb_i

#===========================================================================================================
def req_vz_climb():
    req_vz_climb_i = unit.mps_ftpmin(300.)
    return req_vz_climb_i

#===========================================================================================================
def req_vz_cruise():
    req_vz_cruise_i = unit.mps_ftpmin(0.)
    return req_vz_cruise_i

	
#===========================================================================================================
def cost_mission_disa():
    cost_mission_disa_i = 0.
    return cost_mission_disa_i

#===========================================================================================================
def cost_mission_range(design_range):
    if(design_range < unit.m_NM(400.)): cost_mission_range_i = unit.m_NM(100.)
    elif(design_range < unit.m_NM(1000.)): cost_mission_range_i = unit.m_NM(200.)
    elif(design_range < unit.m_NM(2500.)): cost_mission_range_i = unit.m_NM(400.)
    elif(design_range < unit.m_NM(4500.)): cost_mission_range_i = unit.m_NM(800.)
    elif(design_range < unit.m_NM(6500.)): cost_mission_range_i = unit.m_NM(2000.)
    else:                               cost_mission_range_i = unit.m_NM(4000.)
    return cost_mission_range_i

def fuel_price():
#===========================================================================================================
    fuel_price_i = 2./unit.liter_usgal(1)   # 2 $/USgal
    return fuel_price_i

def elec_price():
#===========================================================================================================
    elec_price_i = 0.15/unit.J_kWh(1)   # 0.15 $/kWh  Assumed de-carbonated
    return elec_price_i

def battery_mass_price():
#===========================================================================================================
    battery_mass_price_i = 20.   # $/kg    installed
    return battery_mass_price_i

#===========================================================================================================
def labor_cost():
    labor_cost_i = 120.   # 120 $/h
    return labor_cost_i

#===========================================================================================================
def irp():
    irp_i = 10.     # 10 years
    return irp_i

#===========================================================================================================
def period():
    period_i = 15.     # 15 years
    return period_i

#===========================================================================================================
def interest_rate():
    interest_rate_i = 0.04     # 4%
    return interest_rate_i

#===========================================================================================================
def utilisation(design_range):
    if(design_range <= unit.m_NM(3500.)): utilisation_i = 1600.
    else:                                utilisation_i = 600.
    return utilisation_i

