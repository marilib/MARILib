#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

from marilib import numpy

from marilib import fsolve
from marilib.tools.math import lin_interp_1d

from marilib.earth import environment as earth

from marilib.aircraft_model.airplane import aerodynamics as airplane_aero

from marilib.airplane.propulsion import jet_models as jet

from marilib.airplane.propulsion.turbofan.turbofan_models import turbofan_thrust

from marilib.airplane.propulsion.hybrid_pte2.hybrid_pte2_models import pte2_thrust

from marilib.airplane.propulsion.hybrid_pte2 import hybrid_pte2_models as hybrid


#===========================================================================================================
def eval_blimp_body_design(aircraft):
    """
    Blimp bodies design
    """

    wing = aircraft.wing

    blimp_body = aircraft.pte2_blimp_body

    blimp_body.x_axe = wing.x_root - 0.5*blimp_body.length
    blimp_body.y_axe = 0.5*wing.y_tip
    blimp_body.z_axe = wing.z_root

    blimp_body.net_wetted_area = 2. * (2.70*blimp_body.length*blimp_body.width)

    blimp_structure_width = 0.20

    usable_section = numpy.pi*(0.5*blimp_body.width - blimp_structure_width)**2

    blimp_body.usable_volume = 2. * (0.85*blimp_body.length*usable_section)

    g = earth.gravity()
    r_air = earth.gas_constant("air")
    r_he = earth.gas_constant(blimp_body.gas_type)

    disa = 0.
    rca = aircraft.design_driver.ref_cruise_altp
    [pamb_ref,tamb_ref,tstd,dtodz] = earth.atmosphere(rca,disa)

    m_he = (pamb_ref*blimp_body.usable_volume) / (r_he*tamb_ref)

    blimp_body.gas_mass = m_he

    m_air = m_he*(r_he/r_air)

    # Air and internal gas are supposed at the same pressure and temperature
    blimp_body.buoyancy_force = g*(m_he - m_air)

    return


#===========================================================================================================
def eval_blimp_body_mass(aircraft):
    """
    Blimp bodies mass & CG
    """

    blimp_body = aircraft.pte2_blimp_body

    blimp_body.mass = 4.0 * blimp_body.net_wetted_area

    blimp_body.c_g = blimp_body.x_axe + 0.5*blimp_body.length

    return


#===========================================================================================================
def eval_pte2_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine

    power_elec = aircraft.pte2_power_elec_chain
    r_engine = aircraft.rear_electric_engine

    low_speed = aircraft.low_speed

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    engine.rating_factor = {"MTO":0.800, "MCN":0.688, "MCL":0.624, "MCR":0.560, "FID":0.100}

    # Propulsion architecture design, definition reference conditions for engine performances
    #-----------------------------------------------------------------------------------------------------------

    # Initialisation
    crm = design_driver.cruise_mach
    toc = design_driver.top_of_climb_altp
    rca = design_driver.ref_cruise_altp
    roa = low_speed.req_oei_altp

    #                      MTO   MCN    MCL  MCR  FID
    fd_disa = {"MTO":15. , "MCN":0.   , "MCL":0. , "MCR":0. , "FID":0. }
    fd_altp = {"MTO":0.  , "MCN":roa  , "MCL":toc, "MCR":rca, "FID":rca}
    fd_mach = {"MTO":0.25, "MCN":crm/2, "MCL":crm, "MCR":crm, "FID":crm}
    fd_nei  = {"MTO":0.  , "MCN":1.   , "MCL":0. , "MCR":0. , "FID":0. }

    propulsion.flight_data = {"disa":fd_disa, "altp":fd_altp, "mach":fd_mach, "nei":fd_nei}

    # Max rear fan shaft power
    #-----------------------------------------------------------------------------------------------------------
    r_shaft_power = numpy.array([r_engine.mto_r_shaft_power,
                                 r_engine.mcn_r_shaft_power,
                                 r_engine.mcl_r_shaft_power,
                                 r_engine.mcr_r_shaft_power,
                                 r_engine.fid_r_shaft_power])

    power_elec.max_power = max(r_shaft_power)
    power_elec.max_power_rating = numpy.argmax(r_shaft_power)

    return


#===========================================================================================================
def eval_pte2_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage

    propulsion = aircraft.propulsion

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    blimp_body = aircraft.pte2_blimp_body

    power_elec = aircraft.pte2_power_elec_chain
    r_engine = aircraft.rear_electric_engine

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    # Thrust factor at take off due to power off take
    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MTO]
    altp = propulsion.flight_data["altp"][MTO]
    mach = propulsion.flight_data["mach"][MTO]
    nei = propulsion.flight_data["nei"][MTO]

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    vsnd = earth.sound_speed(tamb)
    vair = vsnd*mach

    throttle = 1.

    fn,sfc0,data = turbofan_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    power_offtake = throttle * r_engine.mto_r_shaft_power/(nacelle.n_engine-nei) / power_elec.overall_efficiency

    shaft_power1 = shaft_power0 - power_offtake     # Shaft power dedicated to the fan

    fn_fan1 = nacelle.efficiency_prop*shaft_power1/vair     # Effective fan thrust

    engine.kfn_off_take = (fn_core + fn_fan1)/fn0       # Thrust reduction due to power off take for the e-fan

    # Turbofan nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MCR]
    altp = propulsion.flight_data["altp"][MCR]
    mach = propulsion.flight_data["mach"][MCR]
    nei = propulsion.flight_data["nei"][MCR]

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)

    throttle = 0.8  #  REMARK : Design thrust is arbitrary set to 80% of max cruise thrust

    fn,sfc0,data = turbofan_thrust(aircraft,pamb,tamb,mach,MCR,throttle,nei)
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    shaft_power1 = shaft_power0 - r_engine.mcr_r_shaft_power/(nacelle.n_engine - nei)    # Shaft power dedicated to the fan

    nacelle.hub_width = 0.5     # Diameter of the fan hub for blimp bodies

    body_hub_width = nacelle.hub_width

    body_length = blimp_body.length
    body_width = blimp_body.width

    jet.rear_nacelle_design(nacelle,pamb,tamb,mach,shaft_power1,body_hub_width,body_length,body_width)  # Nacelles are at rear end of blimb bodies

    nacelle.y_ext = blimp_body.y_axe
    nacelle.x_ext = blimp_body.x_axe + blimp_body.length
    nacelle.z_ext = blimp_body.z_axe

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    dISA = 0.
    Altp = design_driver.ref_cruise_altp
    Mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(Altp,dISA)

    shaft_power = r_engine.mcr_r_shaft_power
    hub_width = 0.5     # Diameter of the e fan hub

    body_length = fuselage.length
    body_width = fuselage.width

    jet.rear_nacelle_design(r_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width,body_length,body_width)

    r_nacelle.x_axe = fuselage.length + 0.2*r_nacelle.width
    r_nacelle.y_axe = 0.
    r_nacelle.z_axe = 0.91*fuselage.height - 0.55*fuselage.height

    # Rear fan max thrust on each rating
    #-----------------------------------------------------------------------------------------------------------
    r_fan_thrust = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}

    r_shaft_power = {"MTO":r_engine.mto_r_shaft_power,
                     "MCN":r_engine.mcn_r_shaft_power,
                     "MCL":r_engine.mcl_r_shaft_power,
                     "MCR":r_engine.mcr_r_shaft_power,
                     "FID":r_engine.fid_r_shaft_power}

    for rating in propulsion.rating_code:

        altp = propulsion.flight_data["altp"][rating]
        disa = propulsion.flight_data["disa"][rating]
        mach = propulsion.flight_data["mach"][rating]

        (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

        if (r_shaft_power[rating] > 0.):

            if (propulsion.bli_effect>0):
                (fn_fan2,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,pamb,tamb,mach,r_shaft_power[rating])
            else:
                (fn_fan2,q0) = jet.fan_thrust(r_nacelle,pamb,tamb,mach,r_shaft_power[rating])

        else:

            fn_fan2 = 0.

        r_fan_thrust[rating] = fn_fan2

    r_engine.mto_r_fan_thrust = r_fan_thrust[MTO]
    r_engine.mcn_r_fan_thrust = r_fan_thrust[MCN]
    r_engine.mcl_r_fan_thrust = r_fan_thrust[MCL]
    r_engine.mcr_r_fan_thrust = r_fan_thrust[MCR]
    r_engine.fid_r_fan_thrust = r_fan_thrust[FID]

    #-----------------------------------------------------------------------------------------------------------
    if (propulsion.bli_effect>0):
        (eFanFnBli,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,Pamb,Tamb,Mach,shaft_power)
        (eFanFn,q0) = jet.fan_thrust(r_nacelle,Pamb,Tamb,Mach,shaft_power)
        propulsion.bli_r_thrust_factor = eFanFnBli / eFanFn     # Thrust increase due to BLI at iso shaft power for the e-fan
    else:
        propulsion.bli_r_thrust_factor = 1.     # Thrust increase due to BLI at iso shaft power for the e-fan
    #-----------------------------------------------------------------------------------------------------------
    if (propulsion.bli_effect>0):
        (FanFnBli,q1,dVbli) = jet.fan_thrust_with_bli(nacelle,Pamb,Tamb,Mach,shaft_power)
        (FanFn,q0) = jet.fan_thrust(nacelle,Pamb,Tamb,Mach,shaft_power)
        propulsion.bli_thrust_factor = FanFnBli / FanFn     # Thrust increase due to BLI at iso shaft power for the turbofans
    else:
        propulsion.bli_thrust_factor = 1.     # Thrust increase due to BLI at iso shaft power for the e-fan

    return


#===========================================================================================================
def eval_pte2_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    power_elec = aircraft.pte2_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    shaftPowerMax = power_elec.max_power

    turbo_fan_mass0 = 1250. + 0.021*engine.reference_thrust # Statistical regression

    turbo_fan_mass1 = 1250. + 0.021*engine.reference_thrust*engine.kfn_off_take

    k_turbo_fan_mass = turbo_fan_mass1 / turbo_fan_mass0

    k_mass = k_turbo_fan_mass + engine.core_weight_ratio*(1-k_turbo_fan_mass)     # Assuming core mass remains unchanged

    nacelle.mass = nacelle.n_engine * turbo_fan_mass0 * k_mass     # Total engine mass

    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    r_nacelle.mass = (  1./r_nacelle.controller_pw_density + 1./r_nacelle.motor_pw_density \
                      + 1./r_nacelle.nacelle_pw_density \
                      ) * shaftPowerMax

    r_nacelle.c_g = fuselage.length + 0.5*r_nacelle.length

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                       ) * shaftPowerMax

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    return


#===========================================================================================================
def eval_pte2_battery_mass(aircraft):
    """
    Battery predesign
    """

    fuselage = aircraft.fuselage

    weights = aircraft.weights
    c_o_g = aircraft.center_of_gravity
    propulsion = aircraft.propulsion

    battery = aircraft.pte2_battery

    battery.c_g = fuselage.c_g

    if (battery.strategy==0):
        propulsion.battery_energy_density = 0.
        weights.battery_in_owe = 0.
        c_o_g.battery = 0.

    elif (battery.strategy==1):
        battery.mass = (battery.power_feed*battery.time_feed + battery.energy_cruise)/battery.energy_density
        propulsion.battery_energy_density = battery.energy_density
        weights.battery_in_owe = battery.mass
        c_o_g.battery = battery.c_g

    elif (battery.strategy==2):

        battery.energy_cruise = max(0.,battery.mass*battery.energy_density - battery.power_feed*battery.time_feed)
        propulsion.battery_energy_density = battery.energy_density
        weights.battery_in_owe = battery.mass
        c_o_g.battery = battery.c_g

    else:
        raise Exception("battery.strategy index is out of range")


    return


