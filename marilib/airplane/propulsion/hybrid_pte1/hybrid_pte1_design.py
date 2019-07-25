#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy

from scipy.optimize import fsolve
from marilib.tools.math import lin_interp_1d

from marilib.earth import environment as earth

from marilib.aircraft_model.airplane import aerodynamics as airplane_aero

from marilib.airplane.propulsion import jet_models as jet

from marilib.airplane.propulsion.turbofan.turbofan_models import turbofan_thrust

from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_models import pte1_thrust

from marilib.airplane.propulsion.hybrid_pte1 import hybrid_pte1_models as hybrid


#===========================================================================================================
def eval_pte1_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    battery = aircraft.pte1_battery
    power_elec = aircraft.pte1_power_elec_chain
    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    low_speed = aircraft.low_speed

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    # Propulsion architecture design, definition of e-fan power in each fligh t phase
    #-----------------------------------------------------------------------------------------------------------

    # Initialisation
    crm = design_driver.cruise_mach
    toc = design_driver.top_of_climb_altp
    rca = design_driver.ref_cruise_altp
    roa = low_speed.req_oei_altp

    #                      MTO   MCN    MCL  MCR  FIR
    fd_disa = {"MTO":5.  , "MCN":0.   , "MCL":0. , "MCR":0. , "FID":0. }
    fd_altp = {"MTO":0.  , "MCN":roa  , "MCL":toc, "MCR":rca, "FID":rca}
    fd_mach = {"MTO":0.25, "MCN":crm/2, "MCL":crm, "MCR":crm, "FID":crm}
    fd_nei  = {"MTO":0.  , "MCN":1.   , "MCL":0. , "MCR":0. , "FID":0. }

    r_engine.flight_data = {"disa":fd_disa, "altp":fd_altp, "mach":fd_mach, "nei":fd_nei}

    e_fan_power = {"MTO":power_elec.mto,
                   "MCN":power_elec.mcn,
                   "MCL":power_elec.mcl,
                   "MCR":power_elec.mcr,
                   "FID":power_elec.fid}

    # Battery power feed is used in temporary phases only (take off and climb)
    power_factor = battery.power_feed * r_nacelle.controller_efficiency * r_nacelle.motor_efficiency
    battery_power_feed = {"MTO":power_factor,
                          "MCN":0.,
                          "MCL":power_factor,
                          "MCR":0.,
                          "FID":0.}

    e_power_ratio = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}
    e_shaft_power = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}

    for rating in propulsion.rating_code:
        (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(fd_altp[rating],fd_disa[rating])
        (fn,data) = turbofan_thrust(aircraft,Pamb,Tamb,fd_mach[rating],rating,fd_nei[rating])
        (fn_core,fn_fan0,fn0,shaft_power0) = data

        if e_fan_power[rating]>1:       # required eFan shaft power is given, turbofan shaft power ratio is deduced

            # Fraction of the turbofan shaft power dedicated to electric generation
            e_power_ratio[rating] =  ( (e_fan_power[rating] - battery_power_feed[rating] \
                                        )/ power_elec.overall_efficiency \
                                      )/((shaft_power0)*(engine.n_engine-fd_nei[rating]))

            # e-fan shaft power
            e_shaft_power[rating] = e_fan_power[rating]

        else:       # required turbofan shaft power ration is given, absolute shaft power is deduced

            # Shaft power dedicated to electric generator
            shaft_power2 = e_power_ratio[rating]*shaft_power0*(engine.n_engine-fd_nei[rating])

            # Fraction of the shaft power dedicated to the electric generation
            e_power_ratio[rating] = e_fan_power[rating]

            e_shaft_power[rating] =   shaft_power2*power_elec.overall_efficiency \
                                    + battery_power_feed[rating]

    # Storing results
    r_engine.n_engine = 1   # Only one electric fan at rear end of the fuselage

    power_elec.mto_e_power_ratio = e_power_ratio[MTO]
    power_elec.mcn_e_power_ratio = e_power_ratio[MCN]
    power_elec.mcl_e_power_ratio = e_power_ratio[MCL]
    power_elec.mcr_e_power_ratio = e_power_ratio[MCR]
    power_elec.fid_e_power_ratio = e_power_ratio[FID]

    r_engine.mto_e_shaft_power = e_shaft_power[MTO]
    r_engine.mcn_e_shaft_power = e_shaft_power[MCN]
    r_engine.mcl_e_shaft_power = e_shaft_power[MCL]
    r_engine.mcr_e_shaft_power = e_shaft_power[MCR]
    r_engine.fid_e_shaft_power = e_shaft_power[FID]

    # Engine performance update
    #-----------------------------------------------------------------------------------------------------------
    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(fd_altp[MTO],fd_disa[MTO])
    (fn,data) = turbofan_thrust(aircraft,Pamb,Tamb,fd_mach[MTO],MTO,fd_nei[MTO])
    (fn_core,fn_fan0,fn0,shaft_power0) = data

    shaft_power1 = (1.-e_power_ratio[MTO])*shaft_power0     # Shaft power dedicated to the fan at take off

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*fd_mach[MTO]

    fn_fan1 = nacelle.efficiency_prop*shaft_power1/Vair     # Effective fan thrust

    engine.kfn_off_take = (fn_core + fn_fan1)/fn0       # Thrust reduction due to power off take for the e-fan

    return


#===========================================================================================================
def eval_pte1_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    vtp = aircraft.vertical_tail
    wing = aircraft.wing

    propulsion = aircraft.propulsion

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    # Turbofan nacelles geometry adjustment
    #-----------------------------------------------------------------------------------------------------------
    nacWidth0 = 0.49*engine.bpr**0.67 + 4.8e-6*engine.reference_thrust      # Reference dimensions of the nacelle without power off take

    nacLength0 = 0.86*nacWidth0 + engine.bpr**0.37

    kSize = numpy.sqrt(engine.kfn_off_take)      # Diameter decrease due to max thrust decrease

    kSize_eff = (kSize + engine.core_width_ratio * (1.-kSize))      # Diameter decrease considering core is unchanged

    nacelle.width = nacWidth0*kSize_eff     # Real nacelle diameter assuming core section remains unchanged

    nacelle.length = nacLength0*kSize_eff   # Nacelle length is reduced according to the same factor

    knac = numpy.pi*nacelle.width*nacelle.length

    nacelle.net_wetted_area = knac*(1.48 - 0.0076*knac)*engine.n_engine

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1):

        if (engine.n_engine==2):

            nacelle.y_ext = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = - 0.5 * fuselage.height \
                            + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width

        elif (engine.n_engine==4):

            nacelle.y_int = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_int = wing.x_root + (nacelle.y_int-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_int = - 0.5 * fuselage.height \
                            + (nacelle.y_int - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width

            nacelle.y_ext = 2.0 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = - 0.5 * fuselage.height \
                            + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width
        else:
            raise Exception("engine.n_engine, number of engine not supported")

    elif (nacelle.attachment == 2):

        if (engine.n_engine==2):

            nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

            nacelle.x_ext = vtp.x_root - 0.5*nacelle.length

            nacelle.z_ext = 0.5 * fuselage.height

        else:
            raise Exception("engine.n_engine, number of engine not supported")

    else:
        raise Exception("nacelle.attachment, index is out of range")

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    dISA = 0.
    Altp = design_driver.ref_cruise_altp
    Mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(Altp,dISA)

    shaft_power = r_engine.mcr_e_shaft_power
    hub_width = 0.5     # Diameter of the e fan hub

    body_length = fuselage.length
    body_width = fuselage.width

    jet.rear_nacelle_design(r_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width,body_length,body_width)

    r_nacelle.x_axe = fuselage.length + 0.2*r_nacelle.width
    r_nacelle.y_axe = 0.
    r_nacelle.z_axe = 0.91*fuselage.height - 0.55*fuselage.height

    # Engine performance update
    #-----------------------------------------------------------------------------------------------------------
    fd = r_engine.flight_data

    e_fan_thrust = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}

    for rating in propulsion.rating_code:

        altp = fd.get("altp")[rating]
        disa = fd.get("disa")[rating]
        mach = fd.get("mach")[rating]
        nei = fd.get("nei")[rating]

        (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)
        (fn,sec,data) = pte1_thrust(aircraft,Pamb,Tamb,mach,rating,nei)
        (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0) = data

        e_fan_thrust[rating] = fn_fan2

    r_engine.mto_e_fan_thrust = e_fan_thrust[MTO]
    r_engine.mcn_e_fan_thrust = e_fan_thrust[MCN]
    r_engine.mcl_e_fan_thrust = e_fan_thrust[MCL]
    r_engine.mcr_e_fan_thrust = e_fan_thrust[MCR]
    r_engine.fid_e_fan_thrust = e_fan_thrust[FID]

    (eFanFnBli,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,Pamb,Tamb,Mach,shaft_power)

    (eFanFn,q0) = jet.fan_thrust(r_nacelle,Pamb,Tamb,Mach,shaft_power)

    propulsion.bli_e_thrust_factor = eFanFnBli / eFanFn     # Thrust increase due to BLI at iso shaft power for the e-fan

    propulsion.bli_thrust_factor = 1.     # Thrust increase due to BLI at iso shaft power for the turbofans (provision)

    return


#===========================================================================================================
def eval_pte1_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    power_elec = aircraft.pte1_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    e_shaft_power = numpy.array([r_engine.mto_e_shaft_power,
                                 r_engine.mcn_e_shaft_power,
                                 r_engine.mcl_e_shaft_power,
                                 r_engine.mcr_e_shaft_power,
                                 r_engine.fid_e_shaft_power])

    shaftPowerMax = max(e_shaft_power)

    turboFanMass0 = 1250. + 0.021*engine.reference_thrust # Statistical regression

    turboFanMass1 = 1250. + 0.021*engine.reference_thrust*engine.kfn_off_take

    kTurboFanMass = turboFanMass1 / turboFanMass0

    kMass = kTurboFanMass + engine.core_weight_ratio*(1-kTurboFanMass)     # Assuming core mass remains unchanged

    nacelle.mass = engine.n_engine * turboFanMass0 * kMass     # Total engine mass

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                       ) * shaftPowerMax

    r_nacelle.mass = (  1./r_nacelle.controller_pw_density + 1./r_nacelle.motor_pw_density \
                      + 1./r_nacelle.nacelle_pw_density \
                      ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    r_nacelle.c_g = fuselage.length + 0.5*r_nacelle.length

    return


#===========================================================================================================
def eval_pte1_battery_mass(aircraft):
    """
    Battery predesign
    """

    fuselage = aircraft.fuselage

    weights = aircraft.weights
    c_o_g = aircraft.center_of_gravity
    propulsion = aircraft.propulsion

    battery = aircraft.pte1_battery

    battery.c_g = fuselage.c_g

    if (battery.strategy==1):
        battery.mass = (battery.power_feed*battery.time_feed + battery.energy_cruise)/battery.energy_density
        propulsion.battery_energy_density = battery.energy_density
        weights.battery = battery.mass
        c_o_g.battery = battery.c_g

    elif (battery.strategy==2):

        battery.energy_cruise = max(0.,battery.mass*battery.energy_density - battery.power_feed*battery.time_feed)
        propulsion.battery_energy_density = battery.energy_density
        weights.battery = battery.mass
        c_o_g.battery = battery.c_g

    else:
        raise Exception("battery.strategy index is out of range")


    return


