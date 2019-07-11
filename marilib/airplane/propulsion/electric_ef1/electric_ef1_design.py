#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

import numpy

from scipy.optimize import fsolve

from marilib.earth import environment as earth

from marilib.airplane.propulsion import jet_models as jet


#===========================================================================================================
def eval_ef1_pylon_mass(aircraft):
    """
    Electrofan pylon mass & CG estimation
    """

    nacelle = aircraft.electrofan_nacelle
    engine = aircraft.electrofan_engine

    pylon = aircraft.electrofan_pylon

    pylon.mass = 0.0031*engine.reference_thrust*engine.n_engine

    if (engine.n_engine==2):
        pylon.c_g = nacelle.x_ext + 0.75*nacelle.length
    elif (engine.n_engine==4):
        pylon.c_g = 0.5*(nacelle.x_int + nacelle.x_ext) + 0.75*nacelle.length
    else:
        raise Exception("Number of engine is not allowed")

    return


#===========================================================================================================
def eval_ef1_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    propulsion = aircraft.propulsion

    power_elec = aircraft.ef1_power_elec_chain
    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    engine.rating_factor = {"MTO":1.00, "MCN":0.90, "MCL":0.90, "MCR":0.80, "FID":0.05}

    # Propulsion architecture design, definition of e-fan power in each flight phase
    #-----------------------------------------------------------------------------------------------------------
    disa = 15.
    altp = 0.
    mach = 0.25

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    Vsnd = earth.sound_speed(tamb)
    Vair = mach*Vsnd

    ref_shaft_power = engine.reference_thrust * Vair / nacelle.efficiency_prop

    engine.mto_e_shaft_power = ref_shaft_power * engine.rating_factor[MTO]
    engine.mcn_e_shaft_power = ref_shaft_power * engine.rating_factor[MCN]
    engine.mcl_e_shaft_power = ref_shaft_power * engine.rating_factor[MCL]
    engine.mcr_e_shaft_power = ref_shaft_power * engine.rating_factor[MCR]
    engine.fid_e_shaft_power = ref_shaft_power * engine.rating_factor[FID]

    e_shaft_power = numpy.array([engine.mto_e_shaft_power,
                                 engine.mcn_e_shaft_power,
                                 engine.mcl_e_shaft_power,
                                 engine.mcr_e_shaft_power,
                                 engine.fid_e_shaft_power])

    power_elec.max_power = max(e_shaft_power)
    power_elec.max_power_rating = numpy.argmax(e_shaft_power)

    return


#===========================================================================================================
def eval_ef1_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    propulsion = aircraft.propulsion

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    disa = 0.
    altp = design_driver.ref_cruise_altp
    mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)

    shaft_power = engine.mcr_e_shaft_power

    hub_width = 0.2

    eval_efan_nacelle_design(nacelle,Pamb,Tamb,mach,shaft_power,hub_width)

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1) :  # Nacelles are attached under the wing

        nacelle.y_ext = 0.7 * fuselage.width + 1.4 * nacelle.width      # statistical regression

        nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

        nacelle.z_ext = - 0.5 * fuselage.height \
                    + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                    - 0.5*nacelle.width

    elif (nacelle.attachment == 2) :    # Nacelles are attached on rear fuselage

        nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

        nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

        nacelle.z_ext = 0.5 * fuselage.height

    return


#===========================================================================================================
def eval_efan_nacelle_design(this_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width):
    """
    Electrofan nacelle design
    """

    gam = earth.heat_ratio()
    r = earth.gaz_constant()
    Cp = earth.heat_constant(gam,r)

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    # Electrical nacelle geometry : e-nacelle diameter is size by cruise conditions
    #-----------------------------------------------------------------------------------------------------------

    deltaV = 2.*Vair*(this_nacelle.efficiency_fan/this_nacelle.efficiency_prop - 1.)      # speed variation produced by the fan

    PwInput = this_nacelle.efficiency_fan*shaft_power     # kinetic energy produced by the fan

    Vinlet = Vair
    Vjet = Vinlet + deltaV

    q1 = 2.*PwInput / (Vjet**2 - Vinlet**2)

    MachInlet = Mach     # The inlet is in free stream

    Ptot = earth.total_pressure(Pamb,MachInlet)        # Stagnation pressure at inlet position

    Ttot = earth.total_temperature(Tamb,MachInlet)     # Stagnation temperature at inlet position

    MachFan = 0.5       # required Mach number at fan position

    CQoA1 = jet.corrected_air_flow(Ptot,Ttot,MachFan)        # Corrected air flow per area at fan position

    eFanArea = q1/CQoA1     # Fan area around the hub

    fan_width = numpy.sqrt(hub_width**2 + 4*eFanArea/numpy.pi)        # Fan diameter

    TtotJet = Ttot + shaft_power/(q1*Cp)        # Stagnation pressure increases due to introduced work

    Tstat = TtotJet - 0.5*Vjet**2/Cp        # static temperature

    VsndJet = numpy.sqrt(gam*r*Tstat) # Sound velocity at nozzle exhaust

    MachJet = Vjet/VsndJet # Mach number at nozzle output

    PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb)

    CQoA2 = jet.corrected_air_flow(PtotJet,TtotJet,MachJet)     # Corrected air flow per area at nozzle output

    nozzle_area = q1/CQoA2        # Fan area around the hub

    nozzle_width = numpy.sqrt(4*nozzle_area/numpy.pi)       # Nozzle diameter

    this_nacelle.hub_width = hub_width

    this_nacelle.fan_width = fan_width

    this_nacelle.nozzle_width = nozzle_width

    this_nacelle.nozzle_area = nozzle_area

    this_nacelle.width = 1.20*fan_width      # Surrounding structure

    this_nacelle.length = 1.50*this_nacelle.width

    this_nacelle.net_wetted_area = numpy.pi*this_nacelle.width*this_nacelle.length        # Nacelle wetted area

    return


#===========================================================================================================
def eval_ef1_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    power_elec = aircraft.pte1_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    shaftPowerMax = power_elec.max_power

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                      ) * shaftPowerMax

    nacelle.mass = (  1./nacelle.controller_pw_density + 1./nacelle.motor_pw_density \
                    + 1./nacelle.nacelle_pw_density \
                   ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    nacelle.c_g = fuselage.length + 0.5*nacelle.length

    return


#===========================================================================================================
def eval_wing_battery_data(aircraft):
    """
    Wing battery predesign using tank data structure
    """

    propulsion = aircraft.propulsion
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    battery = aircraft.ef1_battery
    tanks = aircraft.tanks

    tanks.cantilever_volume = 0.20 * (wing.area*wing.mac*(0.50*wing.t_o_c_r + 0.30*wing.t_o_c_k + 0.20*wing.t_o_c_t))

    tanks.central_volume = 1.3 * fuselage.width * wing.t_o_c_r * wing.mac**2

    tanks.fuel_density = earth.fuel_density(propulsion.fuel_type)

    tanks.mfw_volume_limited = (tanks.central_volume + tanks.cantilever_volume)*battery.density

    tanks.fuel_cantilever_cg =  0.25*(wing.x_root + 0.40*wing.c_root) \
                              + 0.65*(wing.x_kink + 0.40*wing.c_kink) \
                              + 0.10*(wing.x_tip + 0.40*wing.c_tip)

    tanks.fuel_central_cg = wing.x_root + 0.30*wing.c_root

    tanks.fuel_total_cg =  (tanks.fuel_cantilever_cg*tanks.cantilever_volume + tanks.fuel_central_cg*tanks.central_volume) \
                        / (tanks.central_volume + tanks.cantilever_volume)

    # Batteries will not change their mass during flight
    tanks.fuel_max_fwd_cg = tanks.fuel_total_cg
    tanks.fuel_max_fwd_mass = tanks.mfw_volume_limited * battery.fill_factor

    tanks.fuel_max_bwd_cg = tanks.fuel_total_cg
    tanks.fuel_max_bwd_mass = tanks.mfw_volume_limited * battery.fill_factor

    return


#===========================================================================================================
def eval_ef1_battery_mass(aircraft):
    """
    Wing battery predesign using tank data structure
    """

    tanks = aircraft.tanks

    battery = aircraft.ef1_battery

    battery.mass = tanks.mfw_volume_limited * battery.fill_factor
    battery.c_g = tanks.fuel_total_cg

    return


