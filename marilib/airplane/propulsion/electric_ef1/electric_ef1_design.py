#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib import numpy

from marilib import fsolve

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

    pylon.mass = 0.0031*engine.reference_thrust*nacelle.n_engine

    if (nacelle.n_engine==2):
        pylon.c_g = nacelle.x_ext + 0.75*nacelle.length
    elif (nacelle.n_engine==4):
        pylon.c_g = 0.5*(nacelle.x_int + nacelle.x_ext) + 0.75*nacelle.length
    else:
        raise Exception("Number of engine is not allowed")

    return


#===========================================================================================================
def eval_ef1_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    low_speed = aircraft.low_speed
    propulsion = aircraft.propulsion

    power_elec = aircraft.ef1_power_elec_chain
    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    engine.rating_factor = {"MTO":1.00, "MCN":0.80, "MCL":0.80, "MCR":0.80, "FID":0.05}

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

    # Propulsion architecture design, definition of e-fan power in each flight phase
    #-----------------------------------------------------------------------------------------------------------
    # Main engines
    #-----------------------------------------------------------------------------------------------------------
    disa = fd_disa[MTO]
    altp = fd_altp[MTO]
    mach = fd_mach[MTO]

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
    Vsnd = earth.sound_speed(tamb)
    Vair = mach*Vsnd

    # reference_thrust is driving the design through ref_shat_power but actual thrust is computed later
    ref_shaft_power = 0.8 * engine.reference_thrust * Vair / nacelle.efficiency_prop

    engine.reference_power = ref_shaft_power

    engine.mto_e_shaft_power = ref_shaft_power * engine.rating_factor[MTO]
    engine.mcn_e_shaft_power = ref_shaft_power * engine.rating_factor[MCN]
    engine.mcl_e_shaft_power = ref_shaft_power * engine.rating_factor[MCL]
    engine.mcr_e_shaft_power = ref_shaft_power * engine.rating_factor[MCR]
    engine.fid_e_shaft_power = ref_shaft_power * engine.rating_factor[FID]

    # Max main fan shaft power
    #-----------------------------------------------------------------------------------------------------------
    shaft_power_array = numpy.array([engine.mto_e_shaft_power,
                                     engine.mcn_e_shaft_power,
                                     engine.mcl_e_shaft_power,
                                     engine.mcr_e_shaft_power,
                                     engine.fid_e_shaft_power])

    power_elec.max_power = max(shaft_power_array)
    power_elec.max_power_rating = numpy.argmax(shaft_power_array)

    return


#===========================================================================================================
def eval_ef1_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing
    vtp = aircraft.vertical_tail

    propulsion = aircraft.propulsion

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    nacelle.n_engine = propulsion.n_engine

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    disa = propulsion.flight_data["disa"][MCR]
    altp = propulsion.flight_data["altp"][MCR]
    mach = propulsion.flight_data["mach"][MCR]

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)

    shaft_power = engine.mcr_e_shaft_power

    hub_width = 0.2

    jet.efan_nacelle_design(nacelle,Pamb,Tamb,mach,shaft_power,hub_width)

    tan_phi0 = 0.25*(wing.c_kink-wing.c_tip)/(wing.y_tip-wing.y_kink) + numpy.tan(wing.sweep)

    if (nacelle.attachment == 1):

        if (nacelle.n_engine==2):

            nacelle.y_ext = 0.8 * fuselage.width + 1.5 * nacelle.width      # statistical regression

            nacelle.x_ext = wing.x_root + (nacelle.y_ext-wing.y_root)*tan_phi0 - 0.7*nacelle.length

            nacelle.z_ext = - 0.5 * fuselage.height \
                            + (nacelle.y_ext - 0.5 * fuselage.width) * numpy.tan(wing.dihedral) \
                            - 0.5*nacelle.width

        elif (nacelle.n_engine==4):

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
            raise Exception("nacelle.n_engine, number of engine not supported")

    elif (nacelle.attachment == 2):

        if (nacelle.n_engine==2):

            nacelle.y_ext = 0.5 * fuselage.width + 0.6 * nacelle.width      # statistical regression

            nacelle.x_ext = vtp.x_root - 0.5*nacelle.length

            nacelle.z_ext = 0.5 * fuselage.height

        else:
            raise Exception("nacelle.n_engine, number of engine not supported")

    else:
        raise Exception("nacelle.attachment, index is out of range")

    # Main fan max thrust on each rating
    #-----------------------------------------------------------------------------------------------------------
    fan_thrust = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}

    shaft_power = {"MTO":engine.mto_e_shaft_power,
                   "MCN":engine.mcn_e_shaft_power,
                   "MCL":engine.mcl_e_shaft_power,
                   "MCR":engine.mcr_e_shaft_power,
                   "FID":engine.fid_e_shaft_power}

    for rating in propulsion.rating_code:

        altp = propulsion.flight_data["altp"][rating]
        disa = propulsion.flight_data["disa"][rating]
        mach = propulsion.flight_data["mach"][rating]

        (pamb,tamb,tstd,dtodz) = earth.atmosphere(altp,disa)

        (fn_fan2,q0) = jet.fan_thrust(nacelle,pamb,tamb,mach,shaft_power[rating])

        fan_thrust[rating] = fn_fan2

    engine.mto_e_fan_thrust = fan_thrust[MTO]
    engine.mcn_e_fan_thrust = fan_thrust[MCN]
    engine.mcl_e_fan_thrust = fan_thrust[MCL]
    engine.mcr_e_fan_thrust = fan_thrust[MCR]
    engine.fid_e_fan_thrust = fan_thrust[FID]

    # Eventual rear nacelle
    #-----------------------------------------------------------------------------------------------------------
    if (nacelle.rear_nacelle==1):

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

        (eFanFnBli,q1,dVbli) = jet.fan_thrust_with_bli(r_nacelle,Pamb,Tamb,Mach,shaft_power)

        (eFanFn,q0) = jet.fan_thrust(r_nacelle,Pamb,Tamb,Mach,shaft_power)

        propulsion.bli_r_thrust_factor = eFanFnBli / eFanFn     # Thrust increase due to BLI at iso shaft power for the e-fan

        propulsion.bli_thrust_factor = 1.     # Thrust increase due to BLI at iso shaft power for the turbofans (provision)

    return


#===========================================================================================================
def eval_ef1_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.electrofan_engine
    nacelle = aircraft.electrofan_nacelle

    r_engine = aircraft.rear_electric_engine
    r_nacelle = aircraft.rear_electric_nacelle

    power_elec = aircraft.ef1_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    if (nacelle.rear_nacelle==1):

        r_shaft_power = numpy.array([r_engine.mto_r_shaft_power,
                                     r_engine.mcn_r_shaft_power,
                                     r_engine.mcl_r_shaft_power,
                                     r_engine.mcr_r_shaft_power,
                                     r_engine.fid_r_shaft_power])

        r_shaft_power_max = max(r_shaft_power)

        r_nacelle.mass = (  1./r_nacelle.controller_pw_density + 1./r_nacelle.motor_pw_density \
                          + 1./r_nacelle.nacelle_pw_density \
                          ) * r_shaft_power_max

        r_nacelle.c_g = fuselage.length + 0.5*nacelle.length

    else:

        r_shaft_power_max = 0.
        r_nacelle.mass = 0.
        r_nacelle.c_g = 0.

    shaft_power_max = power_elec.max_power

    nacelle.mass = (  1./nacelle.controller_pw_density + 1./nacelle.motor_pw_density \
                    + 1./nacelle.nacelle_pw_density \
                   ) * shaft_power_max * nacelle.n_engine

    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                      ) * (shaft_power_max * nacelle.n_engine + r_shaft_power_max)

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    return


#===========================================================================================================
def eval_ef1_battery_mass(aircraft):
    """
    Nevertheless, for simplicity reason, battery CG is supposed constant
    """

    tanks = aircraft.tanks
    propulsion = aircraft.propulsion
    battery = aircraft.ef1_battery

    propulsion.battery_energy_density = battery.energy_density

    battery.mass_max = tanks.mfw_volume_limited
    battery.c_g = tanks.fuel_total_cg

    if (battery.stacking=="Variable"):
        aircraft.weights.battery_in_owe = 0.
    elif (battery.stacking=="Max"):
        aircraft.weights.battery_in_owe = battery.mass_max
    elif (battery.stacking=="Given"):
        pass
    else:
        raise Exception("ef1_battery.stacking, index is unknown")

    aircraft.center_of_gravity.battery = battery.c_g

    return


#===========================================================================================================
def eval_battery_cg_range(aircraft):
    """
    Wing battery predesign using tank data structure
    """

    weights = aircraft.weights
    tanks = aircraft.tanks
    battery = aircraft.ef1_battery

    # Need to take account of any possible battery loading
    # WARNING : in early design steps, it may occur that the resulting weight of the airplane would be higher than MTOW
    if (battery.stacking=="Variable"):

        tanks.fuel_max_fwd_cg = tanks.fuel_central_cg    # Battery max forward CG, central volume is forward only within backward swept wing
        tanks.fuel_max_fwd_mass = tanks.central_volume*tanks.fuel_density

        tanks.fuel_max_bwd_cg = tanks.fuel_cantilever_cg    # Battery max Backward CG
        tanks.fuel_max_bwd_mass = tanks.cantilever_volume*tanks.fuel_density

    elif (battery.stacking=="Max"):

        tanks.fuel_max_fwd_cg = tanks.fuel_total_cg    # Battery max forward CG, central volume is forward only within backward swept wing
        tanks.fuel_max_fwd_mass = tanks.mfw_volume_limited*tanks.fuel_density

        tanks.fuel_max_bwd_cg = tanks.fuel_total_cg    # Battery max Backward CG
        tanks.fuel_max_bwd_mass = tanks.mfw_volume_limited*tanks.fuel_density

    elif (battery.stacking=="Given"):

        tanks.fuel_max_fwd_cg = tanks.fuel_total_cg    # Battery max forward CG, central volume is forward only within backward swept wing
        tanks.fuel_max_fwd_mass = weights.battery_in_owe

        tanks.fuel_max_bwd_cg = tanks.fuel_total_cg    # Battery max Backward CG
        tanks.fuel_max_bwd_mass = weights.battery_in_owe

    else:
        raise Exception("ef1_battery.stacking, index is unknown")



    return


