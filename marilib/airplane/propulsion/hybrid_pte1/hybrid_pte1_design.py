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

from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_models import hybrid_thrust

from marilib.airplane.propulsion.hybrid_pte1 import hybrid_pte1_models as hybrid


#===========================================================================================================
def eval_hybrid_engine_design(aircraft):
    """
    Thermal propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage

    propulsion = aircraft.propulsion
    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    battery = aircraft.battery
    power_elec = aircraft.pte1_power_elec_chain
    e_engine = aircraft.rear_electric_engine
    e_nacelle = aircraft.rear_electric_nacelle

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

    e_engine.flight_data = {"disa":fd_disa, "altp":fd_altp, "mach":fd_mach, "nei":fd_nei}

    e_fan_power = {"MTO":power_elec.mto,
                   "MCN":power_elec.mcn,
                   "MCL":power_elec.mcl,
                   "MCR":power_elec.mcr,
                   "FID":power_elec.fid}

    # Battery power feed is used in temporary phases only (take off and climb)
    power_factor = battery.power_feed * e_nacelle.controller_efficiency * e_nacelle.motor_efficiency
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
    e_engine.n_engine = 1   # Only one electric fan at rear end of the fuselage

    power_elec.mto_e_power_ratio = e_power_ratio[MTO]
    power_elec.mcn_e_power_ratio = e_power_ratio[MCN]
    power_elec.mcl_e_power_ratio = e_power_ratio[MCL]
    power_elec.mcr_e_power_ratio = e_power_ratio[MCR]
    power_elec.fid_e_power_ratio = e_power_ratio[FID]

    e_engine.mto_e_shaft_power = e_shaft_power[MTO]
    e_engine.mcn_e_shaft_power = e_shaft_power[MCN]
    e_engine.mcl_e_shaft_power = e_shaft_power[MCL]
    e_engine.mcr_e_shaft_power = e_shaft_power[MCR]
    e_engine.fid_e_shaft_power = e_shaft_power[FID]

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
def eval_hybrid_nacelle_design(aircraft):
    """
    Hybrid propulsive architecture design
    """

    design_driver = aircraft.design_driver
    fuselage = aircraft.fuselage
    wing = aircraft.wing

    propulsion = aircraft.propulsion

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    e_engine = aircraft.rear_electric_engine
    e_nacelle = aircraft.rear_electric_nacelle

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

    # Electric nacelle is design by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    dISA = 0.
    Altp = design_driver.ref_cruise_altp
    Mach = design_driver.cruise_mach

    (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(Altp,dISA)

    shaft_power = e_engine.mcr_e_shaft_power
    hub_width = 0.5     # Diameter of the e fan hub

    body_length = fuselage.length
    body_width = fuselage.width

    eval_bli_nacelle_design(e_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width,body_length,body_width)

    e_nacelle.x_axe = fuselage.length + 0.2*e_nacelle.width
    e_nacelle.y_axe = 0.
    e_nacelle.z_axe = 0.91*fuselage.height - 0.55*fuselage.height

    # Engine performance update
    #-----------------------------------------------------------------------------------------------------------
    fd = e_engine.flight_data

    e_fan_thrust = {"MTO":0., "MCN":0., "MCL":0., "MCR":0., "FID":0.}

    for rating in propulsion.rating_code:

        altp = fd.get("altp")[rating]
        disa = fd.get("disa")[rating]
        mach = fd.get("mach")[rating]
        nei = fd.get("nei")[rating]

        (Pamb,Tamb,Tstd,dTodZ) = earth.atmosphere(altp,disa)
        (fn,sec,data) = hybrid_thrust(aircraft,Pamb,Tamb,mach,rating,nei)
        (fn_core,fn_fan1,fn_fan2,dVbli_o_V,shaft_power2,fn0,shaft_power0) = data

        e_fan_thrust[rating] = fn_fan2

    e_engine.mto_e_fan_thrust = e_fan_thrust[MTO]
    e_engine.mcn_e_fan_thrust = e_fan_thrust[MCN]
    e_engine.mcl_e_fan_thrust = e_fan_thrust[MCL]
    e_engine.mcr_e_fan_thrust = e_fan_thrust[MCR]
    e_engine.fid_e_fan_thrust = e_fan_thrust[FID]

    Vair = Mach*earth.sound_speed(Tamb)

    (eFanFnBli,q1,dVbli) = jet.fan_thrust_with_bli(e_nacelle,Pamb,Tamb,Mach,shaft_power)

    (eFanFn,q0) = jet.fan_thrust(e_nacelle,Pamb,Tamb,Mach,shaft_power)

    propulsion.bli_e_thrust_factor = eFanFnBli / eFanFn     # Thrust increase due to BLI at iso shaft power for the e-fan

    propulsion.bli_thrust_factor = 1.     # Thrust increase due to BLI at iso shaft power for the turbofans (provision)

    return


#===========================================================================================================
def resize_boundary_layer(body_width,hub_width):
    """
    Compute the relation between d0 and d1
    d0 : boundary layer thickness around a tube of constant diameter
    d1 : boundary layer thickness around a the tapered part of the tube
    """

    r0 = 0.5 * body_width   # Radius of the fuselage, supposed constant
    r1 = 0.5 * hub_width    # Radius of the hub of the efan nacelle

    #===========================================================================================================
    def fct_specific_flows(d1,r1,d0,r0):
        (q0s0,q1s0,q2s0,v1s0,dvs0) = jet.specific_air_flows(r0,d0,d0)
        (q0s1,q1s1,q2s1,v1s1,dvs1) = jet.specific_air_flows(r1,d1,d1)
        y = q2s0 - q2s1
        return y
    #-----------------------------------------------------------------------------------------------------------

    n = 25
    yVein = numpy.linspace(0.001,1.50,n)

    body_bnd_layer = numpy.zeros((n,2))

    for j in range (0, n-1):
        fct1s = (r1,yVein[j],r0)
        # computation of d1 theoretical thickness of the boundary layer that passes the same air flow around the hub
        body_bnd_layer[j,0] = yVein[j]
        body_bnd_layer[j,1] = fsolve(fct_specific_flows,yVein[j],fct1s)

    return body_bnd_layer


#===========================================================================================================
def eval_bli_nacelle_design(this_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width,body_length,body_width):
    """
    BLI nacelle design
    """

    gam = earth.heat_ratio()
    r = earth.gaz_constant()
    Cp = earth.heat_constant(gam,r)

    (rho,sig) = earth.air_density(Pamb,Tamb)
    Vsnd = earth.sound_speed(Tamb)
    Re = earth.reynolds_number(Pamb,Tamb,Mach)
    Vair = Vsnd*Mach

    # Precalculation of the relation between d0 and d1
    #-----------------------------------------------------------------------------------------------------------

    body_bnd_layer = resize_boundary_layer(body_width,hub_width)

    # Electrical nacelle geometry : e-nacelle diameter is size by cruise conditions
    #-----------------------------------------------------------------------------------------------------------
    r0 = 0.5*body_width     # Radius of the fuselage, supposed constant
    d0 = jet.boundary_layer(Re,body_length)     # theoritical thickness of the boundary layer without taking account of fuselage tapering
    r1 = 0.5*hub_width        # Radius of the hub of the efan nacelle
    d1 = lin_interp_1d(d0,body_bnd_layer[:,0],body_bnd_layer[:,1])       # Thickness of the BL around the hub

    deltaV = 2.*Vair*(this_nacelle.efficiency_fan/this_nacelle.efficiency_prop - 1.)      # speed variation produced by the fan

    PwInput = this_nacelle.efficiency_fan*shaft_power     # kinetic energy produced by the fan

    #===========================================================================================================
    def fct_power_1(y,PwInput,deltaV,rho,Vair,r1,d1):
        (q0,q1,q2,v1,dVbli) = jet.air_flows(rho,Vair,r1,d1,y)
        Vinlet = Vair - dVbli
        Vjet = Vinlet + deltaV
        Pw = 0.5*q1*(Vjet**2 - Vinlet**2)
        y = PwInput - Pw
        return y
    #-----------------------------------------------------------------------------------------------------------

    fct_arg = (PwInput,deltaV,rho,Vair,r1,d1)

    # Computation of y1 : thickness of the vein swallowed by the inlet
    output_dict = fsolve(fct_power_1, x0=d1, args=fct_arg, full_output=True)

    y1 = output_dict[0][0]

    (q0,q1,q2,v1,dVbli) = jet.air_flows(rho,Vair,r1,d1,y1)

    MachInlet = v1/Vsnd     # Mean Mach number at inlet position

    Ptot = earth.total_pressure(Pamb,MachInlet)        # Stagnation pressure at inlet position

    Ttot = earth.total_temperature(Tamb,MachInlet)     # Stagnation temperature at inlet position

    MachFan = 0.5       # required Mach number at fan position

    CQoA1 = jet.corrected_air_flow(Ptot,Ttot,MachFan)        # Corrected air flow per area at fan position

    eFanArea = q1/CQoA1     # Fan area around the hub

    fan_width = numpy.sqrt(hub_width**2 + 4*eFanArea/numpy.pi)        # Fan diameter

    Vjet = v1 + deltaV      # Jet velocity

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

    this_nacelle.bnd_layer = body_bnd_layer

    this_nacelle.body_length = body_length

    return


#===========================================================================================================
def eval_hybrid_nacelle_mass(aircraft):
    """
    Hybridized propulsive nacelle mass estimations
    """

    fuselage = aircraft.fuselage

    engine = aircraft.turbofan_engine
    nacelle = aircraft.turbofan_nacelle

    e_engine = aircraft.rear_electric_engine
    e_nacelle = aircraft.rear_electric_nacelle

    power_elec = aircraft.pte1_power_elec_chain

    # Propulsion system mass is sized according max power
    # -----------------------------------------------------------------------
    e_shaft_power = numpy.array([e_engine.mto_e_shaft_power,
                                 e_engine.mcn_e_shaft_power,
                                 e_engine.mcl_e_shaft_power,
                                 e_engine.mcr_e_shaft_power,
                                 e_engine.fid_e_shaft_power])

    shaftPowerMax = max(e_shaft_power)

    turboFanMass0 = 1250. + 0.021*engine.reference_thrust # Statistical regression

    turboFanMass1 = 1250. + 0.021*engine.reference_thrust*engine.kfn_off_take

    kTurboFanMass = turboFanMass1 / turboFanMass0

    kMass = kTurboFanMass + engine.core_weight_ratio*(1-kTurboFanMass)     # Assuming core mass remains unchanged

    nacelle.mass = engine.n_engine * turboFanMass0 * kMass     # Total engine mass

    power_elec.mass = (  1./power_elec.generator_pw_density + 1./power_elec.rectifier_pw_density \
                       + 1./power_elec.wiring_pw_density + 1./power_elec.cooling_pw_density \
                       ) * shaftPowerMax

    e_nacelle.mass = (  1./e_nacelle.controller_pw_density + 1./e_nacelle.motor_pw_density \
                      + 1./e_nacelle.nacelle_pw_density \
                      ) * shaftPowerMax

    # Propulsion system CG
    # ------------------------------------------------------------------------
    nacelle.c_g = nacelle.x_ext + 0.70*nacelle.length

    power_elec.c_g = 0.70*nacelle.c_g + 0.30*fuselage.length

    e_nacelle.c_g = fuselage.length + 0.5*e_nacelle.length

    return


#===========================================================================================================
def eval_fuselage_battery_cg(aircraft):
    """
    Body battery predesign
    """

    fuselage = aircraft.fuselage

    battery = aircraft.battery

    battery.c_g = fuselage.c_g

    return


