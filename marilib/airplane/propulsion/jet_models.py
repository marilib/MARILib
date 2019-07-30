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


#===========================================================================================================
def efan_nacelle_design(this_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width):
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

    CQoA1 = corrected_air_flow(Ptot,Ttot,MachFan)        # Corrected air flow per area at fan position

    eFanArea = q1/CQoA1     # Fan area around the hub

    fan_width = numpy.sqrt(hub_width**2 + 4*eFanArea/numpy.pi)        # Fan diameter

    TtotJet = Ttot + shaft_power/(q1*Cp)        # Stagnation pressure increases due to introduced work

    Tstat = TtotJet - 0.5*Vjet**2/Cp        # static temperature

    VsndJet = numpy.sqrt(gam*r*Tstat) # Sound velocity at nozzle exhaust

    MachJet = Vjet/VsndJet # Mach number at nozzle output

    PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb)

    CQoA2 = corrected_air_flow(PtotJet,TtotJet,MachJet)     # Corrected air flow per area at nozzle output

    nozzle_area = q1/CQoA2        # Fan area around the hub

    nozzle_width = numpy.sqrt(4*nozzle_area/numpy.pi)       # Nozzle diameter

    this_nacelle.hub_width = hub_width

    this_nacelle.fan_width = fan_width

    this_nacelle.nozzle_width = nozzle_width

    this_nacelle.nozzle_area = nozzle_area

    this_nacelle.width = 1.20*fan_width      # Surrounding structure

    this_nacelle.length = 1.50*this_nacelle.width

    # Total wetted area of main nacelles
    this_nacelle.net_wetted_area = numpy.pi*this_nacelle.width*this_nacelle.length*this_nacelle.n_engine

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
        (q0s0,q1s0,q2s0,v1s0,dvs0) = specific_air_flows(r0,d0,d0)
        (q0s1,q1s1,q2s1,v1s1,dvs1) = specific_air_flows(r1,d1,d1)
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
def rear_nacelle_design(this_nacelle,Pamb,Tamb,Mach,shaft_power,hub_width,body_length,body_width):
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
    d0 = boundary_layer(Re,body_length)     # theoritical thickness of the boundary layer without taking account of fuselage tapering
    r1 = 0.5*hub_width        # Radius of the hub of the efan nacelle
    d1 = lin_interp_1d(d0,body_bnd_layer[:,0],body_bnd_layer[:,1])       # Thickness of the BL around the hub

    deltaV = 2.*Vair*(this_nacelle.efficiency_fan/this_nacelle.efficiency_prop - 1.)      # speed variation produced by the fan

    PwInput = this_nacelle.efficiency_fan*shaft_power     # kinetic energy produced by the fan

    #===========================================================================================================
    def fct_power_1(y,PwInput,deltaV,rho,Vair,r1,d1):
        (q0,q1,q2,v1,dVbli) = air_flows(rho,Vair,r1,d1,y)
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

    (q0,q1,q2,v1,dVbli) = air_flows(rho,Vair,r1,d1,y1)

    MachInlet = v1/Vsnd     # Mean Mach number at inlet position

    Ptot = earth.total_pressure(Pamb,MachInlet)        # Stagnation pressure at inlet position

    Ttot = earth.total_temperature(Tamb,MachInlet)     # Stagnation temperature at inlet position

    MachFan = 0.5       # required Mach number at fan position

    CQoA1 = corrected_air_flow(Ptot,Ttot,MachFan)        # Corrected air flow per area at fan position

    eFanArea = q1/CQoA1     # Fan area around the hub

    fan_width = numpy.sqrt(hub_width**2 + 4*eFanArea/numpy.pi)        # Fan diameter

    Vjet = v1 + deltaV      # Jet velocity

    TtotJet = Ttot + shaft_power/(q1*Cp)        # Stagnation pressure increases due to introduced work

    Tstat = TtotJet - 0.5*Vjet**2/Cp        # static temperature

    VsndJet = numpy.sqrt(gam*r*Tstat) # Sound velocity at nozzle exhaust

    MachJet = Vjet/VsndJet # Mach number at nozzle output

    PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb)

    CQoA2 = corrected_air_flow(PtotJet,TtotJet,MachJet)     # Corrected air flow per area at nozzle output

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
def fan_thrust_with_bli(nacelle,Pamb,Tamb,Mach,PwShaft):
    """
    Compute the thrust of a fan of a given geometry swallowing
    the boundary layer (BL) of a body of a given geometry
    The amount of swallowed BL depends on the given shaft power and flying
    conditions.
    """

    gam = earth.heat_ratio()
    r = earth.gaz_constant()
    Cp = earth.heat_constant(gam,r)

    #===========================================================================================================
    def fct_power_bli(y,PwShaft,Pamb,rho,Ttot,Vair,r1,d1,nozzle_area):

        (q0,q1,q2,Vinlet,dVbli) = air_flows(rho,Vair,r1,d1,y)
        PwInput = nacelle.efficiency_fan * PwShaft
        Vjet = numpy.sqrt(2.*PwInput/q1 + Vinlet**2)
        TtotJet = Ttot + PwShaft/(q1*Cp)        # Stagnation temperature increases due to introduced work
        Tstat = TtotJet - 0.5*Vjet**2/Cp        # Static temperature
        VsndJet = earth.sound_speed(Tstat)     # Sound speed at nozzle exhaust
        MachJet = Vjet/VsndJet                  # Mach number at nozzle output
        PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb) supposing adapted nozzle
        CQoA1 = corrected_air_flow(PtotJet,TtotJet,MachJet)    # Corrected air flow per area at nozzle position
        q = CQoA1*nozzle_area

        y = q1 - q

        return y
    #-----------------------------------------------------------------------------------------------------------

    nozzle_area = nacelle.nozzle_area
    bnd_layer = nacelle.bnd_layer

    Re = earth.reynolds_number(Pamb,Tamb,Mach)

    d0 = boundary_layer(Re,nacelle.body_length)      # theorical thickness of the boundary layer without taking account of fuselage tapering
    r1 = 0.5*nacelle.hub_width      # Radius of the hub of the eFan nacelle
    d1 = lin_interp_1d(d0,bnd_layer[:,0],bnd_layer[:,1])     # Using the precomputed relation

    Ttot = earth.total_temperature(Tamb,Mach)      # Stagnation temperature at inlet position
    rho,sig = earth.air_density(Pamb,Tamb)
    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    fct_arg = (PwShaft,Pamb,rho,Ttot,Vair,r1,d1,nozzle_area)

    # Computation of y1 : thikness of the vein swallowed by the inlet
    output_dict = fsolve(fct_power_bli, x0=0.50, args=fct_arg, full_output=True)

    y = output_dict[0][0]
    if (output_dict[2]!=1):
        raise Exception("Convergence problem")

    (q0,q1,q2,Vinlet,dVbli) = air_flows(rho,Vair,r1,d1,y)
    PwInput = nacelle.efficiency_fan * PwShaft
    Vjet = numpy.sqrt(2.*PwInput/q1 + Vinlet**2)

    eFn = q1*(Vjet - Vinlet)

    return (eFn,q1,dVbli)


#===========================================================================================================
def fan_thrust(nacelle,Pamb,Tamb,Mach,PwShaft):
    """
    Compute the thrust of a fan of given geometry swallowing free air stream
    """

    gam = earth.heat_ratio()
    r = earth.gaz_constant()
    Cp = earth.heat_constant(gam,r)

    #===========================================================================================================
    def fct_power(q,PwShaft,Pamb,Ttot,Vair,NozzleArea):

        Vinlet = Vair
        PwInput = nacelle.efficiency_fan*PwShaft
        Vjet = numpy.sqrt(2.*PwInput/q + Vinlet**2)         # Supposing isentropic compression
        TtotJet = Ttot + PwShaft/(q*Cp)        # Stagnation temperature increases due to introduced work
        TstatJet = TtotJet - 0.5*Vjet**2/Cp        # Static temperature
        VsndJet = earth.sound_speed(TstatJet)     # Sound speed at nozzle exhaust
        MachJet = Vjet/VsndJet                  # Mach number at nozzle output
        PtotJet = earth.total_pressure(Pamb,MachJet)       # total pressure at nozzle exhaust (P = Pamb)
        CQoA1 = corrected_air_flow(PtotJet,TtotJet,MachJet)       # Corrected air flow per area at fan position
        q0 = CQoA1*NozzleArea

        y = q0 - q

        return y
    #-----------------------------------------------------------------------------------------------------------

    NozzleArea = nacelle.nozzle_area
    FanWidth = nacelle.fan_width

    Ptot = earth.total_pressure(Pamb,Mach)        # Total pressure at inlet position
    Ttot = earth.total_temperature(Tamb,Mach)     # Total temperature at inlet position

    Vsnd = earth.sound_speed(Tamb)
    Vair = Vsnd*Mach

    fct_arg = (PwShaft,Pamb,Ttot,Vair,NozzleArea)

    CQoA0 = corrected_air_flow(Ptot,Ttot,Mach)       # Corrected air flow per area at fan position
    q0init = CQoA0*(0.25*numpy.pi*FanWidth**2)

    # Computation of the air flow swallowed by the inlet
    output_dict = fsolve(fct_power, x0=q0init, args=fct_arg, full_output=True)

    q0 = output_dict[0][0]
    if (output_dict[2]!=1):
        raise Exception("Convergence problem")

    Vinlet = Vair
    PwInput = nacelle.efficiency_fan*PwShaft
    Vjet = numpy.sqrt(2.*PwInput/q0 + Vinlet**2)

    eFn = q0*(Vjet - Vinlet)

    return (eFn,q0)


#===========================================================================================================
def corrected_air_flow(Ptot,Ttot,Mach):
    """
    Computes the corrected air flow per square meter
    """

    R = earth.gaz_constant()
    gam = earth.heat_ratio()

    f_M = Mach*(1. + 0.5*(gam-1)*Mach**2)**(-(gam+1.)/(2.*(gam-1.)))

    CQoA = (numpy.sqrt(gam/R)*Ptot/numpy.sqrt(Ttot))*f_M

    return CQoA


#===========================================================================================================
def air_flows(rho,v_air,r,d,y):
    """
    Air flows and speeds at rear end of a cylinder of radius rear_radius mouving at v_air in the direction of its axes
    y is the elevation upon the surface of the cylinder : 0 < y < inf
    """

    # exponent in the formula of the speed profile inside a turbulent BL of thickness bly : Vy/Vair = (y/d)**(1/7)
    n = 1./7.

    # Cumulated air flow at y_elev, without BL
    q0 = (2.*numpy.pi)*(rho*v_air)*(r*y + 0.5*y**2)

    ym = min(y,d)

    # Cumulated air flow at y_elev, with BL
    q1 = (2.*numpy.pi)*(rho*v_air)*d*( (r/(n+1))*(ym/d)**(n+1) + (d/(n+2))*(ym/d)**(n+2) )

    if (y>d):
        # Add to Q1 the air flow outside the BL
        q1 = q1 + q0 - (2.*numpy.pi)*(rho*v_air)*( r*d + 0.5*d**2 )

    q2 = q1 - q0        # Cumulated air flow at y_elev, inside the BL (going speed wise)

    v1 = v_air*(q1/q0)     # Mean speed of q1 air flow at y_elev

    dv = v_air - v1       # Mean air flow speed variation at y_elev

    return q0,q1,q2,v1,dv


#===========================================================================================================
def specific_air_flows(r,d,y):
    """
    Specific air flows and speeds at rear end of a cylinder of radius R
    mouving at Vair in the direction of Qs = Q/(rho*Vair)     Vs = V/Vair
    its axes, y is the elevation upon the surface of the cylinder :
                              0 < y < inf
    WARNING : even if all mass flows are positive,
    Q0 and Q1 are going backward in fuselage frame, Q2 is going forward
    in ground frame
    """

    n = 1/7 # exponent in the formula of the speed profile inside a turbulent
            # BL of thickness d : Vy/Vair = (y/d)^(1/7)

    q0s = (2.*numpy.pi)*( r*y + 0.5*y**2 )
                            # Cumulated specific air flow at y, without BL
    ym = min(y,d)

    q1s = (2.*numpy.pi)*d*( (r/(n+1))*(ym/d)**(n+1) + (d/(n+2))*(ym/d)**(n+2) )
                            # Cumulated specific air flow at y, without BL
    if y>d:

        q1s = q1s + q0s - (2.*numpy.pi)*( r*d + 0.5*d**2 )
                            # Add to Q1 the specific air flow outside the BL
    q2s = q0s - q1s
        # Cumulated specific air flow at y, inside the BL (going speed wise)
    v1s = (q1s/q0s) # Mean specific speed of Q1 air flow at y

    dVs = (1 - v1s) # Mean specific air flow spped variation at y

    return (q0s,q1s,q2s,v1s,dVs)


#===========================================================================================================
def boundary_layer(re,x_length):
    """
    Thickness of a turbulent boundary layer which developped turbulently from its starting point
    """

    d = (0.385*x_length)/(re*x_length)**(1./5.)

    return d
