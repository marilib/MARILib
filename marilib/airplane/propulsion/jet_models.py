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
