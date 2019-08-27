#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


from marilib import numpy
from marilib.tools.math import maximize_1d, trinome, vander3
from marilib.tools import units as unit

from marilib.earth import environment as earth
from marilib.aircraft_model.airplane import aerodynamics as airplane_aero, regulation as regul
from marilib.airplane.propulsion import propulsion_models as propu
from marilib.aircraft_model.operations import flight_mechanics as flight


#===========================================================================================================
def ceilings(aircraft,toc,oei_ceil):

    design_driver = aircraft.design_driver
    propulsion = aircraft.propulsion
    weights = aircraft.weights

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    disa = 15.

    # propulsion ceilings
    #-----------------------------------------------------------------------------------------------------------
    nei = 0
    altp = toc
    speed_mode = 2                      # WARNING : iso Mach climb mode
    speed = design_driver.cruise_mach

    mass = 0.97*weights.mtow

    rating = MCL    # Max Climb

    slope, vz_clb = flight.air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating)

    rating = MCR    # Max Cruise

    slope, vz_crz = flight.air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating)

    # One engine inoperative ceiling
    #-----------------------------------------------------------------------------------------------------------
    nei = 1
    altp = oei_ceil
    speed_mode = 2                      # WARNING : iso Mach climb mode

    pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)

    mass = 0.95*weights.mtow

    rating = MCN

    oei_slope,vz,oei_mach,cz = flight.max_path(aircraft,nei,altp,disa,speed_mode,mass,rating)

    oei_speed = earth.vcas_from_mach(pamb,oei_mach)

    #-----------------------------------------------------------------------------------------------------------
    return vz_clb,vz_crz,oei_slope,oei_speed


#===========================================================================================================
def time_to_climb(aircraft,toc,disa,mass,vcas1,vcas2,mach):
    """
    Time to climb to initial cruise altitude
    For simplicity reasons, airplane mass is supposed constant
    """

    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    if(vcas1>unit.mps_kt(250.)):
        print("vcas1 = ",unit.kt_mps(vcas1))
        print("time_to_climb_, vcas1 must be lower than or equal to 250kt")
    if(vcas1>vcas2):
        print("vcas1 = ",unit.kt_mps(vcas1))
        print("vcas2 = ",unit.kt_mps(vcas2))
        print("time_to_climb_, vcas1 must be lower than or equal to vcas2")

    cross_over_altp = earth.cross_over_altp(vcas2,mach)

    if(cross_over_altp<unit.m_ft(1500.)):
        print("time_to_climb_, cross over altitude is too low")

    if(toc<cross_over_altp):
        cross_over_altp = toc

    # Duration of initial climb
    #-----------------------------------------------------------------------------------------------------------
    altp0 = unit.m_ft(1500.)
    altp2 = unit.m_ft(10000.)
    altp1 = (altp0+altp2)/2.
    altp = numpy.array([altp0, altp1, altp2])

    nei = 0
    speed_mode = 1    # Constant CAS
    rating = MCL

    [slope,v_z0] = flight.air_path(aircraft,nei,altp[0],disa,speed_mode,vcas1,mass,rating)
    [slope,v_z1] = flight.air_path(aircraft,nei,altp[1],disa,speed_mode,vcas1,mass,rating)
    [slope,v_z2] = flight.air_path(aircraft,nei,altp[2],disa,speed_mode,vcas1,mass,rating)
    v_z = numpy.array([v_z0, v_z1, v_z2])

    if (v_z[0]<0. or v_z[1]<0. or v_z[2]<0.):
        print("time_to_climb_, Climb to acceleration altitude is not possible")

    A = vander3(altp)
    B = 1./v_z
    C = trinome(A,B)

    time1 = ((C[0]*altp[2]/3. + C[1]/2.)*altp[2] + C[2])*altp[2]
    time1 = time1 - ((C[0]*altp[0]/3. + C[1]/2.)*altp[0] + C[2])*altp[0]

    # Acceleration
    #-----------------------------------------------------------------------------------------------------------
    vc0 = vcas1
    vc2 = vcas2
    vc1 = (vc0+vc2)/2.
    vcas = numpy.array([vc0, vc1, vc2])

    acc0 = flight.acceleration(aircraft,nei,altp[2],disa,speed_mode,vcas[0],mass,rating)
    acc1 = flight.acceleration(aircraft,nei,altp[2],disa,speed_mode,vcas[1],mass,rating)
    acc2 = flight.acceleration(aircraft,nei,altp[2],disa,speed_mode,vcas[2],mass,rating)
    acc = numpy.array([acc0, acc1, acc2])

    if(acc[0]<0. or acc[1]<0. or acc[2]<0.):
        print("time_to_climb_, acceleration is not possible")

    A = vander3(vcas)
    B = 1./acc
    C = trinome(A,B)

    time2 = ((C[0]*vcas[2]/3. + C[1]/2.)*vcas[2] + C[2])*vcas[2]
    time2 = time2 - ((C[0]*vcas[0]/3. + C[1]/2.)*vcas[0] + C[2])*vcas[0]

    # Duration of climb to cross over
    #-----------------------------------------------------------------------------------------------------------
    altp0 = unit.m_ft(10000.)
    altp2 = cross_over_altp
    altp1 = (altp0+altp2)/2.
    altp = numpy.array([altp0, altp1, altp2])

    [slope,v_z0] = flight.air_path(aircraft,nei,altp[0],disa,speed_mode,vcas2,mass,rating)
    [slope,v_z1] = flight.air_path(aircraft,nei,altp[1],disa,speed_mode,vcas2,mass,rating)
    [slope,v_z2] = flight.air_path(aircraft,nei,altp[2],disa,speed_mode,vcas2,mass,rating)
    v_z = numpy.array([v_z0, v_z1, v_z2])

    if(v_z[0]<0. or v_z[1]<0. or v_z[2]<0.):
        print("time_to_climb_, Climb to cross over altitude is not possible")

    A = vander3(altp)
    B = 1./v_z
    C = trinome(A,B)

    time3 = ((C[0]*altp[2]/3. + C[1]/2.)*altp[2] + C[2])*altp[2]
    time3 = time3 - ((C[0]*altp[0]/3. + C[1]/2.)*altp[0] + C[2])*altp[0]

    # Duration of climb to altp
    #-----------------------------------------------------------------------------------------------------------
    if(cross_over_altp<toc):
        altp0 = cross_over_altp
        altp2 = toc
        altp1 = (altp0+altp2)/2.
        altp = numpy.array([altp0, altp1, altp2])

        speed_mode = 2    # mach

        [slope,v_z0] = flight.air_path(aircraft,nei,altp[0],disa,speed_mode,mach,mass,rating)
        [slope,v_z1] = flight.air_path(aircraft,nei,altp[1],disa,speed_mode,mach,mass,rating)
        [slope,v_z2] = flight.air_path(aircraft,nei,altp[2],disa,speed_mode,mach,mass,rating)
        v_z = numpy.array([v_z0, v_z1, v_z2])

        if(v_z[0]<0. or v_z[1]<0. or v_z[2]<0.):
            print("time_to_climb_, Climb to top of climb is not possible")

        A = vander3(altp)
        B = 1./v_z
        C = trinome(A,B)

        time4 =  ((C[0]*altp[2]/3. + C[1]/2.)*altp[2] + C[2])*altp[2] \
               - ((C[0]*altp[0]/3. + C[1]/2.)*altp[0] + C[2])*altp[0]
    else:
        time4 = 0.

    #    Total time
    #-----------------------------------------------------------------------------------------------------------
    ttc = time1 + time2 + time3 + time4

    return ttc


#===========================================================================================================
def take_off(aircraft,kvs1g,altp,disa,mass,hld_conf):
    """
    Take off field length and climb path at 35 ft depending on stall margin (kVs1g)
    """

    wing = aircraft.wing
    propulsion = aircraft.propulsion

    (MTO,MCN,MCL,MCR,FID) = propulsion.rating_code

    czmax,cz_0 = airplane_aero.high_lift(wing,hld_conf)

    rating = MTO

    [pamb,tamb,tstd,dtodz] = earth.atmosphere(altp,disa)

    [rho,sig] = earth.air_density(pamb,tamb)

    cz_to = czmax / kvs1g**2

    mach = flight.speed_from_lift(aircraft,pamb,tamb,cz_to,mass)

    throttle = 1.

    nei = 0    # For Magic Line factor computation

    fn,sfc,sec,data = propu.thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    ml_factor = mass**2 / (cz_to*fn*wing.area*sig**0.8 )  # Magic Line factor

    tofl = 15.5*ml_factor + 100.

    nei = 1    # For 2nd segment computation
    speed_mode = 1
    speed = flight.get_speed(pamb,speed_mode,mach)

    seg2path,vz = flight.air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating)

    return seg2path,tofl


#===========================================================================================================
def take_off_field_length(aircraft,altp,disa,mass,hld_conf):
    """
    Take off field length and climb path with eventual kVs1g increase to recover min regulatory slope
    """

    kvs1g = regul.kvs1g_min_take_off()

    [seg2_path,tofl] = take_off(aircraft,kvs1g,altp,disa,mass,hld_conf)

    n_engine = aircraft.propulsion.n_engine

    seg2_min_path = regul.seg2_min_path(n_engine)

    if(seg2_min_path<seg2_path):
        limitation = 1
    else:
        dkvs1g = 0.005
        kvs1g_ = numpy.array([0.,0.])
        kvs1g_[0] = kvs1g
        kvs1g_[1] = kvs1g_[0] + dkvs1g

        seg2_path_ = numpy.array([0.,0.])
        seg2_path_[0] = seg2_path
        seg2_path_[1],trash = take_off(aircraft,kvs1g_[1],altp,disa,mass,hld_conf)

        while(seg2_path_[0]<seg2_path_[1] and seg2_path_[1]<seg2_min_path):
            kvs1g_[0] = kvs1g_[1]
            kvs1g_[1] = kvs1g_[1] + dkvs1g
            seg2_path_[1],trash = take_off(aircraft,kvs1g_[1],altp,disa,mass,hld_conf)

        if(seg2_min_path<seg2_path_[1]):
            kvs1g = kvs1g_[0] + ((kvs1g_[1]-kvs1g_[0])/(seg2_path_[1]-seg2_path_[0]))*(seg2_min_path-seg2_path_[0])
            [seg2_path,tofl] = take_off(aircraft,kvs1g,altp,disa,mass,hld_conf)
            seg2_path = seg2_min_path
            limitation = 2
        else:
            tofl = numpy.nan
            kvs1g = numpy.nan
            seg2_path = 0.
            limitation = 0

    return tofl,seg2_path,kvs1g,limitation


#===========================================================================================================
def approach_speed(aircraft,altp,disa,mass,hld_conf):
    """
    Minimum approach speed (VLS)
    """

    wing = aircraft.wing

    g = earth.gravity()

    czmax,trash = airplane_aero.high_lift(wing,hld_conf)

    stall_margin = regul.kvs1g_min_landing()

    [pamb,tamb,tstd,dtodz] = earth.atmosphere(altp,disa)

    [rho,sig] = earth.air_density(pamb,tamb)

    vapp = numpy.sqrt((mass*g) / (0.5*rho*wing.area*(czmax / stall_margin**2)))

    return vapp



