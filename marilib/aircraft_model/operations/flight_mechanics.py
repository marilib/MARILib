#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy
from scipy.optimize import fsolve
from marilib.tools.math import maximize_1d

from marilib.earth import environment as earth
from marilib.aircraft_model.airplane import aerodynamics as aero
from marilib.airplane.propulsion import propulsion_models as propu


#===========================================================================================================
def lift_from_speed(aircraft,pamb,mach,mass):

    wing = aircraft.wing
    g = earth.gravity()
    gam = earth.heat_ratio()
    c_z = (2.*mass*g)/(gam*pamb*mach**2*wing.area)
    return c_z


#===========================================================================================================
def speed_from_lift(aircraft,pamb,cz,mass):

    wing = aircraft.wing
    g = earth.gravity()
    gam = earth.heat_ratio()
    mach = numpy.sqrt((mass*g)/(0.5*gam*pamb*wing.area*cz))
    return mach


#===========================================================================================================
def get_speed(pamb,speed_mode,mach):
    """
    retrieves CAS or mach from mach depending on speed_mode
    """

    speed = {
            1 : earth.vcas_from_mach(pamb,mach),   # CAS required
            2 : mach                               # mach required
            }.get(speed_mode, "Erreur: select speed_mode equal to 1 or 2")
    return speed


#===========================================================================================================
def get_mach(pamb,speed_mode,speed):
    """
    Retrieves mach from CAS or mach depending on speed_mode
    """

    mach = {
            1 : earth.mach_from_vcas(pamb,speed),   # Input is CAS
            2 : speed                               # Input is mach
            }.get(speed_mode, "Erreur: select speed_mode equal to 1 or 2")
    return mach


#===========================================================================================================
def acceleration(aircraft,nei,altp,disa,speed_mode,speed,mass,rating):
    """
    Aircraft acceleration on level flight
    """

    wing = aircraft.wing
    gam = earth.heat_ratio()

    [pamb,tamb,tstd,dtodz] = earth.atmosphere(altp,disa)

    mach = get_mach(pamb,speed_mode,speed)

    throttle = 1.

    fn,sfc,sec,data = propu.thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    cz = lift_from_speed(aircraft,pamb,mach,mass)

    cx,lod = aero.drag(aircraft,pamb,tamb,mach,cz)

    if(nei>0):
        dcx = propu.oei_drag(aircraft,pamb,mach)
        cx = cx + dcx*nei

    acc = (fn - 0.5*gam*pamb*mach**2*wing.area*cx) / mass

    return acc


#===========================================================================================================
def air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating):
    """
    Retrieves air path in various conditions
    """

    g = earth.gravity()

    [pamb,tamb,tstd,dtodz] = earth.atmosphere(altp,disa)

    mach = get_mach(pamb,speed_mode,speed)

    throttle = 1.

    fn,sfc,sec,data = propu.thrust(aircraft,pamb,tamb,mach,rating,throttle,nei)

    cz = lift_from_speed(aircraft,pamb,mach,mass)

    [cx,lod] = aero.drag(aircraft,pamb,tamb,mach,cz)

    if(nei>0):
        dcx = propu.oei_drag(aircraft,pamb,mach)
        cx = cx + dcx*nei
        lod = cz/cx

    acc_factor = earth.climb_mode(speed_mode,dtodz,tstd,disa,mach)

    slope = ( fn/(mass*g) - 1/lod ) / acc_factor

    vsnd = earth.sound_speed(tamb)

    v_z = mach*vsnd*slope

    return slope,v_z


#===========================================================================================================
def max_path(aircraft,nei,altp,disa,speed_mode,mass,rating):
    """
    Optimizes the speed of the aircraft to maximize the air path
    """

    #=======================================================================================
    def fct_max_path(cz,aircraft,nei,altp,disa,speed_mode,mass,rating,isformax):
        pamb,tamb,tstd,dtodz = earth.atmosphere(altp,disa)
        mach = speed_from_lift(aircraft,pamb,cz,mass)
        speed = get_speed(pamb,speed_mode,mach)
        [slope,vz] = air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating)
        if(isformax==True):
            return slope
        elif(isformax==False):
            return slope,vz,speed
    #---------------------------------------------------------------------------------------

    cz_ini = 0.5
    dcz = 0.05
    isformax=True

    fct = [fct_max_path, aircraft,nei,altp,disa,speed_mode,mass,rating,isformax]

    (cz,slope,rc) = maximize_1d(cz_ini,dcz,fct)

    isformax=False

    [slope,vz,speed] = fct_max_path(cz,aircraft,nei,altp,disa,speed_mode,mass,rating,isformax)

    return slope,vz,speed,cz


#===========================================================================================================
def propulsion_ceiling(aircraft,altp_ini,nei,vzreq,disa,speed_mode,speed,mass,rating):
    """
    Optimizes the speed of the aircraft to maximize the air path
    """

    #=======================================================================================
    def fct_prop_ceiling(altp,aircraft,nei,vzreq,disa,speed_mode,speed,mass,rating):
        [slope,vz] = air_path(aircraft,nei,altp,disa,speed_mode,speed,mass,rating)
        delta_vz = vz - vzreq
        return delta_vz
    #---------------------------------------------------------------------------------------

    fct_arg = (aircraft,nei,vzreq,disa,speed_mode,speed,mass,rating)

    altp_and_infodict = fsolve(fct_prop_ceiling, x0 = altp_ini, args=fct_arg, full_output = True) #fsolve(altp_ini,fct) ;

    altp = altp_and_infodict[0][0]
    rei = altp_and_infodict[2]
    if(rei!=1):
        altp = numpy.NaN

    return altp, rei



