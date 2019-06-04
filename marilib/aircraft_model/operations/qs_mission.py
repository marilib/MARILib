#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""


import numpy
import copy
from scipy.optimize import fsolve,minimize
from marilib.tools.math import maximize_1d, trinome
from marilib.tools import units as unit

from marilib.earth import environment as earth
from marilib.earth import environment_grad as g_earth
from marilib.aircraft_model.airplane import aerodynamics as airplane_aero
from marilib.airplane.propulsion import propulsion_models as propu
from marilib.aircraft_model.operations import flight_mechanics as flight

state_dict = {"time":0, "mass":1, "xg":2, "zg":3, "vgnd":4, "path":5}

xout_dict = {"time":0, "mass":1, "xg":2, "zg":3, "vgnd":4, "path":5,
             "time_d":6, "mass_d":7, "xg_d":8, "zg_d":9, "vgnd_d":10, "path_d":11,
             "vxgnd_d":12, "vzgnd_d":13,
             "vair":14, "vair_d":15, "vxair":16, "vxair_d":17, "vzair":18, "vzair_d":19,
             "mach":20, "mach_d":21, "vcas":22, "vcas_d":23,
             "altp":24, "altp_d":25, "pamb":26, "pamb_d":27,
             "tamb":28, "tamb_d":29, "disa":30, "disa_d":31, "rho":32, "rho_d":33,
             "wx":34, "wx_d":35, "wz":36, "wz_d":37,
             "cz":38, "cx":39, "fn":40}


#===========================================================================================================
def cabin_virtual_altp():
    """
    Virtual pressure altitude of the cabin
    """
    altp_cabin = 2000.
    return altp_cabin


#===========================================================================================================
def cabin_pressure_grad(altp_d):

    disa = 0.
    disa_d = 0.
    altp = 0.

    pamb,pamb_d,tamb,tamb_d,dtodz = g_earth.atmosphere_grad(altp,altp_d,disa,disa_d)

    return pamb_d


#===========================================================================================================
def air(xg,zg, xg_d,zg_d):
    """
    Air field description : pressure, temperature and wind versus position and speed
    """

    disa = 0.
    disa_d = 0.

    pamb,pamb_d,tamb,tamb_d,dtodz,dtodz_d = g_earth.atmosphere_geo_grad(zg,zg_d,disa,disa_d)# Standard atmosphere

    wx = 0.
    wx_d = 0.

    wz = 0.
    wz_d = 0.

    return pamb,pamb_d, tamb,tamb_d, disa,disa_d, wx,wx_d, wz,wz_d


#===========================================================================================================
def state_dot(xin,state,rating,nei,aircraft):

    g = earth.gravity()
    gam = earth.heat_ratio()

    area = aircraft.wing.area

    cz,fn = xin
    t,mass,xg,zg,vgnd,path = state

    xg_d = vgnd*numpy.cos(path)
    zg_d = vgnd*numpy.sin(path)

    pamb,pamb_d, tamb,tamb_d, disa,disa_d, wx,wx_d, wz,wz_d = air(xg,zg, xg_d,zg_d)

    rho,rho_d, sig,sig_d = g_earth.air_density_grad(pamb,pamb_d, tamb,tamb_d)

    vsnd,vsnd_d = g_earth.sound_speed_grad(tamb,tamb_d)

    altp,altp_d = g_earth.pressure_altitude_grad(pamb,pamb_d)

# Compute state_dot
#-----------------------------------------------------------------------------------------------------------
    vxair = xg_d - wx   # Local wind is introduced here
    vzair = zg_d - wz

    vair = numpy.sqrt(vxair**2 + vzair**2)
    mach = vair/vsnd

    cx,lod = airplane_aero.drag(aircraft,pamb,tamb,mach,cz)    # Aerodynamic drag model

    path_d = ((0.5*rho*vair**2)*area*cz - mass*g*numpy.cos(path))/(mass*vair)    # Lift equation
    vgnd_d = (fn - mass*g*numpy.sin(path) - (0.5*rho*vair**2)*area*cx)/mass      # Drag equation

    sfc = propu.sfc(aircraft,pamb,tamb,mach,rating,nei)     # Propulsion consumption model

    mass_d = -sfc*fn

    state_d = numpy.array([1.,mass_d,xg_d,zg_d,vgnd_d,path_d])

# Compute other derivatives
#-----------------------------------------------------------------------------------------------------------
    vxgnd_d = vgnd_d*numpy.cos(path) - vgnd*numpy.sin(path)*path_d
    vzgnd_d = vgnd_d*numpy.sin(path) + vgnd*numpy.cos(path)*path_d

    vxair_d = vxgnd_d - wx_d
    vzair_d = vzgnd_d - wz_d

    vair_d = (vxair*vxair_d + vzair*vzair_d)/vair
    mach_d = vair_d/vsnd - vair*vsnd_d/vsnd**2

    vcas,vcas_d = g_earth.vcas_from_mach_grad(pamb,pamb_d, mach,mach_d)

    xout = numpy.array([t,mass,xg,zg,vgnd,path,
                        1.,mass_d,xg_d,zg_d,vgnd_d,path_d,
                        vxgnd_d,vzgnd_d,
                        vair,vair_d, vxair,vxair_d, vzair,vzair_d,
                        mach,mach_d, vcas,vcas_d,
                        altp,altp_d, pamb,pamb_d,
                        tamb,tamb_d, disa,disa_d, rho,rho_d,
                        wx,wx_d, wz,wz_d,
                        cz,cx,fn])

    return xout,state_d


#===========================================================================================================
def fct_thrust(aircraft,state,rating,nei):

    xg = state[state_dict["xg"]]
    zg = state[state_dict["zg"]]

    xg_d = state[state_dict["vgnd"]]*numpy.cos(state[state_dict["path"]])
    zg_d = state[state_dict["vgnd"]]*numpy.sin(state[state_dict["path"]])

    pamb,pamb_d, tamb,tamb_d, disa,disa_d, wx,wx_d, wz,wz_d = air(xg,zg, xg_d,zg_d)

    vsnd,vsnd_d = g_earth.sound_speed_grad(tamb,tamb_d)

    vair = numpy.sqrt((xg_d - wx)**2 + (zg_d - wz)**2)

    mach = vair/vsnd

    fn,prop_data = propu.thrust(aircraft,pamb,tamb,mach,rating,nei)

    return fn


#===========================================================================================================
def flight_point(aircraft,rating,nei,state,var,input,dat,cmd,dof,cst):

    #===========================================================================================================
    def fct_flight_point(var,state,input,dat,cmd,dof,cst,nei,rating,aircraft):

        for i in range(n_dat):
            param[dat[i][0]] = eval(dat[i][1])

        for i in range(n_cmd):
             param[cmd[i]] = var[i]

        xin = numpy.array([param["cz"],param["fn"]])

        for i in range(n_dof):
            state[state_dict[dof[i]]] = var[n_cmd+i]

        xout,state_d = state_dot(xin,state,rating,nei,aircraft)

        res = [None]*n_cst
        for i in range(n_cst):
            res[i] = xout[xout_dict[cst[i][0]]] - eval(cst[i][1])

        return res
    #-----------------------------------------------------------------------------------------------------------

    # Initialize level flight
    #-----------------------------------------------------------------------------------------------------------
    n_dat = numpy.shape(dat)[0]
    n_cmd = numpy.shape(cmd)[0]
    n_dof = numpy.shape(dof)[0]
    n_cst = numpy.shape(cst)[0]

    param = {"cz":0.,"fn":0.}

    fct_args = (state,input,dat,cmd,dof,cst,nei,rating,aircraft)

    out_dict = fsolve(fct_flight_point, x0=var, args=fct_args, full_output=True)

    rei = out_dict[2]
    if(rei!=1):
        raise Exception("Impossible to converge flight point",unit.NM_m(state[state_dict["xg"]]))

    for i in range(n_dat):
        param[dat[i][0]] = eval(dat[i][1])

    for i in range(n_cmd):
        param[cmd[i]] = out_dict[0][i]

    xin = numpy.array([param["cz"],param["fn"]])

    for i in range(n_dof):
        state[state_dict[dof[i]]] = out_dict[0][n_cmd+i]

    xout,state_d = state_dot(xin,state,rating,nei,aircraft)

    return xout,state,state_d


#===========================================================================================================
def flight_segment(aircraft,rating,nei,state,stop,input,dat,cmd,dof,cst,dt,flight):

    # Initialize level flight
    #-----------------------------------------------------------------------------------------------------------
    xout = flight[-1,:]

    var = init_var(aircraft,state,cmd,dof)

    stop_idx = xout_dict[stop[0]]
    stop_val = eval(stop[2])
    op = stop[1]

    # Play level flight
    #-----------------------------------------------------------------------------------------------------------
    while (eval("xout[stop_idx]"+op+"stop_val")):

        xout,state,state_d = flight_point(aircraft,rating,nei,state,var,input,dat,cmd,dof,cst)

        flight = numpy.vstack([flight,xout])

        state = state + state_d*dt

    # Adjust last point using interpolation
    #-----------------------------------------------------------------------------------------------------------
    xout = flight[-2,:] + (flight[-1,:]-flight[-2,:]) \
                         *(stop_val - flight[-2,xout_dict[stop[0]]]) \
                         /(flight[-1,xout_dict[stop[0]]] - flight[-2,xout_dict[stop[0]]])

    flight[-1,:] = xout
    state = xout[0:6]

    return flight,state


#===========================================================================================================
def iso_cas_climb(aircraft,vcas1,altp1,vcas2,altp2,mach,vz_mcl_min,flight):
    """
    Constant CAS climb :
    vcas1 from altp0 to altp1
    vcas2 from alpt1 to lower altitude from altp2 and crossover altitude with iso mach
    """

    # Initialize climb flight
    #-----------------------------------------------------------------------------------------------------------
    xout = flight[-1,:]
    state = xout[0:6]
    altp0 = xout[xout_dict["altp"]]

    altp_cross_over = earth.cross_over_altp(vcas2,mach)

    altp_stop = min(altp2,altp_cross_over)

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    nei = 0.

    # Data filtering
    #-----------------------------------------------------------------------------------------------------------
    if(vcas1>vcas2):
        raise Exception("iso_cas_climb, vcas2 must be higher or equal than vcas1")

    if(altp0>altp1):
        raise Exception("iso_cas_climb, altp1 must be higher than starting altitude")

    if(altp1>altp2):
        raise Exception("iso_cas_climb, altp2 must be higher than altp1")

    if(altp1>altp_cross_over):
        raise Exception("iso_cas_climb, cross over altitude must be higher than altp1")

    # First climb segment
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([vcas1,altp1])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["vgnd","path"])
    cst = numpy.array([["vcas","input[0]"],["vcas_d","0"],["path_d","0"]])
    stop = numpy.array(["altp","<","input[1]"])
    dt = 30.

    flight,state = flight_segment(aircraft,MCL,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    # Acceleration at constant altitude
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([altp1,vcas2])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["zg","path"])
    cst = numpy.array([["altp","input[0]"],["altp_d","0"],["path_d","0"]])
    stop = numpy.array(["vcas","<","input[1]"])
    dt = 5.

    flight,state = flight_segment(aircraft,MCL,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    # Second climb segment
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([vcas2,altp_stop])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["vgnd","path"])
    cst = numpy.array([["vcas","input[0]"],["vcas_d","0"],["path_d","0"]])
    stop = numpy.array(["altp","<","input[1]"])
    dt = 30.

    flight,state = flight_segment(aircraft,MCL,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    if(flight[-1,xout_dict["vzair"]]<vz_mcl_min):
        raise Exception("iso_cas_climb, final altitude cannot be reached with minimum required vertical speed")

    return flight,state


#===========================================================================================================
def level_flight(aircraft,mach,vz_mcl_min,vz_mcr_min,heading,stop,flight):

    westbound = numpy.concatenate((numpy.arange(2000.,41000.,2000.),numpy.arange(43000.,85000.,4000.)))
    eastbound = numpy.concatenate((numpy.arange(1000.,42000.,2000.),numpy.arange(45000.,87000.,4000.)))
    ladder = {"West":westbound, "East":eastbound}

    # Initialize level flight
    #-----------------------------------------------------------------------------------------------------------
    xout = flight[-1,:]
    state = xout[0:6]
    altp = xout[xout_dict["altp"]]
    altp_up = ladder[heading][numpy.where(ladder[heading]>(altp+1))[0][0]]      # Lower altitude imediatly higher than altp1

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    nei = 0.

    # Test cruise ceiling
    #-----------------------------------------------------------------------------------------------------------
    vz_mcr,vz_mcl = test_ceilings(aircraft,nei,mach,altp,state)

    if (vz_mcl < vz_mcl_min):
        raise Exception("Cannot start level flight because of MCL ceiling")

    if (vz_mcr < vz_mcr_min):
        raise Exception("Cannot start level flight because of MCR ceiling")

    active_ceiling = 0

    # Initialize specific ranges
    #-----------------------------------------------------------------------------------------------------------
    sar,sar_up = test_air_range(aircraft,MCR,nei,mach,altp,altp_up,state)

    dtime = 0.  # Allows to start by a climb by merging transition point to last iso_cas climb

    # Fly step cruise
    #-----------------------------------------------------------------------------------------------------------
    stop_idx = xout_dict[stop[0]]
    stop_val = eval(stop[2])
    op = stop[1]

    while (eval("flight[-1,stop_idx]"+op+"stop_val")):

        sar0 = sar
        sar_up0 = sar_up

        sar,sar_up = test_air_range(aircraft,MCR,nei,mach,altp,altp_up,state)

        if (sar < sar_up):      # Efficiency is better on upper level

            vz_mcr,vz_mcl = test_ceilings(aircraft,nei,mach,altp_up,state)

            if (vz_mcr_min<vz_mcr and vz_mcl_min<vz_mcl):     # Ceilings allow to climb

                if (0<active_ceiling):
                    # Compute transition point according to climb speed
                    #-----------------------------------------------------------------------------------------------------------
                    if (active_ceiling==1):
                        time =  flight[-1,xout_dict["time"]] - dtime/(1+numpy.abs((vz_mcr_min-active_vz)/(vz_mcr-vz_mcr_min)))
                    else:
                        time =  flight[-1,xout_dict["time"]] - dtime/(1+numpy.abs((vz_mcl_min-active_vz)/(vz_mcl-vz_mcl_min)))
                else:
                    # Compute transition point according to SAR
                    #-----------------------------------------------------------------------------------------------------------
                    time =  flight[-1,xout_dict["time"]] - dtime/(1+numpy.abs((sar0-sar_up0)/(sar_up-sar)))

                xout = flight[-2,:] + (flight[-1,:]-flight[-2,:]) \
                                     *(time - flight[-2,xout_dict["time"]]) \
                                     /(flight[-1,xout_dict["time"]] - flight[-2,xout_dict["time"]])

                flight[-1,:] = xout
                state = xout[0:6]

                # Climb to upper altitude
                #-----------------------------------------------------------------------------------------------------------
                dt = 30.
                input = numpy.array([mach,altp_up])
                dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
                cmd = numpy.array(["cz"])
                dof = numpy.array(["vgnd","path"])
                cst = numpy.array([["mach","input[0]"],["mach_d","0"],["path_d","0"]])
                stop = numpy.array(["altp","<","input[1]"])

                flight,state = flight_segment(aircraft,MCL,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

                altp = altp_up
                altp_up = ladder[heading][numpy.where(ladder[heading]>(altp_up+1))[0][0]]

                dtime = 0.  # Allows to continue climbing if best recheable altitude is not found for cruise start

            else:
                # Detect the most constraining condition
                #-----------------------------------------------------------------------------------------------------------
                if ((vz_mcl_min-vz_mcl)<(vz_mcr_min-vz_mcr)):
                    active_ceiling = 1
                    active_vz = vz_mcr
                else:
                    active_ceiling = 2
                    active_vz = vz_mcl

                # Go to next point
                #-----------------------------------------------------------------------------------------------------------
                dt = 120.
                time_stop = state[state_dict["time"]] + dt
                input = numpy.array([altp,mach,time_stop])
                dat = numpy.array([])
                cmd = numpy.array(["cz","fn"])
                dof = numpy.array(["zg","vgnd","path"])
                cst = numpy.array([["altp","input[0]"],["altp_d","0"],["mach","input[1]"],["mach_d","0"],["path_d","0"]])
                stop = numpy.array(["time","<","input[2]"])

                flight,state = flight_segment(aircraft,MCR,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

                dtime = dt

        else:

            # Go to next point
            #-----------------------------------------------------------------------------------------------------------
            dt = 120.
            time_stop = state[state_dict["time"]] + dt
            input = numpy.array([altp,mach,time_stop])
            dat = numpy.array([])
            cmd = numpy.array(["cz","fn"])
            dof = numpy.array(["zg","vgnd","path"])
            cst = numpy.array([["altp","input[0]"],["altp_d","0"],["mach","input[1]"],["mach_d","0"],["path_d","0"]])
            stop = numpy.array(["time","<","input[2]"])

            flight,state = flight_segment(aircraft,MCR,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

            dtime = dt

    # Adjust last point using interpolation
    #-----------------------------------------------------------------------------------------------------------
    n = numpy.where(flight[:,stop_idx]<stop_val)[0][-1]

    xout = flight[n,:] + (flight[n+1,:]-flight[n,:]) \
                        *(stop_val - flight[n,stop_idx]) \
                        /(flight[n+1,stop_idx] - flight[n,stop_idx])

    flight = flight[0:n+2,:]
    flight[n+1,:] = xout
    state = xout[0:6]

    return flight,state


#===========================================================================================================
def standard_descent(aircraft,altp0,vcas1,altp1,vcas2,mach,flight):
    """
    Constant CAS climb :
    vcas1 from altp0 to altp1
    vcas2 from alpt1 to crossover altitude with iso mach
    """

    # Initialize climb flight
    #-----------------------------------------------------------------------------------------------------------
    xout = flight[-1,:]
    state = xout[0:6]
    altp2 = xout[xout_dict["altp"]]

    altp_cross_over = earth.cross_over_altp(vcas2,mach)

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    nei = 0.

    # Data filtering
    #-----------------------------------------------------------------------------------------------------------
    if(vcas1>vcas2):
        raise Exception("iso_cas_climb, vcas2 must be higher or equal than vcas1")

    if(altp0>altp1):
        raise Exception("iso_cas_climb, altp1 must be higher than starting altitude")

    if(altp1>altp2):
        raise Exception("iso_cas_climb, altp2 must be higher than altp1")

    if(altp1>altp_cross_over):
        raise Exception("iso_cas_climb, cross over altitude must be higher than altp1")

    # Constant mach descent segment
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([mach,altp_cross_over])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["vgnd","path"])
    cst = numpy.array([["mach","input[0]"],["mach_d","0"],["path_d","0"]])
    stop = numpy.array(["altp",">","input[1]"])
    dt = 30.

    flight,state = flight_segment(aircraft,FID,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    # First constant CAS descent segment
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([vcas2,altp1])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["vgnd","path"])
    cst = numpy.array([["vcas","input[0]"],["vcas_d","0"],["path_d","0"]])
    stop = numpy.array(["altp","<","input[1]"])
    dt = 30.

    flight,state = flight_segment(aircraft,FID,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    # Deceleration at constant altitude
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([altp1,vcas1])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["zg","path"])
    cst = numpy.array([["altp","input[0]"],["altp_d","0"],["path_d","0"]])
    stop = numpy.array(["vcas",">","input[1]"])
    dt = 5.

    flight,state = flight_segment(aircraft,FID,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    # Second constant CAS descent segment
    #-----------------------------------------------------------------------------------------------------------
    input = numpy.array([vcas1,altp0])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["vgnd","path"])
    cst = numpy.array([["vcas","input[0]"],["vcas_d","0"],["path_d","0"]])
    stop = numpy.array(["altp",">","input[1]"])
    dt = 30.

    flight,state = flight_segment(aircraft,FID,nei,state,stop,input,dat,cmd,dof,cst,dt,flight)

    return flight,state


#===========================================================================================================
def qs_mission(aircraft,disa,tow,range,heading,
               altp0_clb,vcas1_clb,altp1_clb,vcas2_clb,altp2_clb,cruise_mach,
               vcas2_dsc,altp1_dsc,vcas1_dsc,altp0_dsc,vz_mcr_min,vz_mcl_min):

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code
    nei = 0.

    altp_cabin = cabin_virtual_altp()

    t = 0.
    xg = 0.
    zg = earth.altg_from_altp(altp0_clb,disa)
    vgnd = earth.vtas_from_vcas(altp0_clb,disa,vcas1_clb)
    path = 0.
    state = numpy.array([t,tow,xg,zg,vgnd,path])

    xout,state_d = state_dot(xin,state,MCL,nei,aircraft)

    flight = numpy.vstack([[xout]])

    flight,state = iso_cas_climb(aircraft,vcas1_clb,altp1_clb,vcas2_clb,altp2_clb,cruise_mach,vz_mcl_min,flight)


    stop = numpy.array(["xg","<",str(range)])

    flight,state = level_flight(aircraft,cruise_mach,vz_mcl_min,vz_mcr_min,heading,stop,flight)


    flight,state = standard_descent(aircraft,altp0_dsc,vcas1_dsc,altp1_dsc,vcas2_dsc,cruise_mach,flight)






#===========================================================================================================
def init_var(aircraft,state,cmd,dof):

    n_cmd = numpy.shape(cmd)[0]
    n_dof = numpy.shape(dof)[0]

    var = [None]*(n_cmd+n_dof)

    param = {"cz":aircraft.aerodynamics.cz_cruise_lod_max,
             "fn":aircraft.propulsion.mcl_thrust_ref}

    for i in range(n_cmd):
        var[i] = param[cmd[i]]

    for i in range(n_dof):
        if (dof[i] == "path"):
            var[n_cmd+i] = 0.01     # Special init for path variable
        else:
            var[n_cmd+i] = state[state_dict[dof[i]]]

    return var


#===========================================================================================================
def test_ceilings(aircraft,nei,mach,altp,state):

    state_copy = copy.copy(state)

    (MTO,MCN,MCL,MCR,FID) = aircraft.propulsion.rating_code

    input = numpy.array([altp,mach])
    dat = numpy.array([["fn","fct_thrust(aircraft,state,rating,nei)"]])
    cmd = numpy.array(["cz"])
    dof = numpy.array(["zg","vgnd","path"])
    cst = numpy.array([["altp","input[0]"],["mach","input[1]"],["mach_d","0"],["path_d","0"]])

    var = init_var(aircraft,state,cmd,dof)

    xout,state,state_d = flight_point(aircraft,MCR,nei,state_copy,var,input,dat,cmd,dof,cst)
    vz_mcr = xout[xout_dict["vzair"]]

    var = init_var(aircraft,state,cmd,dof)

    xout,state,state_d = flight_point(aircraft,MCL,nei,state_copy,var,input,dat,cmd,dof,cst)
    vz_mcl = xout[xout_dict["vzair"]]

    return vz_mcr,vz_mcl


#===========================================================================================================
def test_air_range(aircraft,rating,nei,mach,altp,altp_up,state):

    state_copy = copy.copy(state)

    dat = numpy.array([])
    cmd = numpy.array(["cz","fn"])
    dof = numpy.array(["zg","vgnd","path"])
    cst = numpy.array([["altp","input[0]"],["altp_d","0"],["mach","input[1]"],["mach_d","0"],["path_d","0"]])

    var = init_var(aircraft,state,cmd,dof)

    input = numpy.array([altp,mach])
    xout,state,state_d = flight_point(aircraft,rating,nei,state_copy,var,input,dat,cmd,dof,cst)
    sar = - state_d[state_dict["xg"]] / state_d[state_dict["mass"]]

    var = init_var(aircraft,state,cmd,dof)

    input = numpy.array([altp_up,mach])
    xout_up,state_up,state_d_up = flight_point(aircraft,rating,nei,state_copy,var,input,dat,cmd,dof,cst)
    sar_up = - state_d_up[state_dict["xg"]] / state_d_up[state_dict["mass"]]

    return sar,sar_up

