#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

from marilib import numpy

from marilib.earth import environment as earth


#===========================================================================================================
def air_density_grad(pamb,pamb_d,tamb,tamb_d):
    """
    Ideal gas density
    """
    R = earth.gas_constant()
    rho0 = earth.sea_level_density()
    rho = pamb/(R*tamb)
    rho_d = pamb_d/(R*tamb) - tamb_d*pamb/(R*tamb**2)
    sig = rho / rho0
    sig_d = rho_d / rho0
    return rho,rho_d, sig,sig_d


#===========================================================================================================
def sound_speed_grad(tamb,tamb_d):
    """
    Sound speed for ideal gas
    """
    R = earth.gas_constant()
    gam = earth.heat_ratio()
    vsnd = numpy.sqrt(gam*R*tamb)
    vsnd_d = 0.5*numpy.sqrt(gam*R/tamb)*tamb_d
    return vsnd,vsnd_d


#===========================================================================================================
def atmosphere_grad(altp,altp_d,disa,disa_d):
    """
    Pressure from pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gas_constant()

    Zs = 100.   # Smoothing altitude shift

    Z = numpy.array([0., 11000., 20000.,32000., 47000., 50000.])
    Z_d = numpy.zeros_like(Z)

    dtodz = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])
    dtodz_d = numpy.zeros_like(dtodz)

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    P_d = numpy.zeros_like(P)

    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])
    T_d = numpy.zeros_like(T)

    #===========================================================================================================
    def fct_atm_grad_t_p(j,T,T_d,P,P_d):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        T_d[j+1] = T_d[j] + dtodz_d[j]*(Z[j+1]-Z[j]) + dtodz[j]*(Z_d[j+1]-Z_d[j])
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(Z[1+j]-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(Z_d[1+j]-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(Z[1+j]-Z[j]))/(T[j]+disa) - (dtodz[j]*(Z[1+j]-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            P[1+j] = P[j]*B1**B2
            P_d[1+j] = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((Z[1+j]-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((Z_d[1+j]-Z_d[j])/(T[j]+disa)) + (g/R)*((Z[1+j]-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            P[j+1] = P[j]*numpy.exp(B3)
            P_d[j+1] = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        return T,T_d,P,P_d
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_atm_grad_pamb(j,altp,altp_d,disa,disa_d):
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(altp-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(altp_d-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(altp-Z[j]))/(T[j]+disa) - (dtodz[j]*(altp-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            pamb = P[j]*B1**B2
            pamb_d = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((altp-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((altp_d-Z_d[j])/(T[j]+disa)) + (g/R)*((altp-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            pamb = P[j]*numpy.exp(B3)
            pamb_d = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        return pamb,pamb_d
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_atm_grad_tamb(j,altp,altp_d,disa,disa_d):
        tamb = T[j] + dtodz[j]*(altp-Z[j]) + disa
        tamb_d = T_d[j] + dtodz_d[j]*(altp-Z[j]) + dtodz[j]*(altp_d-Z_d[j]) + disa_d
        return tamb,tamb_d
    #-----------------------------------------------------------------------------------------------------------

    if (Z[-1]<=altp):
        raise Exception("atmosphere_grad, altitude cannot exceed 50km")

    j = 0

    while (Z[1+j]-Zs<=altp):
        T,T_d,P,P_d = fct_atm_grad_t_p(j,T,T_d,P,P_d)
        j = j + 1

    if (Z[j]-Zs<=altp and altp<=Z[j]+Zs):
        t0,t0_d = fct_atm_grad_tamb(j-1,Z[j]-Zs,altp_d,disa,disa_d)
        t1,t1_d = fct_atm_grad_tamb(j-1,Z[j],altp_d,disa,disa_d)
        t2,t2_d = fct_atm_grad_tamb(j,Z[j]+Zs,altp_d,disa,disa_d)
        s = (altp-Z[j]+Zs)/(Zs*2.)
        s_d = altp_d/(Zs*2.)
        tamb = (1.-s)**2*t0 + 2.*s*(1.-s)*t1 + s**2*t2
        tamb_d = 2.*((t1-t0) + (t0-t1*2.+t2)*s)*s_d + (1.-s)**2*t0_d + 2.*s*(1.-s)*t1_d + s**2*t2_d
        dt_o_dz = ((t1-t0) + (t0-t1*2.+t2)*s)/Zs
        dt_o_dz_d = (t0-t1*2.+t2)/(2.*Zs**2)
        if (altp<=Z[j]):
            pamb,pamb_d = fct_atm_grad_pamb(j-1,altp,altp_d,disa,disa_d)
        else:
            pamb,pamb_d = fct_atm_grad_pamb(j,altp,altp_d,disa,disa_d)
    else:
        pamb,pamb_d = fct_atm_grad_pamb(j,altp,altp_d,disa,disa_d)
        tamb,tamb_d = fct_atm_grad_tamb(j,altp,altp_d,disa,disa_d)
        dt_o_dz = dtodz[j]
        dt_o_dz_d = dtodz_d[j]

    return pamb,pamb_d,tamb,tamb_d,dt_o_dz,dt_o_dz_d


#===========================================================================================================
def atmosphere_grad_old(altp,altp_d,disa,disa_d):
    """
    Pressure from pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gas_constant()

    Z = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    Z_d = numpy.zeros_like(Z)

    dtodz = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])
    dtodz_d = numpy.zeros_like(dtodz)

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    P_d = numpy.zeros_like(P)

    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])
    T_d = numpy.zeros_like(T)

    #===========================================================================================================
    def fct_atm_grad_t_p(j,T,T_d,P,P_d):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        T_d[j+1] = T_d[j] + dtodz_d[j]*(Z[j+1]-Z[j]) + dtodz[j]*(Z_d[j+1]-Z_d[j])
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(Z[1+j]-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(Z_d[1+j]-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(Z[1+j]-Z[j]))/(T[j]+disa) - (dtodz[j]*(Z[1+j]-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            P[1+j] = P[j]*B1**B2
            P_d[1+j] = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((Z[1+j]-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((Z_d[1+j]-Z_d[j])/(T[j]+disa)) + (g/R)*((Z[1+j]-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            P[j+1] = P[j]*numpy.exp(B3)
            P_d[j+1] = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        return T,T_d,P,P_d
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_atm_grad_amb(j,altp,altp_d,disa,disa_d):
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(altp-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(altp_d-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(altp-Z[j]))/(T[j]+disa) - (dtodz[j]*(altp-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            pamb = P[j]*B1**B2
            pamb_d = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((altp-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((altp_d-Z_d[j])/(T[j]+disa)) + (g/R)*((altp-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            pamb = P[j]*numpy.exp(B3)
            pamb_d = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        tamb = T[j] + dtodz[j]*(altp-Z[j]) + disa
        tamb_d = T_d[j] + dtodz_d[j]*(altp-Z[j]) + dtodz[j]*(altp_d-Z_d[j]) + disa_d
        return tamb,tamb_d,pamb,pamb_d
    #-----------------------------------------------------------------------------------------------------------

    if (Z[-1]<=altp):
        raise Exception("atmosphere_grad, altitude cannot exceed 50km")

    j = 0

    while (Z[1+j]<=altp):
        T,T_d,P,P_d = fct_atm_grad_t_p(j,T,T_d,P,P_d)
        j = j + 1

    tamb,tamb_d,pamb,pamb_d = fct_atm_grad_amb(j,altp,altp_d,disa,disa_d)

    return pamb,pamb_d,tamb,tamb_d,dtodz[j],dtodz_d[j]


#===========================================================================================================
def atmosphere_geo_grad(altg,altg_d,disa,disa_d):
    """
    Pressure from pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gas_constant()

    Zs = 100.

    Zi = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    dtodzi = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])

    Z = numpy.zeros_like(Zi)
    Z_d = numpy.zeros_like(Zi)

    dtodz = numpy.zeros_like(dtodzi)
    dtodz_d = numpy.zeros_like(dtodzi)

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    P_d = numpy.zeros_like(P)

    T = numpy.array([earth.sea_level_temperature()+disa, 0., 0., 0., 0., 0.])
    T_d = numpy.zeros_like(T)

    #===========================================================================================================
    def fct_atm_geo_grad_t_p(j,Z,Z_d,T,T_d,P,P_d):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        T_d[j+1] = T_d[j] + dtodz_d[j]*(Z[j+1]-Z[j]) + dtodz[j]*(Z_d[j+1]-Z_d[j])
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(Z[1+j]-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(Z_d[1+j]-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(Z[1+j]-Z[j]))/(T[j]+disa) - (dtodz[j]*(Z[1+j]-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            P[1+j] = P[j]*B1**B2
            P_d[1+j] = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((Z[1+j]-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((Z_d[1+j]-Z_d[j])/(T[j]+disa)) + (g/R)*((Z[1+j]-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            P[j+1] = P[j]*numpy.exp(B3)
            P_d[j+1] = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        K = 1 + disa/T[j+1]
        K_d = disa_d/T[j+1] - disa*T_d[j+1]/T[j+1]**2
        dtodz[j+1] = dtodzi[j+1]/K
        dtodz_d[j+1] = -dtodzi[j+1]*K_d/K**2
        Z[j+2] = Z[j+1] + (Zi[j+2]-Zi[j+1])*K
        Z_d[j+2] = Z_d[j+1] + (Zi[j+2]-Zi[j+1])*K_d
        return Z,Z_d,T,T_d,P,P_d,dtodz,dtodz_d
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_atm_geo_grad_pamb(j,altg,altg_d,disa,disa_d):
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(altg-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(altg_d-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(altg-Z[j]))/(T[j]+disa) - (dtodz[j]*(altg-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            pamb = P[j]*B1**B2
            pamb_d = P[j]*(B2*B1_dz*B1**(B2-1.) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((altg-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((altg_d-Z_d[j])/(T[j]+disa)) + (g/R)*((altg-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            pamb = P[j]*numpy.exp(B3)
            pamb_d = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        return pamb,pamb_d
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_atm_geo_grad_tamb(j,altg,altg_d,disa,disa_d):
        tamb = T[j] + dtodz[j]*(altg-Z[j]) + disa
        tamb_d = T_d[j] + dtodz_d[j]*(altg-Z[j]) + dtodz[j]*(altg_d-Z_d[j]) + disa_d
        return tamb,tamb_d
    #-----------------------------------------------------------------------------------------------------------

    K = 1 + disa/T[0]
    K_d = disa_d/T[0]
    dtodz[0] = dtodzi[0]/K
    dtodz_d[0] = -dtodzi[0]*K_d/K**2
    Z[1] = Z[0] + (Zi[1]-Zi[0])*K
    Z_d[1] = (Zi[1]-Zi[0])*K_d

    n = len(P)-1
    j = 0

    while (j<n and Z[1+j]-Zs<=altg):
        Z,Z_d,T,T_d,P,P_d,dtodz,dtodz_d = fct_atm_geo_grad_t_p(j,Z,Z_d,T,T_d,P,P_d)
        j = j + 1

    if (Z[1+j]<altg):
        raise Exception("atmosphere_geo_grad, altitude cannot exceed 50km")

    s = 0.

    if (Z[j]-Zs<=altg and altg<=Z[j]+Zs):
        t0,t0_d = fct_atm_geo_grad_tamb(j-1,Z[j]-Zs,altg_d,disa,disa_d)
        t1,t1_d = fct_atm_geo_grad_tamb(j-1,Z[j],altg_d,disa,disa_d)
        t2,t2_d = fct_atm_geo_grad_tamb(j,Z[j]+Zs,altg_d,disa,disa_d)
        s = (altg-Z[j]+Zs)/(Zs*2.)
        s_d = (altg_d - Z_d[j])/(Zs*2.)
        tamb = (1.-s)**2*t0 + 2.*s*(1.-s)*t1 + s**2*t2
        tamb_d = 2.*((t1-t0) + (t0-t1*2.+t2)*s)*s_d + (1.-s)**2*t0_d + 2.*s*(1.-s)*t1_d + s**2*t2_d
        dt_o_dz = ((t1-t0) + (t0-t1*2.+t2)*s)/Zs
        dt_o_dz_d = (t0-t1*2.+t2)/(2.*Zs**2)
        if (altg<=Z[j]):
            pamb,pamb_d = fct_atm_geo_grad_pamb(j-1,altg,altg_d,disa,disa_d)
        else:
            pamb,pamb_d = fct_atm_geo_grad_pamb(j,altg,altg_d,disa,disa_d)
    else:
        pamb,pamb_d = fct_atm_geo_grad_pamb(j,altg,altg_d,disa,disa_d)
        tamb,tamb_d = fct_atm_geo_grad_tamb(j,altg,altg_d,disa,disa_d)
        dt_o_dz = dtodz[j]
        dt_o_dz_d = dtodz_d[j]

    return pamb,pamb_d,tamb,tamb_d,dt_o_dz,dt_o_dz_d,s


#===========================================================================================================
def pressure_altitude_grad(pamb,pamb_d):
    """
    Pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gas_constant()

    Z = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    dtodz = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])

    #===========================================================================================================
    def fct_altp_grad_t_p(j,T,P):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        if (0.<numpy.abs(dtodz[j])):
            P[j+1] = P[j]*(1 + (dtodz[j]/T[j])*(Z[j+1]-Z[j]))**(-g/(R*dtodz[j]))
        else:
            P[j+1] = P[j]*numpy.exp(-(g/R)*((Z[j+1]-Z[j])/T[j]))
        return T,P
    #-----------------------------------------------------------------------------------------------------------

    #===========================================================================================================
    def fct_altp_grad_amb(j,pamb,T,P):
        if (0.<numpy.abs(dtodz[j])):
            altp = Z[j] + ((pamb/P[j])**(-(R*dtodz[j])/g) - 1.)*(T[j]/dtodz[j])
            altp_d = (T[j]/dtodz[j])*(-(R*dtodz[j])/(g*P[j]))*pamb_d*(pamb/P[j])**(-(R*dtodz[j])/g-1.)
        else:
            altp = Z[j] - (T[j]/(g/R))*numpy.log(pamb/P[j])
            altp_d =  - ((T[j]*R)/g)*(pamb_d/pamb)
        return altp,altp_d
    #-----------------------------------------------------------------------------------------------------------

    j = 0
    n = len(P)-1
    P[1] = P[0]*(1 + (dtodz[0]/T[0])*(Z[1]-Z[0]))**(-g/(R*dtodz[0]))
    T[1] = T[0] + dtodz[0]*(Z[1]-Z[0])

    while (j<n-1 and pamb<P[j+1]):
        j = j + 1
        T,P = fct_altp_grad_t_p(j,T,P)

    if (pamb<P[n]):
        raise Exception("pressure_altitude_grad, altitude cannot exceed 50km")

    altp,altp_d = fct_altp_grad_amb(j,pamb,T,P)

    return altp,altp_d


#===========================================================================================================
def vcas_from_mach_grad(pamb,pamb_d, mach,mach_d):
    """
    Calibrated air speed from Mach number, subsonic only
    """
    gam = earth.heat_ratio()
    P0 = earth.sea_level_pressure()
    vc0 = earth.sea_level_sound_speed()

    fac = gam/(gam-1.)
    Bp = pamb/P0
    Bp_d = pamb_d/P0
    Bm = ((gam-1.)/2.)*mach**2
    Bm_d = (gam-1.)*mach*mach_d
    B1 = Bp*((1.+Bm)**fac-1.)
    B1_dp = Bp_d*((1.+Bm)**fac-1.)
    B1_dm = Bp*Bm_d*fac*(1.+Bm)**(fac-1.)
    fac1 = 2./(gam-1.)
    B2 = fac1*((1.+B1)**(1./fac)-1.)
    B2_dp = fac1*B1_dp*(1./fac)*(1.+B1)**(1./fac-1.)
    B2_dm = fac1*B1_dm*(1./fac)*(1.+B1)**(1./fac-1.)

    vcas = vc0*numpy.sqrt(B2)
    vcas_d = (vc0/2.)*(B2_dp + B2_dm)/numpy.sqrt(B2)

    return vcas,vcas_d

