#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

import numpy

from marilib.earth import environment as earth


#===========================================================================================================
def air_density_grad(pamb,pamb_d,tamb,tamb_d):
    """
    Ideal gaz density
    """
    R = earth.gaz_constant()
    rho0 = earth.sea_level_density()
    rho = pamb/(R*tamb)
    rho_d = pamb_d/(R*tamb) - tamb_d*pamb/(R*tamb**2)
    sig = rho / rho0
    sig_d = rho_d / rho0
    return rho,rho_d, sig,sig_d


#===========================================================================================================
def sound_speed_grad(tamb,tamb_d):
    """
    Sound speed for ideal gaz
    """
    R = earth.gaz_constant()
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
    R = earth.gaz_constant()

    Z = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    Z_d = numpy.array([0., 0., 0.,0., 0., 0.])

    dtodz = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])
    dtodz_d = numpy.array([0., 0., 0., 0., 0.])

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    P_d = numpy.array([0., 0., 0., 0., 0., 0.])

    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])
    T_d = numpy.array([0., 0., 0., 0., 0., 0.])

    if (Z[-1]<altp):
        raise Exception("atmosphere_grad, altitude cannot exceed 50km")

    j = 0

    while (Z[1+j]<=altp):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        T_d[j+1] = T_d[j] + dtodz_d[j]*(Z[j+1]-Z[j]) + dtodz[j]*(Z_d[j+1]-Z_d[j])
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(Z[1+j]-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(Z_d[1+j]-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(Z[1+j]-Z[j]))/(T[j]+disa) - (dtodz[j]*(Z[1+j]-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            P[1+j] = P[j]*B1**B2
            P_d[1+j] = P[j]*(B2*B1_dz*B1**(B2-1) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((Z[1+j]-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((Z_d[1+j]-Z_d[j])/(T[j]+disa)) + (g/R)*((Z[1+j]-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            P[j+1] = P[j]*numpy.exp(B3)
            P_d[j+1] = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        j = j + 1

    if (0.<numpy.abs(dtodz[j])):
        B1 = 1 + (dtodz[j]*(altp-Z[j]))/(T[j]+disa)
        B1_dz = (dtodz[j]*(altp_d-Z_d[j]))/(T[j]+disa)
        B1_dt = (dtodz_d[j]*(altp-Z[j]))/(T[j]+disa) - (dtodz[j]*(altp-Z[j])*disa_d)/(T[j]+disa)**2
        B2 = -g/(R*dtodz[j])
        B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
        pamb = P[j]*B1**B2
        pamb_d = P[j]*(B2*B1_dz*B1**(B2-1) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
    else:
        B3 = -(g/R)*((altp-Z[j])/(T[j]+disa))
        B3_d = -(g/R)*((altp_d-Z_d[j])/(T[j]+disa)) + (g/R)*((altp-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
        pamb = P[j]*numpy.exp(B3)
        pamb_d = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)

    tamb = T[j] + dtodz[j]*(altp-Z[j]) + disa
    tamb_d = T_d[j] + dtodz_d[j]*(altp-Z[j]) + dtodz[j]*(altp_d-Z_d[j]) + disa_d

    return pamb,pamb_d,tamb,tamb_d,dtodz[j]


#===========================================================================================================
def atmosphere_geo_grad(altg,altg_d,disa,disa_d):
    """
    Pressure from pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gaz_constant()

    Zi = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    dtodzi = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])

    Z = numpy.array([0., 0., 0., 0., 0., 0.])
    Z_d = numpy.array([0., 0., 0.,0., 0., 0.])

    dtodz = numpy.array([0., 0., 0., 0., 0.])
    dtodz_d = numpy.array([0., 0., 0., 0., 0.])

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    P_d = numpy.array([0., 0., 0., 0., 0., 0.])

    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])
    T_d = numpy.array([0., 0., 0., 0., 0., 0.])

    K = 1 + disa/T[0]
    K_d = disa_d/T[0]
    dtodz[0] = dtodzi[0]/K
    dtodz_d[0] = -dtodzi[0]*K_d/K**2
    Z[1] = Z[0] + (Zi[1]-Zi[0])*K
    Z_d[1] = (Zi[1]-Zi[0])*K_d

    n = len(P)-1
    j = 0

    while (Z[1+j]<=altg):
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        T_d[j+1] = T_d[j] + dtodz_d[j]*(Z[j+1]-Z[j]) + dtodz[j]*(Z_d[j+1]-Z_d[j])
        if (0.<numpy.abs(dtodz[j])):
            B1 = 1 + (dtodz[j]*(Z[1+j]-Z[j]))/(T[j]+disa)
            B1_dz = (dtodz[j]*(Z_d[1+j]-Z_d[j]))/(T[j]+disa)
            B1_dt = (dtodz_d[j]*(Z[1+j]-Z[j]))/(T[j]+disa) - (dtodz[j]*(Z[1+j]-Z[j])*disa_d)/(T[j]+disa)**2
            B2 = -g/(R*dtodz[j])
            B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
            P[1+j] = P[j]*B1**B2
            P_d[1+j] = P[j]*(B2*B1_dz*B1**(B2-1) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
        else:
            B3 = -(g/R)*((Z[1+j]-Z[j])/(T[j]+disa))
            B3_d = -(g/R)*((Z_d[1+j]-Z_d[j])/(T[j]+disa)) + (g/R)*((Z[1+j]-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
            P[j+1] = P[j]*numpy.exp(B3)
            P_d[j+1] = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)
        j = j + 1
        K = 1 + disa/T[j]
        K_d = disa_d/T[j] - disa*T_d[j]/T[j]**2
        dtodz[j] = dtodzi[j]/K
        dtodz_d[j] = -dtodzi[j]*K_d/K**2
        Z[j+1] = Z[j] + (Zi[j+1]-Zi[j])*K
        Z_d[j+1] = Z_d[j] + (Zi[j+1]-Zi[j])*K_d

    if (Z[1+j]<altg):
        raise Exception("atmosphere_geo_grad, altitude cannot exceed 50km")

    if (0.<numpy.abs(dtodz[j])):
        B1 = 1 + (dtodz[j]*(altg-Z[j]))/(T[j]+disa)
        B1_dz = (dtodz[j]*(altg_d-Z_d[j]))/(T[j]+disa)
        B1_dt = (dtodz_d[j]*(altg-Z[j]))/(T[j]+disa) - (dtodz[j]*(altg-Z[j])*disa_d)/(T[j]+disa)**2
        B2 = -g/(R*dtodz[j])
        B2_dt = (g*dtodz_d[j])/(R*dtodz[j]**2)
        pamb = P[j]*B1**B2
        pamb_d = P[j]*(B2*B1_dz*B1**(B2-1) + (B1**B2)*(B2_dt*numpy.log(B1)+B2*B1_dt/B1)) + P_d[j]*B1**B2
    else:
        B3 = -(g/R)*((altg-Z[j])/(T[j]+disa))
        B3_d = -(g/R)*((altg_d-Z_d[j])/(T[j]+disa)) + (g/R)*((altg-Z[j])*(T_d[j]+disa_d)/(T[j]+disa)**2)
        pamb = P[j]*numpy.exp(B3)
        pamb_d = P_d[j]*numpy.exp(B3) + P[j]*B3_d*numpy.exp(B3)

    tamb = T[j] + dtodz[j]*(altg-Z[j]) + disa
    tamb_d = T_d[j] + dtodz_d[j]*(altg-Z[j]) + dtodz[j]*(altg_d-Z_d[j]) + disa_d

    return pamb,pamb_d,tamb,tamb_d,dtodz[j],dtodz_d[j]


#===========================================================================================================
def pressure_altitude_grad(pamb,pamb_d):
    """
    Pressure altitude from ground to 50 km
    """
    g = earth.gravity()
    R = earth.gaz_constant()

    Z = numpy.array([0., 10999., 19999.,31999., 46999., 49999.])
    dtodz = numpy.array([-0.0065, 0., 0.0010, 0.0028, 0.])

    P = numpy.array([earth.sea_level_pressure(), 0., 0., 0., 0., 0.])
    T = numpy.array([earth.sea_level_temperature(), 0., 0., 0., 0., 0.])

    j = 0
    n = len(P)-1
    P[1] = P[0]*(1 + (dtodz[0]/T[0])*(Z[1]-Z[0]))**(-g/(R*dtodz[0]))
    T[1] = T[0] + dtodz[0]*(Z[1]-Z[0])

    while (j<n and pamb<P[j+1]):
        j = j + 1
        T[j+1] = T[j] + dtodz[j]*(Z[j+1]-Z[j])
        if (0.<numpy.abs(dtodz[j])):
            P[j+1] = P[j]*(1 + (dtodz[j]/T[j])*(Z[j+1]-Z[j]))**(-g/(R*dtodz[j]))
        else:
            P[j+1] = P[j]*numpy.exp(-(g/R)*((Z[j+1]-Z[j])/T[j]))

    if (pamb<P[n]):
        raise Exception("pressure_altitude_grad, altitude cannot exceed 50km")

    if (0.<numpy.abs(dtodz[j])):
        altp = Z[j] + ((pamb/P[j])**(-(R*dtodz[j])/g) - 1)*(T[j]/dtodz[j])
        altp_d = (T[j]/dtodz[j])*(-(R*dtodz[j])/(g*P[j]))*pamb_d*(pamb/P[j])**(-(R*dtodz[j])/g-1)
    else:
        altp = Z[j] - (T[j]/(g/R))*numpy.log(pamb/P[j])
        altp_d =  - ((T[j]*R)/g)*(pamb_d/pamb)

    return altp,altp_d


#===========================================================================================================
def vcas_from_mach_grad(pamb,pamb_d, mach,mach_d):
    """
    Calibrated air speed from Mach number, subsonic only
    """
    gam = earth.heat_ratio()
    P0 = earth.sea_level_pressure()
    vc0 = earth.sea_level_sound_speed()

    fac = gam/(gam-1)
    Bp = pamb/P0
    Bp_d = pamb_d/P0
    Bm = ((gam-1)/2)*mach**2
    Bm_d = (gam-1)*mach*mach_d
    B1 = Bp*((1+Bm)**fac-1)
    B1_dp = Bp_d*((1+Bm)**fac-1)
    B1_dm = Bp*Bm_d*fac*(1+Bm)**(fac-1)
    fac1 = 2/(gam-1)
    B2 = fac1*((1+B1)**(1/fac)-1)
    B2_dp = fac1*B1_dp*(1/fac)*(1+B1)**(1/fac-1)
    B2_dm = fac1*B1_dm*(1/fac)*(1+B1)**(1/fac-1)

    vcas = vc0*numpy.sqrt(B2)
    vcas_d = (vc0/2)*(B2_dp + B2_dm)/numpy.sqrt(B2)

    return vcas,vcas_d

