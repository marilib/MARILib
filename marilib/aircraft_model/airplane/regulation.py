#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""


from marilib.tools import units as unit


#===========================================================================================================
def kvs1g_min_take_off():
    k_vs1g_to = 1.13
    return k_vs1g_to

#===========================================================================================================
def kvs1g_min_landing():
    k_vs1g_ld = 1.23
    return k_vs1g_ld

#===========================================================================================================
def seg2_min_path(n_engine):
    """
    Regulatory min climb path versus number of engine
    """
    if(n_engine == 2):
        seg2_min_path = 0.024
    elif(n_engine == 3):
        seg2_min_path = 0.027
    elif(n_engine >= 4):
        seg2_min_path = 0.030
    return seg2_min_path

#===========================================================================================================
def ceil_oei_min_path(n_engine):
    """
    Regulatory min climb path depending on the number of engine
    """
    if(n_engine == 2):
        oei_min_path = 0.011
    elif(n_engine == 3):
        oei_min_path = 0.013
    elif(n_engine >= 4):
        oei_min_path = 0.016
    return oei_min_path

#===========================================================================================================
def diversion_range():
    diversion_range = unit.m_NM(200.)
    return diversion_range

#===========================================================================================================
def holding_time():
    holding_time = unit.s_min(30.)
    return holding_time

#===========================================================================================================
def reserve_fuel_ratio(design_range):
    if (design_range<unit.m_NM(6500.)):
        reserve_fuel_ratio = 0.05 ;
    else:
        reserve_fuel_ratio = 0.03 ;
    return reserve_fuel_ratio

#===========================================================================================================
def static_stability_margin():
    static_margin = 0.05
    return static_margin

