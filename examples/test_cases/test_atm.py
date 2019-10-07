#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib import numpy
from marilib.earth import environment_grad as earth

disa = 0.
disa_d = 1.

altg_d = 1.


Zg = numpy.linspace(10850.,11150.,31)

P = numpy.linspace(101325.,50.,2000)

for altg in Zg:
    pamb,pamb_d,tamb,tamb_d,dt_o_dz,dt_o_dz_d,s       = earth.atmosphere_geo_grad(altg,altg_d,disa,disa_d)
    pamb0,pamb_d0,tamb0,tamb_d0,dt_o_dz0,dt_o_dz_d0,S0 = earth.atmosphere_geo_grad(altg,altg_d,disa,0.)
    pamb1,pamb_d1,tamb1,tamb_d1,dt_o_dz1,dt_o_dz_d1,S1 = earth.atmosphere_geo_grad(altg,0.,disa,disa_d)
    pamb2,pamb2_d,tamb2,tamb2_d,dt_o_dz2,dt_o_dz_d2,S2 = earth.atmosphere_geo_grad(altg+0.0001,0.,disa,0.)
    pamb3,pamb3_d,tamb3,tamb3_d,dt_o_dz3,dt_o_dz_d3,S3 = earth.atmosphere_geo_grad(altg,0.,disa+0.1,0.)
    #print (pamb_d0, tamb_d1)
    print(altg, pamb_d0 - (pamb2-pamb0)/0.0001, tamb_d1 - (tamb3-tamb1)/0.1, dt_o_dz, dt_o_dz_d,s)










