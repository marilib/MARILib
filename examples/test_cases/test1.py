#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

from marilib import numpy
from marilib.earth import environment_grad as earth

disa = 15.
disa_d = 1.

altg_d = 1.

pamb_d = 1.

Zp = numpy.linspace(0.,49999.,100)
P = numpy.linspace(101325.,50.,2000)

for pamb in P:
    altp,altp_d = earth.pressure_altitude_grad_2(pamb,pamb_d)
    altp2,altp_d2 = earth.pressure_altitude_grad_2(pamb,pamb_d)
    print(pamb, altp-altp2, altp_d-altp_d2)










