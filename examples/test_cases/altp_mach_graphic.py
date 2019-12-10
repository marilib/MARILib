#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry
"""

import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

from marilib.tools import units as unit

tab = numpy.genfromtxt("scan_result_1.txt",delimiter = ";")[:,1:]

mach = list(set(tab[0,:]))
mach.sort()
altp = list(set(tab[1,:]))
altp.sort()

nx = len(mach)
ny = len(altp)

mtow = tab[5,:]
fuel = tab[17,:]
cost = tab[18,:]
feff = tab[19,:]

X, Y = np.meshgrid(altp, mach)

fig = plt.figure()
#ax = fig.gca(projection='3d')

fig.suptitle("Cruise Mach & Altp effect \n All airplanes optimized on same requirements \n 800NM evaluation mission", fontsize=14)

ax1 = fig.add_subplot(2, 2, 2, projection='3d')
ax1.set_title("Fuel (kg/trip)")
Z = fuel.reshape(nx,ny)
surf1 = ax1.plot_surface(X, Y, Z, cmap=cm.cool, linewidth=0, antialiased=False)

ax2 = fig.add_subplot(2, 2, 1, projection='3d')
ax2.set_title("MTOW (kg)")
Z = mtow.reshape(nx,ny)
surf2 = ax2.plot_surface(X, Y, Z, cmap=cm.cool, linewidth=0, antialiased=False)

ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax3.set_title("COC ($/trip)")
Z = cost.reshape(nx,ny)
surf3 = ax3.plot_surface(X, Y, Z, cmap=cm.cool, linewidth=0, antialiased=False)

ax4 = fig.add_subplot(2, 2, 4, projection='3d')
ax4.set_title("CO2 metric (10e-3kg/km/m0.48)")
Z = feff.reshape(nx,ny)
surf4 = ax4.plot_surface(X, Y, Z, cmap=cm.cool, linewidth=0, antialiased=False)


plt.show()
