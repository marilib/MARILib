#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import numpy

#===========================================================================================================
def lin_interp_1d(x,X,Y):
    """
    linear interpolation without any control
        x current position
        X array of the abscissa of the known points
        Y array of the known values at given abscissa
    """

    n = numpy.size(X)
    for j in range(1,n):
        if x<X[j] :
           y = Y[j-1]+(Y[j]-Y[j-1])*(x-X[j-1])/(X[j]-X[j-1])
           return y
    y = Y[n-2]+(Y[n-1]-Y[n-2])*(x-X[n-2])/(X[n-1]-X[n-2])
    return y


#===========================================================================================================
def vander3(X):
    """
    Return the vandermonde matrix of a dim 3 array
    A = [X^2, X, 1]
    """
    V = numpy.array([[X[0]**2, X[0], 1.],
                     [X[1]**2, X[1], 1.],
                     [X[2]**2, X[2], 1.]])
    return V


#===========================================================================================================
def trinome(A,Y):
    """
    calculates trinome coefficients from 3 given points
    A = [X2, X, 1] (Vandermonde matrix)
    """
    X = numpy.array([A[0][1], A[1][1], A[2][1]])
    X2 = numpy.array([A[0][0], A[1][0], A[2][0]])

    det = X2[0]*(X[1]-X[2])-X2[1]*(X[0]-X[2])+X2[2]*(X[0]-X[1])

    adet = Y[0]*(X[1]-X[2])-Y[1]*(X[0]-X[2])+Y[2]*(X[0]-X[1])

    bdet = X2[0]*(Y[1]-Y[2])-X2[1]*(Y[0]-Y[2])+X2[2]*(Y[0]-Y[1])

    cdet =  X2[0]*(X[1]*Y[2]-X[2]*Y[1])-X2[1]*(X[0]*Y[2]-X[2]*Y[0]) \
          + X2[2]*(X[0]*Y[1]-X[1]*Y[0])

    if det!=0:
        C = numpy.array([adet/det, bdet/det, cdet/det])
    elif X[0]!=X[2]:
        C = numpy.array([0., Y[0]-Y[2], Y[2]*X[0]-Y[0]*X[2]/(X[0]-X[2])])
    else:
        C = numpy.array([0., 0., (Y[0]+Y[1]+Y[2])/3.])

    return C


#===========================================================================================================
def maximize_1d(xini,dx,*fct):
    """
    Optimize 1 single variable, no constraint
    xini : initial value of the variable
    dx : fixed search step
    fct : function with the signature : ['function_name',a1,a2,a3,...,an]
          and function_name(x,a1,a2,a3,...,an)
    """

    n = len(fct[0])

    X0 = xini
    Y0 = fct[0][0](X0,*fct[0][1:n])

    X1 = X0+dx
    Y1 = fct[0][0](X1,*fct[0][1:n])

    if Y0>Y1:
        dx = -dx
        X0,X1 = X1,X0

    X2 = X1+dx
    Y2 = fct[0][0](X2,*fct[0][1:n])

    while Y1<Y2:
        X0 = X1
        X1 = X2
        X2 = X2+dx
        Y0 = Y1
        Y1 = Y2
        Y2 = fct[0][0](X2,*fct[0][1:n])

    X = numpy.array([X0,X1,X2])
    Y = numpy.array([Y0,Y1,Y2])

    A = vander3(X)     # [X**2, X, numpy.ones(3)]
    C = trinome(A,Y)

    xres = -C[1]/(2.*C[0])
    yres = fct[0][0](xres,*fct[0][1:n])

    rc = 1

    return (xres,yres,rc)


