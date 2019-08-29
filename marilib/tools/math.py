#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019
@author: DRUOT Thierry : original Scilab implementation
         ROCHES Pascal : portage to Python
"""

import warnings

from numpy.linalg import solve
from numpy.linalg.linalg import LinAlgError

from marilib import numpy
import marilib


# from autograd import jacobian
# from autograd.numpy.numpy_boxes import ArrayBox
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


def newton_solve(res_func, x0, dres_dy=None, args=(),
                 res_max=1e-12, relax=0.995, max_iter=50, full_output=False):
    if res_max <= 0:
        raise ValueError("res_max too small (%g <= 0)" % res_max)
    if max_iter < 1:
        raise ValueError("max_iter must be greater than 0")
    elif dres_dy is None:
        res_func, dres_dy = get_jac_func(res_func, args)

    msg = ""
    ier = 1
    k = 0
    y_curr = numpy.atleast_1d(1.0 * x0)
    myargs = (y_curr,) + args
    curr_res = res_func(*myargs)
    n0 = norm(curr_res)
    if n0 == 0.:
        return y_curr, n0, 1
    if not numpy.any(dres_dy(*myargs)):
        msg = "jacobian was zero."
        if marilib.is_driven_by_gems:
            warnings.warn(msg, RuntimeWarning)
        if full_output:
            infodict = {"stop_crit":n0, "niter":1}
            return y_curr, infodict, ier, msg
        else:
            return y_curr

    stop_crit = norm(curr_res) / n0
    while k < max_iter and stop_crit > res_max:
        drdy = dres_dy(*myargs)
        if marilib.is_using_autograd:
            if isinstance(drdy, ArrayBox):
                drdy = drdy._value
            if isinstance(curr_res, ArrayBox):
                curr_res = curr_res._value
        newt_step = solve(numpy.atleast_2d(drdy).astype("float64"),
                          curr_res.astype("float64"))
        y_curr -= relax * newt_step
        myargs = (y_curr,) + args
        curr_res = res_func(*myargs)
        stop_crit = norm(curr_res) / n0
        k += 1
    if k == max_iter and stop_crit > res_max:
        ier = 0
        msg = "Failed to converge Newton solver within " + str(max_iter) + " iterations "
        if marilib.is_driven_by_gems:
            raise LinAlgError(msg)

    if full_output:
        infodict = {"stop_crit":stop_crit, "niter":k}
        return y_curr, infodict, ier, msg
    else:
        return y_curr


def norm(x):
    return numpy.dot(x, x.T)**0.5


def approx_jac(res, y, args=(), step=1e-7):
    if not hasattr(y, "__len__"):
        y = numpy.array([y])
        n = 1
    else:
        n = len(y)
    jac = []
    myargs = (y,) + args
    res_ref = res(*myargs)
    pert = numpy.zeros_like(y)
    for i in xrange(n):
        yi = y[i]
        if hasattr(yi, "_value"):
            yi = yi._value
        pert[i] = step * 10**int(numpy.log10(abs(yi)))
        y += pert
        myargs_loc = (y,) + args
        res_pert = res(*myargs_loc)
        jac.append((res_pert - res_ref) / pert[i])
        y -= pert
        pert[i] = 0.
    if len(jac) == 1:
        return jac[0]
    jac = numpy.concatenate(jac)
    return jac


def get_jac_func(res, args=(), step=1e-6):
    if marilib.is_using_autograd:
        jac = jacobian(res, argnum=0)

        def j_func(*args):
            val = jac(*args)
            if isinstance(val, numpy.ndarray):
                return val
            if isinstance(val, ArrayBox):
                return val._value

        def my_res(*args):
            val = res(*args)
            if isinstance(val, numpy.ndarray):
                return val
            if isinstance(val, ArrayBox):
                return val._value
        return my_res, j_func

    def apprx(y, *args):
        return approx_jac(res, y, args, step)

    return res, apprx
