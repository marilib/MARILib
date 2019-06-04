#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class TurbofanPylon(object):
    """
    Turbofan pylon data
    """
    def __init__(self, mass = None,
                        c_g = None):
        """
        Constructor :
            :param mass: Uu kg - OoM 10^3 - Equipped mass of the pylons
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the pylons
        """
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbofanNacelle(object):
    """
    Turbofan nacelle data
    """
    def __init__(self, attachment = None,
                       width = None,
                       length = None,
                       x_ext = None,
                       y_ext = None,
                       z_ext = None,
                       x_int = None,
                       y_int = None,
                       z_int = None,
                       net_wetted_area = None,
                       efficiency_fan = None,
                       efficiency_prop = None,
                       hub_width = None,
                       fan_width = None,
                       nozzle_width = None,
                       nozzle_area = None,
                       body_length = None,
                       bnd_layer = None,
                       mass = None,
                       c_g = None):
        """
        Constructor

            :param attachment: Uu int - OoM 10^0 - Nacelle attachment (1= under wing, 2= rear fuselage)
            :param width: Uu m - OoM 10^0 - Maximum width of the nacelles
            :param length: Uu m - OoM 10^0 - Length of the fan cowl
            :param x_ext: Uu m - OoM 10^1 - Longitudinal position of the center of the air inlet of the external nacelle
            :param y_ext: Uu m - OoM 10^1 - Span wise position of the center of the air inlet of the external nacelle
            :param z_ext: Uu m - OoM 10^0 - Vertical position of the center of the air inlet of the external nacelle
            :param x_int: Uu m - OoM 10^1 - Longitudinal position of the center of the air inlet of the internal nacelle
            :param y_int: Uu m - OoM 10^1 - Span wise position of the center of the air inlet of the internal nacelle
            :param z_int: Uu m - OoM 10^0 - Vertical position of the center of the air inlet of the internal nacelle
            :param net_wetted_area: Uu m2 - OoM 10^1 - Total net wetted area of the nacelles (fan cowls)
            :param efficiency_fan: Uu no_dim - OoM 10^0 - Fan efficiency for turbofan (capability to turn shaft power into kinetic energy)
            :param efficiency_prop: Uu no_dim - OoM 10^0 - "Propeller like" Fan+Cowl efficiency for turbofan (FanThrust.Vair)/(Shaft power)
            :param hub_width: Uu m - OoM 10^0 - Diameter of the hub of the turbofan nacelle (for pusher fan only)
            :param fan_width: Uu m - OoM 10^0 - Diameter of the fan of the turbofan nacelle
            :param nozzle_width: Uu m - OoM 10^0 - Diameter of the nozzle of the turbofan nacelle
            :param nozzle_area: Uu m2 - OoM 10^0 - Exhaust nozzle area of the turbofan nacelle
            :param body_length: Uu m - OoM 10^0 - Length of the body in front of the turbofan nacelle
            :param bnd_layer: Uu structure - OoM 10^0 - Boundary layer thickness law in front of the e-fan, 2d array
            :param mass: Uu kg - OoM 10^3 - Equipped mass of the nacelles (including engine mass)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the nacelles
        """
        self.attachment = attachment
        self.width = width
        self.length = length
        self.x_ext = x_ext
        self.y_ext = y_ext
        self.z_ext = z_ext
        self.x_int = x_int
        self.y_int = y_int
        self.z_int = z_int
        self.net_wetted_area = net_wetted_area
        self.efficiency_fan = efficiency_fan
        self.efficiency_prop = efficiency_prop
        self.hub_width = hub_width
        self.fan_width = fan_width
        self.nozzle_width = nozzle_width
        self.nozzle_area = nozzle_area
        self.body_length = body_length
        self.bnd_layer = bnd_layer
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbofanEngine(object):
    """
    Turbofan engine data
    """
    def __init__(self, n_engine = None,
                       bpr = None,
                       reference_thrust = None,
                       rating_factor = None,
                       core_thrust_ratio = None,
                       core_width_ratio = None,
                       core_weight_ratio = None,
                       kfn_off_take = None):
        """
        Constructor :
            :param n_engine: Uu int - OoM 10^0 - Number of turbofan
            :param bpr: Uu no_dim - OoM 10^0 - By Pass Ratio of the turbofan
            :param reference_thrust: Uu daN - OoM 10^4 - Design Reference Thrust of the engines
            :param rating_factor: Uu int - OoM 10^0 - Array of rating factors versus reference thrust
            :param core_thrust_ratio: Uu no_dim - OoM 10^0 - Fraction of the total thrust of a turbofan which is due to the core (typically between 10% & 16% for BPR>5)
            :param core_width_ratio: Uu no_dim - OoM 10^0 - Fraction of the total nacelle diameter which is taken by the core
            :param core_weight_ratio: Uu no_dim - OoM 10^0 - Fraction of the total nacelle mass which is taken by the core
            :param kfn_off_take: Uu no_dim - OoM 10^0 - reference_thrust factor due to power off take (if any)
        """
        self.n_engine = n_engine
        self.bpr = bpr
        self.reference_thrust = reference_thrust
        self.rating_factor = rating_factor
        self.core_thrust_ratio = core_thrust_ratio
        self.core_width_ratio = core_width_ratio
        self.core_weight_ratio = core_weight_ratio
        self.kfn_off_take = kfn_off_take

