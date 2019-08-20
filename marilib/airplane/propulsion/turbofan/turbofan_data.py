#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python

The TF architecture corresponds to the following features :

- Twin or quad turbofan. In case of twin, possibility to attached engine on rear fuselage
- Tube & wing architecture with possibility to attach the wing on top of the fuselage
- Classical tail planes with the possibility to attach horizontal plane on top of the vertical tail
"""

#--------------------------------------------------------------------------------------------------------------------------------
class TurbofanPylon(object):
    """
    Turbofan pylon data
    """
    INFO = {\
    "mass":{"unit":"kg", "om":1.e3, "txt":"Equipped mass of the pylons"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the pylons"}
    }
    def __init__(self, mass = None,
                        c_g = None):
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbofanNacelle(object):
    """
    Turbofan nacelle data
    """
    INFO = {\
    "n_engine":{"unit":"int", "om":1.e0, "txt":"Number of turbofan"},
    "attachment":{"unit":"int", "om":1.e0, "txt":"Nacelle attachment (1= under wing, 2= rear fuselage)"},
    "rear_nacelle":{"unit":"int", "om":1.e0, "txt":"Rear nacelle (0= no, 1= yes)"},
    "width":{"unit":"m", "om":1.e0, "txt":"Maximum width of the nacelles"},
    "length":{"unit":"m", "om":1.e0, "txt":"Length of the fan cowl"},
    "x_ext":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the center of the air inlet of the external nacelle"},
    "y_ext":{"unit":"m", "om":1.e1, "txt":"Span wise position of the center of the air inlet of the external nacelle"},
    "z_ext":{"unit":"m", "om":1.e0, "txt":"Vertical position of the center of the air inlet of the external nacelle"},
    "x_int":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the center of the air inlet of the internal nacelle"},
    "y_int":{"unit":"m", "om":1.e1, "txt":"Span wise position of the center of the air inlet of the internal nacelle"},
    "z_int":{"unit":"m", "om":1.e0, "txt":"Vertical position of the center of the air inlet of the internal nacelle"},
    "net_wetted_area":{"unit":"m2", "om":1.e1, "txt":"Total net wetted area of the nacelles (fan cowls)"},
    "efficiency_fan":{"unit":"no_dim", "om":1.e0, "txt":"Fan efficiency for turbofan (capability to turn shaft power into kinetic energy)"},
    "efficiency_prop":{"unit":"no_dim", "om":1.e0, "txt":"Propeller like Fan+Cowl efficiency for turbofan (FanThrust.Vair)/(Shaft power)"},
    "hub_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the hub of the turbofan nacelle (for pusher fan only)"},
    "fan_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the fan of the turbofan nacelle"},
    "nozzle_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the nozzle of the turbofan nacelle"},
    "nozzle_area":{"unit":"m2", "om":1.e0, "txt":"Exhaust nozzle area of the turbofan nacelle"},
    "body_length":{"unit":"m", "om":1.e0, "txt":"Length of the body in front of the turbofan nacelle"},
    "bnd_layer":{"unit":"array", "om":1.e0, "txt":"Boundary layer thickness law in front of the e-fan, 2d array"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Equipped mass of the nacelles (including engine mass)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the nacelles"}
    }
    def __init__(self, n_engine = None,
                       attachment = None,
                       rear_nacelle = None,
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
        self.n_engine = n_engine
        self.attachment = attachment
        self.rear_nacelle = rear_nacelle
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
    INFO = {\
    "reference_thrust":{"unit":"daN", "om":1.e4, "txt":"Design Reference Thrust of the engines"},
    "bpr":{"unit":"no_dim", "om":1.e0, "txt":"By Pass Ratio of the turbofan"},
    "rating_factor":{"unit":"int", "om":1.e0, "txt":"Array of rating factors versus reference thrust"},
    "core_thrust_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Fraction of the total thrust of a turbofan which is due to the core (typically between 10% & 16% for BPR>5)"},
    "core_width_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Fraction of the total nacelle diameter which is taken by the core"},
    "core_weight_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Fraction of the total nacelle mass which is taken by the core"},
    "kfn_off_take":{"unit":"no_dim", "om":1.e0, "txt":"reference_thrust factor due to power off take (if any)"}
    }
    def __init__(self, reference_thrust = None,
                       bpr = None,
                       rating_factor = None,
                       core_thrust_ratio = None,
                       core_width_ratio = None,
                       core_weight_ratio = None,
                       kfn_off_take = None):
        self.reference_thrust = reference_thrust
        self.bpr = bpr
        self.rating_factor = rating_factor
        self.core_thrust_ratio = core_thrust_ratio
        self.core_width_ratio = core_width_ratio
        self.core_weight_ratio = core_weight_ratio
        self.kfn_off_take = kfn_off_take

