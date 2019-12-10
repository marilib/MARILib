#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

The TP architecture corresponds to the following features :

- Tube & wing architecture with default upper wing
- Twin or quad turboprop attached on the wing
- T-tail by default
"""

#--------------------------------------------------------------------------------------------------------------------------------
class TurbopropNacelle(object):
    """
    Turboprop nacelle data
    """
    INFO = {\
    "n_engine":{"unit":"int", "om":1.e0, "txt":"Number of turboprop"},
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
    "efficiency_prop":{"unit":"no_dim", "om":1.e0, "txt":"Propeller like Fan+Cowl efficiency for turboprop (FanThrust.Vair)/(Shaft power)"},
    "hub_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the hub of the turboprop nacelle (for pusher prop only)"},
    "propeller_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the propeller"},
    "body_length":{"unit":"m", "om":1.e0, "txt":"Length of the body in front of the turboprop nacelle"},
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
                       efficiency_prop = None,
                       hub_width = None,
                       propeller_width = None,
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
        self.efficiency_prop = efficiency_prop
        self.hub_width = hub_width
        self.propeller_width = propeller_width
        self.body_length = body_length
        self.bnd_layer = bnd_layer
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class TurbopropEngine(object):
    """
    Turboprop engine data
    """
    INFO = {\
    "reference_thrust":{"unit":"daN", "om":1.e4, "txt":"Design Reference Thrust of the engines"},
    "reference_power":{"unit":"kW", "om":1.e3, "txt":"Design reference power of the engine"},
    "rating_factor":{"unit":"int", "om":1.e0, "txt":"Array of rating factors versus reference thrust"},
    "kfn_off_take":{"unit":"no_dim", "om":1.e0, "txt":"reference_thrust factor due to power off take (if any)"}
    }
    def __init__(self, reference_thrust = None,
                       reference_power = None,
                       rating_factor = None,
                       kfn_off_take = None):
        self.reference_thrust = reference_thrust
        self.reference_power = reference_power
        self.rating_factor = rating_factor
        self.kfn_off_take = kfn_off_take

