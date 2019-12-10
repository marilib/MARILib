#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry
"""

#--------------------------------------------------------------------------------------------------------------------------------
class RearElectricNacelle(object):
    """
    Electric nacelle data
    """
    INFO = {\
    "width":{"unit":"m", "om":1.e0, "txt":"Maximum width of the electric fan cowl"},
    "length":{"unit":"m", "om":1.e0, "txt":"Length of the electric fan cowl"},
    "x_axe":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the center of the electric nacelle air inlet"},
    "y_axe":{"unit":"m", "om":1.e1, "txt":"Span wise position of the center of the electric nacelle air inlet"},
    "z_axe":{"unit":"m", "om":1.e0, "txt":"Vertical position of the center of the electric nacelle air inlet"},
    "net_wetted_area":{"unit":"m2", "om":1.e1, "txt":"Total net wetted area of the electric fan nacelle (fan cowl)"},
    "efficiency_fan":{"unit":"no_dim", "om":1.e0, "txt":"Fan efficiency for turbofan (capability to turn shaft power into kinetic energy)"},
    "efficiency_prop":{"unit":"no_dim", "om":1.e0, "txt":"Propeller like Fan+Cowl efficiency for turbofan (FanThrust.Vair)/(Shaft power)"},
    "motor_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Motor efficiency"},
    "controller_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Controller electric efficiency"},
    "controller_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of controller"},
    "motor_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of electric motor"},
    "nacelle_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of e-fan nacelle and mountings"},
    "hub_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the hub of the electric nacelle"},
    "fan_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the fan of the electric nacelle"},
    "nozzle_width":{"unit":"m", "om":1.e0, "txt":"Diameter of the nozzle of the electric nacelle"},
    "nozzle_area":{"unit":"m2", "om":1.e0, "txt":"Exhaust nozzle area of the electric nacelle"},
    "body_length":{"unit":"m", "om":1.e0, "txt":"Length of the body behind the electric nacelle"},
    "bnd_layer":{"unit":"array", "om":1.e0, "txt":"Boundary layer thickness law in front of the e-fan, 2d array"},
    "mass":{"unit":"kg", "om":1.e2, "txt":"Equipped mass of the nacelle of the electric fan (including the controller, motor and nacelle)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the electric nacelle"}
    }
    def __init__(self, width = None,
                       length = None,
                       x_axe = None,
                       y_axe = None,
                       z_axe = None,
                       net_wetted_area = None,
                       efficiency_fan = None,
                       efficiency_prop = None,
                       motor_efficiency = None,
                       controller_efficiency = None,
                       controller_pw_density = None,
                       motor_pw_density = None,
                       nacelle_pw_density = None,
                       hub_width = None,
                       fan_width = None,
                       nozzle_width = None,
                       nozzle_area = None,
                       body_length = None,
                       bnd_layer = None,
                       mass = None,
                       c_g = None):
        self.width = width
        self.length = length
        self.x_axe = x_axe
        self.y_axe = y_axe
        self.z_axe = z_axe
        self.net_wetted_area = net_wetted_area
        self.efficiency_fan = efficiency_fan
        self.efficiency_prop = efficiency_prop
        self.motor_efficiency = motor_efficiency
        self.controller_efficiency = controller_efficiency
        self.controller_pw_density = controller_pw_density
        self.motor_pw_density = motor_pw_density
        self.nacelle_pw_density = nacelle_pw_density
        self.hub_width = hub_width
        self.fan_width = fan_width
        self.nozzle_width = nozzle_width
        self.nozzle_area = nozzle_area
        self.body_length = body_length
        self.bnd_layer = bnd_layer
        self.mass = mass
        self.c_g = c_g


#--------------------------------------------------------------------------------------------------------------------------------
class RearElectricEngine(object):
    """
    Electric motor rating power in given conditions
    """
    INFO = {\
    "mto_r_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in take off rating (one engine), Sea Level, ISA+15, Mach 0,25"},
    "mto_r_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in take off rating (one engine), Sea Level, ISA+15, Mach 0,25"},
    "mcn_r_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach"},
    "mcn_r_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach"},
    "mcl_r_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach"},
    "mcl_r_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach"},
    "mcr_r_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "mcr_r_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "fid_r_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "fid_r_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    }
    def __init__(self, mto_r_shaft_power = None,
                       mto_r_fan_thrust = None,
                       mcn_r_shaft_power = None,
                       mcn_r_fan_thrust = None,
                       mcl_r_shaft_power = None,
                       mcl_r_fan_thrust = None,
                       mcr_r_shaft_power = None,
                       mcr_r_fan_thrust = None,
                       fid_r_shaft_power = None,
                       fid_r_fan_thrust = None):
        self.mto_r_shaft_power = mto_r_shaft_power
        self.mto_r_fan_thrust = mto_r_fan_thrust
        self.mcn_r_shaft_power = mcn_r_shaft_power
        self.mcn_r_fan_thrust = mcn_r_fan_thrust
        self.mcl_r_shaft_power = mcl_r_shaft_power
        self.mcl_r_fan_thrust = mcl_r_fan_thrust
        self.mcr_r_shaft_power = mcr_r_shaft_power
        self.mcr_r_fan_thrust = mcr_r_fan_thrust
        self.fid_r_shaft_power = fid_r_shaft_power
        self.fid_r_fan_thrust = fid_r_fan_thrust

