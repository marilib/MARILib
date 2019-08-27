#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python

The PTE2 architecture corresponds to the following features :

- Tube & wing architecture with blimb body mounted on each wing
- Thermal main fans with electrical generators
- Rear electrical fan with use defined power on each rating
- Eventual additional battery (attribute aircraft.pte2_battery.battery_strategy)

IMPORTANT REMARKS :
If an additional battery is installed, two modes are available :
1- Battery mass is driven by the necessary amount of energy to ensure a power boost at take off and climb of "pte2_battery.power_feed" for a cumulated duration of "pte2_battery.time_feed"
   and(or) an additional energy "pte2_battery.energy_cruise" delivered all along the cruise
2- battery mass is given by the user (attribute aircraft.pte2_battery.mass)
"""

#--------------------------------------------------------------------------------------------------------------------------------
class Pte2BlimpBody(object):
    """
    Blimp body data
    """
    INFO = {\
    "width":{"unit":"m", "om":1.e0, "txt":"Diameter of the blimp bodies"},
    "length":{"unit":"m", "om":1.e1, "txt":"Length of the blimp bodies"},
    "x_axe":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the center of the blimp body nose"},
    "y_axe":{"unit":"m", "om":1.e1, "txt":"Span wise position of the center of the right blimp body nose"},
    "z_axe":{"unit":"m", "om":1.e0, "txt":"Vertical position of the center of the blimp body nose"},
    "net_wetted_area":{"unit":"m2", "om":1.e2, "txt":"Total net wetted area of the blimp bodies"},
    "usable_volume":{"unit":"m3", "om":1.e3, "txt":"Usable volume for helium into the blimp bodies"},
    "gas_type":{"unit":"string", "om":0., "txt":"type of buoyancy gas : helium or hydrogen"},
    "gas_mass":{"unit":"kg", "om":1.e2, "txt":"Mass of gas on board"},
    "buoyancy_force":{"unit":"daN", "om":1.e4, "txt":"Buoyancy force"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Equipped mass of the blimp bodies without engines"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the blimp bodies"}
    }
    def __init__(self, width = None,
                       length = None,
                       x_axe = None,
                       y_axe = None,
                       z_axe = None,
                       net_wetted_area = None,
                       usable_volume = None,
                       gas_type = None,
                       gas_mass = None,
                       buoyancy_force = None,
                       mass = None,
                       c_g = None):
        self.length = length
        self.width = width
        self.x_axe = x_axe
        self.y_axe = y_axe
        self.z_axe = z_axe
        self.net_wetted_area = net_wetted_area
        self.usable_volume = usable_volume
        self.gas_type = gas_type
        self.gas_mass = gas_mass
        self.buoyancy_force = buoyancy_force
        self.mass = mass
        self.c_g = c_g


#--------------------------------------------------------------------------------------------------------------------------------
class Pte2PowerElectricChain(object):
    """
    Electric chain data
    """
    INFO = {\
    "max_power":{"unit":"kW", "om":1.e4, "txt":"E-fan motor maximum power"},
    "max_power_rating":{"unit":"int", "om":1.e0, "txt":"Engine rating of e-fan motor maximum power"},
    "overall_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Power efficiency of the electric chain"},
    "generator_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of electric generation"},
    "rectifier_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of rectifiers"},
    "wiring_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of wiring"},
    "cooling_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of cooling system"},
    "mass":{"unit":"kg", "om":1.e2, "txt":"Mass of the electric chain (generator, rectifier, wires, cooling)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the electric chain"}
    }
    def __init__(self, max_power = None,
                       max_power_rating = None,
                       overall_efficiency = None,
                       generator_pw_density = None,
                       rectifier_pw_density = None,
                       wiring_pw_density = None,
                       cooling_pw_density = None,
                       mass = None,
                       c_g = None):
        self.max_power = max_power
        self.max_power_rating = max_power_rating
        self.overall_efficiency = overall_efficiency
        self.generator_pw_density = generator_pw_density
        self.rectifier_pw_density = rectifier_pw_density
        self.wiring_pw_density = wiring_pw_density
        self.cooling_pw_density = cooling_pw_density
        self.mass = mass
        self.c_g = c_g


#--------------------------------------------------------------------------------------------------------------------------------
class Pte2Battery(object):
    """
    Battery data
    """
    INFO = {\
    "strategy":{"unit":"int", "om":1.e0, "txt":"Battery sizing strategy, 0= no battery, 1= power_feed & energy_cruise driven, 2= battery mass driven"},
    "power_feed":{"unit":"kW", "om":1.e4, "txt":"Power delivered to e-fan(s) at take off and(or) climb during a total of time_feed"},
    "time_feed":{"unit":"min", "om":1.e1, "txt":"Maximum duration of the power_feed delivered to e-fan(s)"},
    "energy_cruise":{"unit":"kWh", "om":1.e1, "txt":"Total battery energy dedicated to cruise"},
    "energy_density":{"unit":"kWh/kg", "om":1.e0, "txt":"Battery energy density"},
    "power_density":{"unit":"kW/kg", "om":1.e0, "txt":"Battery power density (capability to release power per mass unit"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Total battery mass"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Global CG of batteries"}
    }
    def __init__(self, strategy = None,
                       power_feed = None,
                       time_feed = None,
                       energy_cruise = None,
                       energy_density = None,
                       power_density = None,
                       mass = None,
                       c_g = None):
        self.strategy = strategy
        self.power_feed = power_feed
        self.time_feed = time_feed
        self.energy_cruise = energy_cruise
        self.energy_density = energy_density
        self.power_density = power_density
        self.mass = mass
        self.c_g = c_g

