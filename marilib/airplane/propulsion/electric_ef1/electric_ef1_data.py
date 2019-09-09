#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry

The EF1 architecture corresponds to the following features :

- Electrical main fans
- Eventual rear electrical fan  (attribute aircraft.ef1_power_elec_chain.rear_nacelle)
- Battery powered

IMPORTANT REMARKS :
Batteries are supposed to be stacked in the wing box in place of liquid fuel, so available volume is computed the same way
Three modes of battery stacking are possible :
1- "Variable"
   For each typical missions (nominal, max payload, max fuel, zero payload, cost) the mass of the battery
   which is put on board is exactly what is required for the mission, which supposes that the battery
   stacking system allows this ...
   In this case, battery mass is not accounted into OWE
2- "Max"
   Battery occupies the maximum available volume into the wing box. In this case, battery mass is fixed whatever the mission
   and is accounted into the OWE
3- "Fixed"
   Battery mass is fixed by the user through the attribute aircraft.weights.battery_in_owe
"""

#--------------------------------------------------------------------------------------------------------------------------------
class Ef1PowerElectricChain(object):
    """
    Electric chain data
    """
    INFO = {\
    "max_power":{"unit":"kW", "om":1.e4, "txt":"Maximum shaft power of optional rear engine"},
    "max_power_rating":{"unit":"int", "om":1.e0, "txt":"Rating of maximum shaft power of optional rear engine"},
    "overall_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Overall power efficiency of the electric chain"},
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
class ElectrofanPylon(object):
    """
    Electrofan pylon data
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
class ElectrofanNacelle(object):
    """
    Electrofan nacelle data
    """
    INFO = {\
    "n_engine":{"unit":"int", "om":1.e0, "txt":"Number of electric engine"},
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
    "motor_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Motor efficiency"},
    "controller_efficiency":{"unit":"no_dim", "om":1.e0, "txt":"Controller electric efficiency"},
    "controller_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of controller"},
    "motor_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of electric motor"},
    "nacelle_pw_density":{"unit":"kW/kg", "om":1.e0, "txt":"Power density of e-fan nacelle and mountings"},
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
        self.motor_efficiency = motor_efficiency
        self.controller_efficiency = controller_efficiency
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
class ElectrofanEngine(object):
    """
    Electric motor rating power in given conditions
    """
    INFO = {\
    "reference_thrust":{"unit":"daN", "om":1.e4, "txt":"Design Reference Thrust of main engines"},
    "reference_power":{"unit":"kW", "om":1.e3, "txt":"Design Reference Shaft Power of main engines"},
    "rating_factor":{"unit":"dict", "om":1.e0, "txt":"Dictinonary of rating factors versus reference thrust"},
    "mto_e_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in take off rating (one engine), Sea Level, ISA+15, Mach 0,25"},
    "mto_e_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in take off rating (one engine), Sea Level, ISA+15, Mach 0,25"},
    "mcn_e_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach"},
    "mcn_e_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach"},
    "mcl_e_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach"},
    "mcl_e_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach"},
    "mcr_e_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "mcr_e_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "fid_e_shaft_power":{"unit":"kW", "om":1.e3, "txt":"E-fan shaft power in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    "fid_e_fan_thrust":{"unit":"daN", "om":1.e3, "txt":"E-fan thrust in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach"},
    }
    def __init__(self, reference_thrust = None,
                       reference_power = None,
                       rating_factor = None,
                       mto_e_shaft_power = None,
                       mto_e_fan_thrust = None,
                       mcn_e_shaft_power = None,
                       mcn_e_fan_thrust = None,
                       mcl_e_shaft_power = None,
                       mcl_e_fan_thrust = None,
                       mcr_e_shaft_power = None,
                       mcr_e_fan_thrust = None,
                       fid_e_shaft_power = None,
                       fid_e_fan_thrust = None):
        self.reference_thrust = reference_thrust
        self.reference_power = reference_power
        self.rating_factor = rating_factor
        self.mto_e_shaft_power = mto_e_shaft_power
        self.mto_e_fan_thrust = mto_e_fan_thrust
        self.mcn_e_shaft_power = mcn_e_shaft_power
        self.mcn_e_fan_thrust = mcn_e_fan_thrust
        self.mcl_e_shaft_power = mcl_e_shaft_power
        self.mcl_e_fan_thrust = mcl_e_fan_thrust
        self.mcr_e_shaft_power = mcr_e_shaft_power
        self.mcr_e_fan_thrust = mcr_e_fan_thrust
        self.fid_e_shaft_power = fid_e_shaft_power
        self.fid_e_fan_thrust = fid_e_fan_thrust

#--------------------------------------------------------------------------------------------------------------------------------
class Ef1Battery(object):
    """
    Battery data
    """
    INFO = {\
    "stacking":{"unit":"string", "om":1.e0, "txt":"Battery mass strategy, can be Variable, Max or Given"},
    "energy_density":{"unit":"kWh/kg", "om":1.e0, "txt":"Battery energy density"},
    "power_density":{"unit":"kW/kg", "om":1.e0, "txt":"Battery power density (capability to release power per mass unit)"},
    "density":{"unit":"kg/m3", "om":1.e3, "txt":"Battery density (mass per volume unit)"},
    "mass_max":{"unit":"kg", "om":1.e3, "txt":"Maximum battery mass"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Global CG of batteries"}
    }
    def __init__(self, stacking = None,
                       energy_density = None,
                       power_density = None,
                       density = None,
                       mass_max = None,
                       c_g = None):
        self.stacking = stacking
        self.energy_density = energy_density
        self.power_density = power_density
        self.density = density
        self.mass_max = mass_max
        self.c_g = c_g

