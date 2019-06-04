#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class PowerElectricChain(object):
    """
    Electric chain data
    """
    def __init__(self, mto = None,
                       mcn = None,
                       mcl = None,
                       mcr = None,
                       fid = None,
                       max_power = None,
                       max_power_rating = None,
                       overall_efficiency = None,
                       generator_pw_density = None,
                       rectifier_pw_density = None,
                       wiring_pw_density = None,
                       cooling_pw_density = None,
                       mass = None,
                       c_g = None):
        """
            :param mto: Uu uc - OoM 10^0 - Take off power, mto<1: turbofan shaft power off take ratio, mto>1: e-fan motor power
            :param mcn: Uu uc - OoM 10^0 - Maxi continuous power, mcn<1: turbofan shaft power off take ratio, mcn>1: e-fan motor power
            :param mcl: Uu uc - OoM 10^0 - Max climb power, mcl<1: turbofan shaft power off take ratio, mcl>1: e-fan motor power
            :param mcr: Uu uc - OoM 10^0 - Max cruise power, mcr<1: turbofan shaft power off take ratio, mcr>1: e-fan motor power
            :param fid: Uu uc - OoM 10^0 - Flight idle power, fid<1: turbofan shaft power off take ratio, fid>1: e-fan motor power
            :param max_power: Uu kW - OoM 10^4 - E-fan motor maximum power
            :param max_power_rating: Uu int - OoM 10^0 - Engine rating of e-fan motor maximum power
            :param overall_efficiency: Uu no_dim - OoM 10^0 - Power efficiency of the electric chain
            :param generator_pw_density: Uu kW/kg - OoM 10^0 - Power density of electric generation
            :param rectifier_pw_density: Uu kW/kg - OoM 10^0 - Power density of rectifiers
            :param wiring_pw_density: Uu kW/kg - OoM 10^0 - Power density of wiring
            :param cooling_pw_density: Uu kW/kg - OoM 10^0 - Power density of cooling system
            :param mass: Uu kg - OoM 10^2 - Mass of the electric chain (generator, rectifier, wires, cooling)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the electric chain
        """
        self.mto = mto
        self.mcn = mcn
        self.mcl = mcl
        self.mcr = mcr
        self.fid = fid
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
class ElectricNacelle(object):
    """
    Electric nacelle data
    """
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
        """
        Constructor :
            :param width: Uu m - OoM 10^0 - Maximum width of the electric fan cowl
            :param length: Uu m - OoM 10^0 - Length of the electric fan cowl
            :param x_axe: Uu m - OoM 10^1 - Longitudinal position of the center of the electric nacelle air inlet
            :param y_axe: Uu m - OoM 10^1 - Span wise position of the center of the electric nacelle air inlet
            :param z_axe: Uu m - OoM 10^0 - Vertical position of the center of the electric nacelle air inlet
            :param net_wetted_area: Uu m2 - OoM 10^1 - Total net wetted area of the electric fan nacelle (fan cowl)
            :param efficiency_fan: Uu no_dim - OoM 10^0 - Fan efficiency for turbofan (capability to turn shaft power into kinetic energy)
            :param efficiency_prop: Uu no_dim - OoM 10^0 - "Propeller like" Fan+Cowl efficiency for turbofan (FanThrust.Vair)/(Shaft power)
            :param motor_efficiency: Uu no_dim - OoM 10^0 - Motor efficiency
            :param controller_efficiency: Uu no_dim - OoM 10^0 - Controller electric efficiency
            :param controller_pw_density: Uu kW/kg - OoM 10^0 - Power density of controller
            :param motor_pw_density: Uu kW/kg - OoM 10^0 - Power density of electric motor
            :param nacelle_pw_density: Uu kW/kg - OoM 10^0 - Power density of e-fan nacelle and mountings
            :param hub_width: Uu m - OoM 10^0 - Diameter of the hub of the electric nacelle
            :param fan_width: Uu m - OoM 10^0 - Diameter of the fan of the electric nacelle
            :param nozzle_width: Uu m - OoM 10^0 - Diameter of the nozzle of the electric nacelle
            :param nozzle_area: Uu m2 - OoM 10^0 - Exhaust nozzle area of the electric nacelle
            :param body_length: Uu m - OoM 10^0 - Length of the body behind the electric nacelle
            :param bnd_layer: Uu structure - OoM 10^0 - Boundary layer thickness law in front of the e-fan, 2d array
            :param mass: Uu kg - OoM 10^2 - Equipped mass of the nacelle of the electric fan (including the controller, motor and nacelle)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the electric nacelle
        """
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
class ElectricEngine(object):
    """
    Electric motor rating power in given conditions
    """
    def __init__(self, mto_e_power_ratio = None,
                       mto_e_shaft_power = None,
                       mto_e_fan_thrust = None,
                       mcn_e_power_ratio = None,
                       mcn_e_shaft_power = None,
                       mcn_e_fan_thrust = None,
                       mcl_e_power_ratio = None,
                       mcl_e_shaft_power = None,
                       mcl_e_fan_thrust = None,
                       mcr_e_power_ratio = None,
                       mcr_e_shaft_power = None,
                       mcr_e_fan_thrust = None,
                       fid_e_power_ratio = None,
                       fid_e_shaft_power = None,
                       fid_e_fan_thrust = None,
                       flight_data = None):
        """
        Constructor :
            :param mto_e_power_ratio: Uu no_dim - OoM 10^0 - Turbofan off take power ratio in take off rating (one engine), Sea Level, ISA+15, Mach 0,25
            :param mto_e_shaft_power: Uu kW - OoM 10^3 - E-fan shaft power in take off rating (one engine), Sea Level, ISA+15, Mach 0,25
            :param mto_e_fan_thrust : Uu daN - OoM 10^3 - E-fan thrust in take off rating (one engine), Sea Level, ISA+15, Mach 0,25
            :param mcn_e_power_ratio: Uu no_dim - OoM 10^0 - Turbofan off take power ratio in maxi continuous rating (one engine), required ceiling altitude, ISA, half cruise Mach
            :param mcn_e_shaft_power: Uu kW - OoM 10^3 - E-fan shaft power in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach
            :param mcn_e_fan_thrust : Uu daN - OoM 10^3 - E-fan thrust in maxi continuous rating (one engine), required ceiling altitude, ISA, cruise Mach
            :param mcl_e_power_ratio: Uu no_dim - OoM 10^0 - Turbofan off take power ratio in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach
            :param mcl_e_shaft_power: Uu kW - OoM 10^3 - E-fan shaft power in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach
            :param mcl_e_fan_thrust : Uu daN - OoM 10^3 - E-fan thrust in max climb rating (one engine), required Top of Climb altitude, ISA, cruise Mach
            :param mcr_e_power_ratio: Uu no_dim - OoM 10^0 - Turbofan off take power ratio in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param mcr_e_shaft_power: Uu kW - OoM 10^3 - E-fan shaft power in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param mcr_e_fan_thrust : Uu daN - OoM 10^3 - E-fan thrust in max cruise rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param fid_e_power_ratio: Uu no_dim - OoM 10^0 - Turbofan off take power ratio in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param fid_e_shaft_power: Uu kW - OoM 10^3 - E-fan shaft power in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param fid_e_fan_thrust : Uu daN - OoM 10^3 - E-fan thrust in flight idle rating (one engine), reference cruise altitude, ISA, cruise Mach
            :param flight_data : Uu dict - Dictionary of flying conditions for each rating {"disa":array, "altp":array, "mach":array, "nei":array}
        """
        self.mto_e_power_ratio = mto_e_power_ratio
        self.mto_e_shaft_power = mto_e_shaft_power
        self.mto_e_fan_thrust = mto_e_fan_thrust
        self.mcn_e_power_ratio = mcn_e_power_ratio
        self.mcn_e_shaft_power = mcn_e_shaft_power
        self.mcn_e_fan_thrust = mcn_e_fan_thrust
        self.mcl_e_power_ratio = mcl_e_power_ratio
        self.mcl_e_shaft_power = mcl_e_shaft_power
        self.mcl_e_fan_thrust = mcl_e_fan_thrust
        self.mcr_e_power_ratio = mcr_e_power_ratio
        self.mcr_e_shaft_power = mcr_e_shaft_power
        self.mcr_e_fan_thrust = mcr_e_fan_thrust
        self.fid_e_power_ratio = fid_e_power_ratio
        self.fid_e_shaft_power = fid_e_shaft_power
        self.fid_e_fan_thrust = fid_e_fan_thrust
        self.flight_data = flight_data

#--------------------------------------------------------------------------------------------------------------------------------
class Battery(object):
    """
    Battery data
    """
    def __init__(self, strategy = None,
                       power_feed = None,
                       time_feed = None,
                       energy_cruise = None,
                       energy_density = None,
                       power_density = None,
                       mass = None,
                       c_g = None):
        """
        Constructor :
            :param strategy: Uu int - OoM 10^0 - Battery sizing strategy, 1: power_feed & energy_cruise driven, 2: battery mass driven
            :param power_feed: Uu kW - OoM 10^4 - Power delivered to e-fan(s) at take off and(or) climb during a total of time_feed
            :param time_feed: Uu min - OoM 10^1 - Maximum duration of the power_feed delivered to e-fan(s)
            :param energy_cruise: Uu kWh - OoM 10^1 - Total battery energy dedicated to cruise
            :param energy_density: Uu kWh/kg - OoM 10^0 - Battery energy density
            :param power_density: Uu kW/kg - OoM 10^0 - Battery power density (capability to release power per mass unit
            :param mass: Uu kg - OoM 10^3 - Total battery mass
            :param c_g: Uu m - OoM 10^1 - Global CG of batteries
        """
        self.strategy = strategy
        self.power_feed = power_feed
        self.time_feed = time_feed
        self.energy_cruise = energy_cruise
        self.energy_density = energy_density
        self.power_density = power_density
        self.mass = mass
        self.c_g = c_g

