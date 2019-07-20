#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry, DELMIRO Thales, GALLARD Francois

"""
from collections import OrderedDict

from marilib.tools import units as unit

from marilib.aircraft_data.operational_performances \
    import DesignDriver, LowSpeed, HighSpeed, MaxPayloadMission, \
           NominalMission, MaxFuelMission, ZeroPayloadMission, \
           CostMission, Economics, Environmental_Impact

from marilib.aircraft_data.physical_performances \
    import Aerodynamics, Propulsion, CharacteristicWeight, CenterOfGravity

from marilib.airplane.airframe.airframe_data \
    import Cabin, Payload, Fuselage, Wing, Tanks, LandingGears, \
           Systems, HorizontalTail, VerticalTail

from marilib.airplane.propulsion.turbofan.turbofan_data \
    import TurbofanPylon, TurbofanNacelle, TurbofanEngine

from marilib.airplane.propulsion.hybrid_pte1.hybrid_pte1_data \
    import Pte1PowerElectricChain, Pte1Battery, RearElectricNacelle, RearElectricEngine

from marilib.airplane.propulsion.electric_ef1.electric_ef1_data \
    import Ef1PowerElectricChain, Ef1Battery, ElectrofanPylon, ElectrofanNacelle, ElectrofanEngine


#--------------------------------------------------------------------------------------------------------------------------------
class Aircraft(object):
    """
    Assembling all aircraft data branches
    """
    def __init__(self, name = None):
        """
            Data structure branches, no ramification
        """
        self.name = name
        self.design_driver = DesignDriver()
        self.low_speed = LowSpeed()
        self.high_speed = HighSpeed()
        self.max_payload_mission = MaxPayloadMission()
        self.nominal_mission = NominalMission()
        self.max_fuel_mission = MaxFuelMission()
        self.zero_payload_mission = ZeroPayloadMission()
        self.cost_mission = CostMission()
        self.economics = Economics()
        self.environmental_impact = Environmental_Impact()

        self.aerodynamics = Aerodynamics()
        self.propulsion = Propulsion()
        self.weights = CharacteristicWeight()
        self.center_of_gravity = CenterOfGravity()

        self.cabin = Cabin()
        self.payload = Payload()
        self.fuselage = Fuselage()
        self.wing = Wing()
        self.landing_gears = LandingGears()
        self.horizontal_tail = HorizontalTail()
        self.vertical_tail = VerticalTail()
        self.tanks = Tanks()
        self.systems = Systems()

        self.turbofan_pylon = TurbofanPylon()
        self.turbofan_nacelle = TurbofanNacelle()
        self.turbofan_engine = TurbofanEngine()

        self.rear_electric_nacelle = RearElectricNacelle()
        self.rear_electric_engine = RearElectricEngine()
        self.pte1_power_elec_chain = Pte1PowerElectricChain()
        self.pte1_battery = Pte1Battery()

        self.electrofan_pylon = ElectrofanPylon()
        self.electrofan_nacelle = ElectrofanNacelle()
        self.electrofan_engine = ElectrofanEngine()
        self.ef1_power_elec_chain = Ef1PowerElectricChain()
        self.ef1_battery = Ef1Battery()

    def export_to_file(self, file="Aircraft.ini", def_order = True,\
                           user_format = True):
        """
        Build  ini file :
            Data tree file
            "out_file_path: File path - default value "Aircraft.ini"
            "def_order: parameters' order - default value True for class definition order
                                                  alternative value False for alphabetical order
        """
        from configobj import ConfigObj

        out_parser = ConfigObj(indent_type="    ")
        out_parser.filename = file

        if def_order: # class definition order
            data_dict = self.get_ordered_data_dict()
            write_ordered_data_dict(data_dict, 'Aircraft', out_parser, user_format)

        else: # alphabetical order
            data_dict = self.get_data_dict()
            write_data_dict(data_dict, 'Aircraft', out_parser, user_format)

        out_parser.write()


    def get_data_dict(self):
        return get_data_dict(self, "Aircraft", {})

    def get_ordered_data_dict(self):
        return get_ordered_data_dict(self, "Aircraft", OrderedDict())



#--------------------------------------------------------------------------------------------------------------------------------
def write_data_dict(data_dict, section, out_parser, user_format):

    for key in data_dict.keys():
        value = data_dict[key]

        if isinstance(value, dict):
            out_parser[key] = {}
            write_data_dict(value, key, out_parser[key], user_format)

        else:
            if user_format:
                out_parser[key] = unit.user_format(value)
            else:
                out_parser[key] = value

#--------------------------------------------------------------------------------------------------------------------------------
def write_ordered_data_dict(data_dict, section, out_parser, user_format):

    for key in data_dict.keys():
        value = data_dict[key]

        if isinstance(value, OrderedDict):
            out_parser[key] = OrderedDict()
            write_ordered_data_dict(value, key, out_parser[key], user_format)

        else:
            if user_format:
                out_parser[key] = unit.user_format(value)
            else:
                out_parser[key] = value

#--------------------------------------------------------------------------------------------------------------------------------
def is_basetype(obj):

    if obj is None or not hasattr(obj, "__dict__"):
        return True

    return False

#--------------------------------------------------------------------------------------------------------------------------------
def get_data_dict(obj, obj_name, data_dict):

    curr_data_d = {}
    data_dict[obj_name] = curr_data_d

    if not hasattr(obj, "__dict__"):
        return

    for attr_name in obj.__dict__.keys():
        attribute = getattr(obj, attr_name)

        if is_basetype(attribute):
            curr_data_d[attr_name] = attribute

        else:
            get_data_dict(attribute, attr_name, curr_data_d)

    return data_dict

#--------------------------------------------------------------------------------------------------------------------------------
def get_ordered_data_dict(obj, obj_name, ord_dict):
    from inspect import getsource
    ord_dict_d = OrderedDict()
    ord_dict[obj_name] = ord_dict_d
    if not hasattr(obj, "__dict__"):
        return
    l = getsource(obj.__class__)
    i = l.find("__init__")
    while i <= len(l):
        occur=l.find("self.",i)
        if occur == -1: # No attribute or end of class
            return ord_dict
        else:
            equal = l.find("=", occur)
            eol = l.find("\n", occur)
            if equal == -1 or equal >= eol: # no "=" inline
                return ord_dict
            else:
                spc = l.find(" ", occur,equal)
                if spc != -1:
                    attr_name = l[occur+5:spc]
                else:
                    attr_name = l[occur+5:equal]
                if attr_name != "info":
                    if attr_name in obj.__dict__.keys():
                        attribute = getattr(obj, attr_name)
                        if is_basetype(attribute):
                            ord_dict_d[attr_name] = attribute
                        else:
                            get_ordered_data_dict(attribute, attr_name, ord_dict_d)
                i = eol + 1
    return ord_dict


#===============================================================================
# main
#===============================================================================
if __name__ == "__main__":
    ac = Aircraft()
    ordered_dict = get_ordered_data_dict(ac, "Aircraft", OrderedDict())
    print(ordered_dict)
