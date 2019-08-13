#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class Cabin(object):
    """
    Cabin data
    """
    INFO = {\
    "n_pax_ref":{"unit":"int", "om":1.e2, "txt":"Reference Number of passengers (often 2 class layout)"},
    "n_aisle":{"unit":"int", "om":1.e0, "txt":"Number of aisle in economic section"},
    "n_pax_front":{"unit":"int", "om":1.e0, "txt":"Number of seats in a row in economic section"},
    "fwd_limit":{"unit":"m", "om":1.e0, "txt":"Distance between aircraft nose and cabin forward limit"},
    "width":{"unit":"m", "om":1.e0, "txt":"Maximum width of the cabin (volume width, not floor width)"},
    "length":{"unit":"m", "om":1.e1, "txt":"Total length of the cabin"},
    "projected_area":{"unit":"m2", "om":1.e2, "txt":"Area of the cabin taking into account its maximum width (not floor area)"},
    "m_furnishing":{"unit":"kg", "om":1.e3, "txt":"Total mass of furnishing equipements"},
    "m_op_item":{"unit":"kg", "om":1.e3, "txt":"Total mass of operator items"},
    "cg_furnishing":{"unit":"m", "om":1.e1, "txt":"Center of gravity of furnishing equipements"},
    "cg_op_item":{"unit":"m", "om":1.e1, "txt":"Center of gravity of operator items"}
    }
    def __init__(self, n_pax_ref = None,
                       n_aisle = None,
                       n_pax_front = None,
                       fwd_limit = None,
                       width = None,
                       length = None,
                       projected_area = None,
                       m_furnishing = None,
                       m_op_item = None,
                       cg_furnishing = None,
                       cg_op_item = None):
        self.n_pax_ref = n_pax_ref
        self.n_aisle = n_aisle
        self.n_pax_front = n_pax_front
        self.fwd_limit = fwd_limit
        self.width = width
        self.length = length
        self.projected_area = projected_area
        self.m_furnishing = m_furnishing
        self.m_op_item = m_op_item
        self.cg_furnishing = cg_furnishing
        self.cg_op_item = cg_op_item

#--------------------------------------------------------------------------------------------------------------------------------
class Payload(object):
    """
    Payload CG range
    """
    INFO = {\
    "m_pax_nominal":{"unit":"kg", "om":1.e2, "txt":"Mass allowance per passenger to compute nominal payload"},
    "m_pax_max":{"unit":"kg", "om":1.e2, "txt":"Mass allowance per passenger to compute maximum payload"},
    "m_container_pallet":{"unit":"kg", "om":1.e3, "txt":"Mass of containers or pallets empty"},
    "nominal":{"unit":"kg", "om":1.e4, "txt":"Mass of nominal payload"},
    "maximum":{"unit":"kg", "om":1.e4, "txt":"Mass of maximum payload"},
    "max_fwd_mass":{"unit":"kg", "om":1.e2, "txt":"Payload mass at maximum forward payload CG"},
    "max_fwd_req_cg":{"unit":"m", "om":1.e1, "txt":"Required maximum forward payload CG"},
    "max_bwd_mass":{"unit":"kg", "om":1.e2, "txt":"Payload mass at maximum backward payload CG"},
    "max_bwd_req_cg":{"unit":"m", "om":1.e1, "txt":"Required maximum backward payload CG"},
    "cg_container_pallet":{"unit":"m", "om":1.e1, "txt":"Center of gravity of containers or pallets empty"}
    }
    def __init__(self, m_pax_nominal = None,
                       m_pax_max = None,
                       m_container_pallet = None,
                       nominal = None,
                       maximum = None,
                       max_fwd_mass = None,
                       max_fwd_req_cg = None,
                       max_bwd_mass = None,
                       max_bwd_req_cg = None,
                       cg_container_pallet = None):
        self.m_pax_nominal = m_pax_nominal
        self.m_pax_max = m_pax_max
        self.m_container_pallet = m_container_pallet
        self.nominal = nominal
        self.maximum = maximum
        self.max_fwd_mass = max_fwd_mass
        self.max_fwd_req_cg = max_fwd_req_cg
        self.max_bwd_mass = max_bwd_mass
        self.max_bwd_req_cg = max_bwd_req_cg
        self.cg_container_pallet = cg_container_pallet

#--------------------------------------------------------------------------------------------------------------------------------
class Fuselage(object):
    """
    Fuselage data
    """
    INFO = {\
    "width":{"unit":"m", "om":1.e0, "txt":"Fuselage width of the cylindrical part"},
    "height":{"unit":"m", "om":1.e0, "txt":"Fuselage height of the cylindrical part"},
    "length":{"unit":"m", "om":1.e1, "txt":"Total fuselage length"},
    "tail_cone_length":{"unit":"m", "om":1.e0, "txt":"Length of rear evolutive part of the fuselage"},
    "net_wetted_area":{"unit":"m2", "om":1.e2, "txt":"Fuselage total net wetted area"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Equipped fuselage mass (without systems)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the fuselage CG"}
    }
    def __init__(self, width = None,
                       height = None,
                       length = None,
                       tail_cone_length = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None):
        self.width = width
        self.height = height
        self.length = length
        self.tail_cone_length = tail_cone_length
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class Wing(object):
    """
    Wing data
    """
    INFO = {\
    "attachment":{"unit":"int", "om":1.e0, "txt":"Wing attachment, 1: low wing, 2: high wing"},
    "morphing":{"unit":"int", "om":1.e0, "txt":"Wing deformation driver, 1: aspect ratio, 2: span"},
    "hld_type":{"unit":"int", "om":1.e0, "txt":"Type of high lift devices"},
    "t_o_c_r":{"unit":"%", "om":1.e1, "txt":"Thickness over chord ratio of the wing at root"},
    "t_o_c_k":{"unit":"%", "om":1.e1, "txt":"Thickness over chord ratio of the wing at main kink"},
    "t_o_c_t":{"unit":"%", "om":1.e1, "txt":"Thickness over chord ratio at wing tip"},
    "sweep":{"unit":"deg", "om":1.e0, "txt":"Wing sweep angle at 25% of the chord"},
    "dihedral":{"unit":"deg", "om":1.e0, "txt":"Mean dihedral of the wing"},
    "setting":{"unit":"deg", "om":1.e0, "txt":"Setting angle of the wing at root"},
    "taper_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Wing taper ratio"},
    "aspect_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Wing aspect ratio"},
    "area":{"unit":"m2", "om":1.e2, "txt":"Wing reference area"},
    "span":{"unit":"m", "om":1.e1, "txt":"Wing span"},
    "mac":{"unit":"m", "om":1.e0, "txt":"Mean aerodynamic chord of the wing"},
    "net_wetted_area":{"unit":"m2", "om":1.e2, "txt":"Wing total net wetted area"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Equipped wing mass (without systems)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the wing CG"},
    "x_root":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of 0% of wing root chord"},
    "y_root":{"unit":"m", "om":1.e0, "txt":"Span wise position of 0% of the wing root chord"},
    "z_root":{"unit":"m", "om":1.e0, "txt":"Vertical position of 0% of the wing root chord"},
    "c_root":{"unit":"m", "om":1.e0, "txt":"Wing root chord length"},
    "x_kink":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of 0% of wing kink chord"},
    "y_kink":{"unit":"m", "om":1.e0, "txt":"Span wise position of 0% of the wing kink chord"},
    "z_kink":{"unit":"m", "om":1.e0, "txt":"Vertical position of 0% of the wing kink chord"},
    "c_kink":{"unit":"m", "om":1.e0, "txt":"Wing kink chord length"},
    "x_tip":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of 0% of wing tip chord"},
    "y_tip":{"unit":"m", "om":1.e0, "txt":"Span wise position of 0% of the wing tip chord"},
    "z_tip":{"unit":"m", "om":1.e0, "txt":"Vertical position of 0% of the wing tip chord"},
    "c_tip":{"unit":"m", "om":1.e0, "txt":"Wing tip chord length"},
    "x_mac":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of wing mean aerodynamic chord"},
    "y_mac":{"unit":"m", "om":1.e0, "txt":"Span wise position of wing mean aerodynamic chord"}
    }
    def __init__(self, attachment = None,
                       morphing = None,
                       hld_type = None,
                       t_o_c_r = None,
                       t_o_c_k = None,
                       t_o_c_t = None,
                       sweep = None,
                       dihedral = None,
                       setting = None,
                       taper_ratio = None,
                       aspect_ratio = None,
                       area = None,
                       span = None,
                       mac = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None,
                       x_root = None,
                       y_root = None,
                       z_root = None,
                       c_root = None,
                       x_kink = None,
                       y_kink = None,
                       z_kink = None,
                       c_kink = None,
                       x_tip = None,
                       y_tip = None,
                       z_tip = None,
                       c_tip = None,
                       x_mac = None,
                       y_mac = None):
        self.attachment = attachment
        self.morphing = morphing
        self.hld_type = hld_type
        self.t_o_c_r = t_o_c_r
        self.t_o_c_k = t_o_c_k
        self.t_o_c_t = t_o_c_t
        self.sweep = sweep
        self.dihedral = dihedral
        self.setting = setting
        self.taper_ratio = taper_ratio
        self.aspect_ratio = aspect_ratio
        self.area = area
        self.span = span
        self.mac = mac
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g
        self.x_root = x_root
        self.y_root = y_root
        self.z_root = z_root
        self.c_root = c_root
        self.x_kink = x_kink
        self.y_kink = y_kink
        self.z_kink = z_kink
        self.c_kink = c_kink
        self.x_tip = x_tip
        self.y_tip = y_tip
        self.z_tip = z_tip
        self.c_tip = c_tip
        self.x_mac = x_mac
        self.y_mac = y_mac

#--------------------------------------------------------------------------------------------------------------------------------
class LandingGears(object):
    """
    Landing gear data
    """
    INFO = {\
    "mass":{"unit":"kg", "om":1.e3, "txt":"Mass of landing gears (nose and main)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the landing gears CG"}
    }
    def __init__(self, mass = None,
                        c_g = None):
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class HorizontalTail(object):
    """
    Horizontal tail plane data
    """
    INFO = {\
    "attachment":{"unit":"int", "om":1.e0, "txt":"Configuration of horizontal tail, 1: classical, 2: T-tail"},
    "sweep":{"unit":"deg", "om":1.e0, "txt":"Horizontal tail sweep angle at 25% of the chords"},
    "taper_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Taper ratio of the horizontal tail"},
    "aspect_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Aspect ratio of the horizontal tail"},
    "t_o_c":{"unit":"no_dim", "om":1.e0, "txt":"Thickness to chord ratio of the vertical tail"},
    "dihedral":{"unit":"deg", "om":1.e0, "txt":"Mean dihedral of the horizontal tail"},
    "volume":{"unit":"no_dim", "om":1.e0, "txt":"Volume coefficient of the  horizontal tail"},
    "lever_arm":{"unit":"m", "om":1.e1, "txt":"Lever arm of the horizontal tail (from 25% wing MAC to 25% HTTP MAC"},
    "area":{"unit":"m2", "om":1.e2, "txt":"Horizontal tail reference area"},
    "span":{"unit":"m", "om":1.e1, "txt":"Horizontal tail span"},
    "mac":{"unit":"m", "om":1.e0, "txt":"Mean aerodynamic part of the horizontal tail"},
    "net_wetted_area":{"unit":"m2", "om":1.e2, "txt":"Total net wetted area of the horizontal tail"},
    "mass":{"unit":"kg", "om":1.e2, "txt":"Equipped mass of the horizontal tail"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the horizontal tail"},
    "x_axe":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the horizontal tail central chord"},
    "z_axe":{"unit":"m", "om":1.e1, "txt":"Vertical position of the horizontal tail central chord"},
    "c_axe":{"unit":"m", "om":1.e1, "txt":"Horizontal tail central chord"},
    "x_tip":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the horizontal tail tip chord"},
    "y_tip":{"unit":"m", "om":1.e1, "txt":"Lateral position of the horizontal tail tip chord"},
    "z_tip":{"unit":"m", "om":1.e1, "txt":"Vertical position of the horizontal tail tip chord"},
    "c_tip":{"unit":"m", "om":1.e1, "txt":"Horizontal tail tip chord"},
    "x_mac":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the horizontal tail mean aerodynamic chord"},
    "y_mac":{"unit":"m", "om":1.e1, "txt":"Lateral position of the horizontal tail mean chord"}
    }
    def __init__(self, attachment = None,
                       sweep = None,
                       taper_ratio = None,
                       aspect_ratio = None,
                       t_o_c = None,
                       dihedral = None,
                       volume = None,
                       lever_arm = None,
                       area = None,
                       span = None,
                       mac = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None,
                       x_axe = None,
                       z_axe = None,
                       c_axe = None,
                       x_tip = None,
                       y_tip = None,
                       z_tip = None,
                       c_tip = None,
                       x_mac = None,
                       y_mac = None):
        self.attachment = attachment
        self.sweep = sweep
        self.taper_ratio = taper_ratio
        self.aspect_ratio = aspect_ratio
        self.t_o_c = t_o_c
        self.dihedral = dihedral
        self.volume = volume
        self.lever_arm = lever_arm
        self.area = area
        self.span = span
        self.mac = mac
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g
        self.x_axe = x_axe
        self.z_axe = z_axe
        self.c_axe = c_axe
        self.x_tip = x_tip
        self.y_tip = y_tip
        self.z_tip = z_tip
        self.c_tip = c_tip
        self.x_mac = x_mac
        self.y_mac = y_mac

#--------------------------------------------------------------------------------------------------------------------------------
class VerticalTail(object):
    """
    Vertical tail plane data
    """
    INFO = {\
    "sweep":{"unit":"deg", "om":1.e0, "txt":"Vertical tail sweep angle at 25% of the chords"},
    "taper_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Taper ratio of the vertical tail"},
    "aspect_ratio":{"unit":"no_dim", "om":1.e0, "txt":"Aspect ratio of the vertical tail"},
    "t_o_c":{"unit":"no_dim", "om":1.e0, "txt":"Thickness to chord ratio of the vertical tail"},
    "volume":{"unit":"m2/kN", "om":1.e0, "txt":"Volume coefficient of the  vertical tail"},
    "lever_arm":{"unit":"m", "om":1.e1, "txt":"Lever arm of the vertical tail (from 25% wing MAC to 25% HTTP MAC"},
    "area":{"unit":"m2", "om":1.e2, "txt":"Vertical tail reference area"},
    "height":{"unit":"m2", "om":1.e2, "txt":"Vertical tail height"},
    "mac":{"unit":"m", "om":1.e0, "txt":"Mean aerodynamic part of the vertical tail"},
    "net_wetted_area":{"unit":"m2", "om":1.e2, "txt":"Total net wetted area of the vertical tail"},
    "mass":{"unit":"kg", "om":1.e2, "txt":"Equipped mass of the vertical tail"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the CG of the vertical tail"},
    "x_root":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the vertical tail root chord"},
    "z_root":{"unit":"m", "om":1.e1, "txt":"Vertical position of the vertical tail root chord"},
    "c_root":{"unit":"m", "om":1.e1, "txt":"Vertical tail root chord"},
    "x_tip":{"unit":"m", "om":1.e2, "txt":"Longitudinal position of the vertical tail tip chord"},
    "z_tip":{"unit":"m", "om":1.e1, "txt":"Vertical position of the vertical tail tip chord"},
    "c_tip":{"unit":"m", "om":1.e1, "txt":"Vertical tail tip chord"},
    "x_mac":{"unit":"m", "om":1.e2, "txt":"Longitudinal position of the vertical tail mean aerodynamic chord"}
    }
    def __init__(self, sweep = None,
                       taper_ratio = None,
                       aspect_ratio = None,
                       t_o_c = None,
                       volume = None,
                       lever_arm = None,
                       area = None,
                       height =None,
                       mac = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None,
                       x_root = None,
                       z_root = None,
                       c_root = None,
                       x_tip = None,
                       z_tip = None,
                       c_tip = None,
                       x_mac = None):
        self.sweep = sweep
        self.taper_ratio = taper_ratio
        self.aspect_ratio = aspect_ratio
        self.t_o_c = t_o_c
        self.volume = volume
        self.lever_arm = lever_arm
        self.area = area
        self.height = height
        self.mac = mac
        self.net_wetted_area = net_wetted_area
        self.mass = mass
        self.c_g = c_g
        self.x_root = x_root
        self.z_root = z_root
        self.c_root = c_root
        self.x_tip = x_tip
        self.z_tip = z_tip
        self.c_tip = c_tip
        self.x_mac = x_mac

#--------------------------------------------------------------------------------------------------------------------------------
class Tanks(object):
    """
    Fuel volume data
    """
    INFO = {\
    "cantilever_volume":{"unit":"m3", "om":1.e1, "txt":"Volume of tanks in the cantilever wing"},
    "central_volume":{"unit":"m3", "om":1.e1, "txt":"Volume of tanks in the central part of the wing (inside the fuselage)"},
    "mfw_volume_limited":{"unit":"kg", "om":1.e4, "txt":"Maximum geometrical fuel volume"},
    "fuel_density":{"unit":"kg/m3", "om":1.e2, "txt":"Fuel density"},
    "fuel_cantilever_cg":{"unit":"m", "om":1.e1, "txt":"Center of gravity of tanks in the cantilever wing"},
    "fuel_central_cg":{"unit":"m", "om":1.e1, "txt":"Center of gravity of tanks in the central part of the wing (inside the fuselage)"},
    "fuel_total_cg":{"unit":"m", "om":1.e1, "txt":"Center of gravity of wing tanks"},
    "fuel_max_fwd_mass":{"unit":"kg", "om":1.e3, "txt":"Fuel mass of max forward fuel cg"},
    "fuel_max_fwd_cg":{"unit":"m", "om":1.e1, "txt":"Max forward fuel cg"},
    "fuel_max_bwd_mass":{"unit":"kg", "om":1.e3, "txt":"Fuel mass of max backward fuel cg"},
    "fuel_max_bwd_cg":{"unit":"m", "om":1.e1, "txt":"Max backward fuel cg"}
    }
    def __init__(self, cantilever_volume = None,
                       central_volume = None,
                       mfw_volume_limited = None,
                       fuel_density = None,
                       fuel_cantilever_cg = None,
                       fuel_central_cg = None,
                       fuel_total_cg = None,
                       fuel_max_fwd_mass = None,
                       fuel_max_fwd_cg = None,
                       fuel_max_bwd_mass = None,
                       fuel_max_bwd_cg = None,
                 ):
        self.cantilever_volume = cantilever_volume
        self.central_volume = central_volume
        self.mfw_volume_limited = mfw_volume_limited
        self.fuel_density = fuel_density
        self.fuel_cantilever_cg = fuel_cantilever_cg
        self.fuel_central_cg = fuel_central_cg
        self.fuel_total_cg = fuel_total_cg
        self.fuel_max_fwd_mass = fuel_max_fwd_mass
        self.fuel_max_fwd_cg = fuel_max_fwd_cg
        self.fuel_max_bwd_mass = fuel_max_bwd_mass
        self.fuel_max_bwd_cg = fuel_max_bwd_cg

#--------------------------------------------------------------------------------------------------------------------------------
class Systems(object):
    """
    System data
    """
    INFO = {\
    "mass":{"unit":"kg", "om":1.e3, "txt":"Mass of all airplane systems"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of the system CG"}
    }
    def __init__(self, mass = None,
                        c_g = None):
        self.mass = mass
        self.c_g = c_g

