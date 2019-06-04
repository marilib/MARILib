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
        """
        Constructor :
            :param n_pax_ref: Uu int - OoM 10^2 - Reference Number of passengers (often 2 class layout)
            :param n_aisle: Uu int - OoM 10^0 - Number of aisle in economic section
            :param n_pax_front: Uu int - OoM 10^0 - Number of seats in a row in economic section
            :param fwd_limit: Uu m - OoM 10^0 - Distance between aircraft nose and cabin forward limit
            :param width: Uu m - OoM 10^0 - Maximum width of the cabin (volume width, not floor width)
            :param length: Uu m - OoM 10^1 - Total length of the cabin
            :param projected_area: Uu m2 - OoM 10^2 - Area of the cabin taking into account its maximum width (not floor area)
            :param m_furnishing: Uu kg - OoM 10^3 - Total mass of furnishing equipements
            :param m_op_item: Uu kg - OoM 10^3 - Total mass of operator items
            :param cg_furnishing: Uu m - OoM 10^1 - Center of gravity of furnishing equipements
            :param cg_op_item: Uu m - OoM 10^1 - Center of gravity of operator items
        """
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
        """
            :param m_pax_nominal: Uu kg - OoM 10^2 - Mass allowance per passenger to compute nominal payload
            :param m_pax_max: Uu kg - OoM 10^2 - Mass allowance per passenger to compute maximum payload
            :param m_container_pallet: Uu kg - OoM 10^3 - Mass of containers or pallets empty
            :param nominal: Uu kg - OoM 10^4 - Mass of nominal payload
            :param maximum: Uu kg - OoM 10^4 - Mass of maximum payload
            :param max_fwd_mass: Uu kg - OoM 10^2 - Payload mass at maximum forward payload CG
            :param max_fwd_req_cg: Uu m - OoM 10^1 - Required maximum forward payload CG
            :param max_bwd_mass: Uu kg - OoM 10^2 - Payload mass at maximum backward payload CG
            :param max_bwd_req_cg: Uu m - OoM 10^1 - Required maximum backward payload CG
            :param cg_container_pallet: Uu m - OoM 10^1 - Center of gravity of containers or pallets empty
        """
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
    def __init__(self, width = None,
                       height = None,
                       length = None,
                       tail_cone_length = None,
                       net_wetted_area = None,
                       mass = None,
                       c_g = None):
        """
        Constructor :
            :param width: Uu m - OoM 10^0 - Fuselage width of the cylindrical part
            :param height: Uu m - OoM 10^0 - Fuselage height of the cylindrical part
            :param length: Uu m - OoM 10^1 - Total fuselage length
            :param tail_cone_length: Uu m - OoM 10^0 - Length of rear evolutive part of the fuselage
            :param net_wetted_area: Uu m2 - OoM 10^2 - Fuselage total net wetted area
            :param mass: Uu kg - OoM 10^3 - Equipped fuselage mass (without systems)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the fuselage CG 
        """
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
        """
        Constructor :
            :param attachment: Uu int - OoM 10^0 - Wing attachment, 1: low wing, 2: high wing
            :param morphing: Uu int - OoM 10^0 - Wing deformation driver, 1: aspect ratio, 2: span
			:param hld_type: Uu int - OoM 10^0 - Type of high lift devices
            :param t_o_c_r: Uu % - OoM 10^1 - Thickness over chord ratio of the wing at root
            :param t_o_c_k: Uu % - OoM 10^1 - Thickness over chord ratio of the wing at main kink
            :param t_o_c_t: Uu % - OoM 10^1 - Thickness over chord ratio at wing tip
            :param sweep: Uu deg - OoM 10^0 - Wing sweep angle at 25% of the chord
            :param dihedral: Uu deg - OoM 10^0 - Mean dihedral of the wing
            :param setting: Uu deg - OoM 10^0 - Setting angle of the wing at root
            :param taper_ratio: Uu no_dim - OoM 10^0 - Wing taper ratio
            :param aspect_ratio: Uu no_dim - OoM 10^0 - Wing aspect ratio
            :param area: Uu m2 - OoM 10^2 - Wing reference area
            :param span: Uu m - OoM 10^1 - Wing span
            :param mac: Uu m - OoM 10^0 - Mean aerodynamic chord of the wing
            :param net_wetted_area: Uu m2 - OoM 10^2 - Wing total net wetted area
            :param mass: Uu kg - OoM 10^3 - Equipped wing mass (without systems)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the wing CG
            :param x_root: Uu m - OoM 10^1 - Longitudinal position of 0% of wing root chord
            :param y_root: Uu m - OoM 10^0 - Span wise position of 0% of the wing root chord
            :param z_root: Uu m - OoM 10^0 - Vertical position of 0% of the wing root chord
            :param c_root: Uu m - OoM 10^0 - Wing root chord length
            :param x_kink: Uu m - OoM 10^1 - Longitudinal position of 0% of wing kink chord
            :param y_kink: Uu m - OoM 10^0 - Span wise position of 0% of the wing kink chord
            :param z_kink: Uu m - OoM 10^0 - Vertical position of 0% of the wing kink chord
            :param c_kink: Uu m - OoM 10^0 - Wing kink chord length
            :param x_tip: Uu m - OoM 10^1 - Longitudinal position of 0% of wing tip chord
            :param y_tip: Uu m - OoM 10^0 - Span wise position of 0% of the wing tip chord
            :param z_tip: Uu m - OoM 10^0 - Vertical position of 0% of the wing tip chord
            :param c_tip: Uu m - OoM 10^0 - Wing tip chord length
            :param x_mac: Uu m - OoM 10^1 - Longitudinal position of wing mean aerodynamic chord
            :param y_mac: Uu m - OoM 10^0 - Span wise position of wing mean aerodynamic chord
        """
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
    def __init__(self, mass = None,
                        c_g = None):
        """
        Constructor :
            :param mass: Uu kg - OoM 10^3 - Mass of landing gears (nose and main)
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the landing gears CG
        """
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class HorizontalTail(object):
    """
    Horizontal tail plane data
    """
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
        """
        Constructor :
            :param attachment: Uu int - OoM 10^0 - Configuration of horizontal tail, 1: classical, 2: T-tail
            :param sweep: Uu deg - OoM 10^0 - Horizontal tail sweep angle at 25% of the chords
            :param taper_ratio: Uu no_dim - OoM 10^0 - Taper ratio of the horizontal tail
            :param aspect_ratio: Uu no_dim - OoM 10^0 - Aspect ratio of the horizontal tail
            :param t_o_c: Uu no_dim - OoM 10^0 - Thickness to chord ratio of the vertical tail
            :param dihedral: Uu deg - OoM 10^0 - Mean dihedral of the horizontal tail
            :param volume: Uu no_dim - OoM 10^0 - Volume coefficient of the  horizontal tail
            :param lever_arm: Uu m - OoM 10^1 - Lever arm of the horizontal tail (from 25% wing MAC to 25% HTTP MAC
            :param area: Uu m2 - OoM 10^2 - Horizontal tail reference area
            :param span: Uu m - OoM 10^1 - Horizontal tail span
            :param mac: Uu m - OoM 10^0 - Mean aerodynamic part of the horizontal tail
            :param net_wetted_area: Uu m2 - OoM 10^2 - Total net wetted area of the horizontal tail
            :param mass: Uu kg - OoM 10^2 - Equipped mass of the horizontal tail
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the horizontal tail
            :param x_axe: Uu m - OoM 10^1 - Longitudinal position of the horizontal tail central chord
            :param z_axe: Uu m - OoM 10^1 - Vertical position of the horizontal tail central chord
            :param c_axe: Uu m - OoM 10^1 - Horizontal tail central chord
            :param x_tip: Uu m - OoM 10^1 - Longitudinal position of the horizontal tail tip chord
            :param y_tip: Uu m - OoM 10^1 - Lateral position of the horizontal tail tip chord
            :param z_tip: Uu m - OoM 10^1 - Vertical position of the horizontal tail tip chord
            :param c_tip: Uu m - OoM 10^1 - Horizontal tail tip chord
            :param x_mac: Uu m - OoM 10^1 - Longitudinal position of the horizontal tail mean aerodynamic chord
            :param y_mac: Uu m - OoM 10^1 - Lateral position of the horizontal tail mean chord
        """
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
    def __init__(self, sweep = None,
                       taper_ratio = None,
                       aspect_ratio = None,
                       t_o_c = None,
                       dihedral = None,
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
        """
        Constructor :
            :param sweep: Uu deg - OoM 10^0 - Vertical tail sweep angle at 25% of the chords
            :param taper_ratio: Uu no_dim - OoM 10^0 - Taper ratio of the vertical tail
            :param aspect_ratio: Uu no_dim - OoM 10^0 - Aspect ratio of the vertical tail
            :param t_o_c: Uu no_dim - OoM 10^0 - Thickness to chord ratio of the vertical tail
            :param dihedral: Uu deg - OoM 10^0 - Mean dihedral of the vertical tail
            :param volume: Uu m2/kN - OoM 10^0 - Volume coefficient of the  vertical tail
            :param lever_arm: Uu m - OoM 10^1 - Lever arm of the vertical tail (from 25% wing MAC to 25% HTTP MAC
            :param area: Uu m2 - OoM 10^2 - Vertical tail reference area
            :param height
            :param mac: Uu m - OoM 10^0 - Mean aerodynamic part of the vertical tail
            :param net_wetted_area: Uu m2 - OoM 10^2 - Total net wetted area of the vertical tail
            :param mass: Uu kg - OoM 10^2 - Equipped mass of the vertical tail
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the CG of the vertical tail
            :param x_root: Uu m - OoM 10^1 -
            :param x_root: Uu m - OoM 10^2 - Longitudinal position of the vertical tail root chord
            :param z_root: Uu m - OoM 10^1 - Vertical position of the vertical tail root chord
            :param c_root: Uu m - OoM 10^1 - Vertical tail root chord
            :param x_tip : Uu m - OoM 10^2 - Longitudinal position of the vertical tail tip chord
            :param z_tip : Uu m - OoM 10^1 - Vertical position of the vertical tail tip chord
            :param c_tip : Uu m - OoM 10^1 - Vertical tail tip chord
            :param x_mac : Uu m - OoM 10^2 - Longitudinal position of the vertical tail mean aerodynamic chord
        """
        self.sweep = sweep
        self.taper_ratio = taper_ratio
        self.aspect_ratio = aspect_ratio
        self.t_o_c = t_o_c
        self.dihedral = dihedral
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
    def __init__(self, cantilever_volume = None,
                       central_volume = None,
                       body_volume = None,
                       mfw_volume_limited = None,
                       fuel_density = None,
                       fuel_cantilever_cg = None,
                       fuel_central_cg = None,
                       fuel_body_cg = None,
                       fuel_max_fwd_mass = None,
                       fuel_max_fwd_cg = None,
                       fuel_max_bwd_mass = None,
                       fuel_max_bwd_cg = None,
                 ):
        """
        Constructor :
            :param cantilever_volume: Uu m3 - OoM 10^1 - Volume of tanks in the cantilever wing
            :param central_volume: Uu m3 - OoM 10^1 - Volume of tanks in the central part of the wing (inside the fuselage)
            :param body_volume: Uu m3 - OoM 10^1 - Volume of tanks in the external bodies
            :param mfw_volume_limited: Uu kg - OoM 10^4 - Maximum geometrical fuel volume
            :param fuel_density: Uu kg/m3 - OoM 10^2 - Fuel density
            :param fuel_cantilever_cg: Uu m - OoM 10^1 - Center of gravity of tanks in the cantilever wing
            :param fuel_central_cg: Uu m - OoM 10^1 - Center of gravity of tanks in the central part of the wing (inside the fuselage)
            :param fuel_body_cg: Uu m - OoM 10^1 - Center of gravity of tanks in the external bodies
            :param fuel_max_fwd_mass: Uu kg - OoM 10^3 - Fuel mass of max forward fuel cg
            :param fuel_max_fwd_cg: Uu m - OoM 10^1 - Max forward fuel cg
            :param fuel_max_bwd_mass: Uu kg - OoM 10^3 - Fuel mass of max backward fuel cg
            :param fuel_max_bwd_cg: Uu m - OoM 10^1 - Max backward fuel cg
        """
        self.cantilever_volume = cantilever_volume
        self.central_volume = central_volume
        self.body_volume = body_volume
        self.mfw_volume_limited = mfw_volume_limited
        self.fuel_density = fuel_density
        self.fuel_cantilever_cg = fuel_cantilever_cg
        self.fuel_central_cg = fuel_central_cg
        self.fuel_body_cg = fuel_body_cg
        self.fuel_max_fwd_mass = fuel_max_fwd_mass
        self.fuel_max_fwd_cg = fuel_max_fwd_cg
        self.fuel_max_bwd_mass = fuel_max_bwd_mass
        self.fuel_max_bwd_cg = fuel_max_bwd_cg

#--------------------------------------------------------------------------------------------------------------------------------
class Systems(object):
    """
    System data
    """
    def __init__(self, mass = None,
                        c_g = None):
        """
        Constructor :
            :param mass: Uu kg - OoM 10^3 - Mass of all airplane systems
            :param c_g: Uu m - OoM 10^1 - Longitudinal position of the system CG
        """
        self.mass = mass
        self.c_g = c_g

