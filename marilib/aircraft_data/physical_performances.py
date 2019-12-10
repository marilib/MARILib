#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 23:22:21 2019

@author: DRUOT Thierry : original Scilab implementation
         PETEILH Nicolas : portage to Python
"""

#--------------------------------------------------------------------------------------------------------------------------------
class Aerodynamics(object):
    """
    Aerodynamic data
    """
    INFO = {\
    "cruise_lod_max":{"unit":"no_dim", "om":1.e1, "txt":"Maximum lift over drag ratio in cruise condition"},
    "cz_cruise_lod_max":{"unit":"no_dim", "om":1.e0, "txt":"Lift coefficient corresponding to maximum lift over drag"},
    "hld_conf_clean":{"unit":"no_dim", "om":1.e0, "txt":"High lift device setting in clean configuration (0 by definition)"},
    "cz_max_clean":{"unit":"no_dim", "om":1.e0, "txt":"Maximum lift coefficient in clean wing configuration"},
    "hld_conf_to":{"unit":"no_dim", "om":1.e0, "txt":"High lift device setting in take off configuration (0 < hld_conf < 0,6)"},
    "cz_max_to":{"unit":"no_dim", "om":1.e0, "txt":"Maximum lift coefficient in take off configuration"},
    "hld_conf_ld":{"unit":"no_dim", "om":1.e0, "txt":"High lift device setting in landing configuration (nominal value is 1)"},
    "cz_max_ld":{"unit":"no_dim", "om":1.e0, "txt":"Maximum lift coefficient in landing configuration"}
    }
    def __init__(self, cruise_lod_max = None,
                       cz_cruise_lod_max = None,
                       hld_conf_clean = None,
                       cz_max_clean = None,
                       hld_conf_to = None,
                       cz_max_to = None,
                       hld_conf_ld = None,
                       cz_max_ld = None):
        self.cruise_lod_max = cruise_lod_max
        self.cz_cruise_lod_max = cz_cruise_lod_max
        self.hld_conf_clean = hld_conf_clean
        self.cz_max_clean = cz_max_clean
        self.hld_conf_to = hld_conf_to
        self.cz_max_to = cz_max_to
        self.hld_conf_ld = hld_conf_ld
        self.cz_max_ld = cz_max_ld

#--------------------------------------------------------------------------------------------------------------------------------
class Propulsion(object):
    """
    Propulsion data
    """
    INFO = {\
    "architecture":{"unit":"string", "om":1.e0, "txt":"Propulsion architecture, TF:turbofan, PTE1:partial turbo electric n°1, EF1, electric with batteries"},
    "n_engine":{"unit":"int", "om":1.e0, "txt":"Number of main engines"},
    "fuel_type":{"unit":"string", "om":1.e0, "txt":"Type of fuel, Kerosene, Hydrogene, Methane, Battery"},
    "battery_energy_density":{"unit":"kWh/kg", "om":1.e0, "txt":"Battery energy density"},
    "reference_thrust":{"unit":"daN", "om":1.e5, "txt":"Reference_thrust that drives the engine size and assocoated systems"},
    "reference_thrust_effective":{"unit":"daN", "om":1.e5, "txt":"Effective reference_thrust computed as max thrust(Mach = 0.25, ISA+15, Sea Level) / 0.8"},
    "y_ext_nacelle":{"unit":"m", "om":1.e1, "txt":"Span wise position of the out most right engine nacelle"},
    "sfc_cruise_ref":{"unit":"kg/daN/h", "om":1.e0, "txt":"Specific Fuel Consumption in cruise condition, isa, ref_cruise_altp, cruise_mach"},
    "sec_cruise_ref":{"unit":"kW/daN/h", "om":1.e0, "txt":"Specific Energy Consumption of the electric chain (if any) in cruise condition, isa, ref_cruise_altp, cruise_mach"},
    "bli_effect":{"unit":"int", "om":1.e0, "txt":"BLI effect switch, 0: without, 1: with"},
    "bli_r_thrust_factor":{"unit":"no_dim", "om":1.e0, "txt":"Thrust factor at constant power due to boundary layer ingestion of the e-fan in cruise condition"},
    "bli_thrust_factor":{"unit":"no_dim", "om":1.e0, "txt":"Thrust factor at constant power due to boundary layer ingestion of other fans in cruise condition"},
    "rating_code":{"unit":"int", "om":1.e0, "txt":"Array of rating codes [0:MTO, 1:MCN, 2:MCL, 3:MCR, 4:FID]"},
    "flight_data":{"unit":"dict", "txt":"Dictionary of flying conditions for each rating {'disa':array, 'altp':array, 'mach':array, 'nei':array}"},
    "mto_thrust_ref":{"unit":"daN", "om":1.e4, "txt":"Turbofan thrust in take off rating (one engine), Sea Level, ISA+15, Mach 0,25"},
    "mcn_thrust_ref":{"unit":"daN", "om":1.e4, "txt":"Turbofan thrust in maxi continuous rating (one engine), Required ceiling altitude, ISA, cruise Mach"},
    "mcl_thrust_ref":{"unit":"daN", "om":1.e4, "txt":"Turbofan thrust in max climb rating (one engine), Required Top of Climb altitude, ISA, cruise Mach"},
    "mcr_thrust_ref":{"unit":"daN", "om":1.e4, "txt":"Turbofan thrust in max cruise rating (one engine), Reference cruise altitude, ISA, cruise Mach"},
    "fid_thrust_ref":{"unit":"daN", "om":1.e4, "txt":"Turbofan thrust in flight idle rating (one engine), Reference cruise altitude, ISA, cruise Mach"},
    "mass":{"unit":"kg", "om":1.e3, "txt":"Total mass of the propulsion system (pylons, nacelles, engines, ...)"},
    "c_g":{"unit":"m", "om":1.e1, "txt":"Global CG position for the whole propulsion system (pylons, nacelles, engines, ...)"}
    }
    def __init__(self, architecture = None,
                       n_engine = None,
                       fuel_type = None,
                       battery_energy_density = None,
                       reference_thrust = None,
                       reference_thrust_effective = None,
                       y_ext_nacelle = None,
                       sfc_cruise_ref = None,
                       sec_cruise_ref = None,
                       bli_effect = None,
                       bli_r_thrust_factor = None,
                       bli_thrust_factor = None,
                       rating_code = None,
                       flight_data = None,
                       mto_thrust_ref = None,
                       mcn_thrust_ref = None,
                       mcl_thrust_ref = None,
                       mcr_thrust_ref = None,
                       fid_thrust_ref = None,
                       mass = None,
                       c_g = None):
        self.architecture = architecture
        self.n_engine = n_engine
        self.fuel_type = fuel_type
        self.battery_energy_density = battery_energy_density
        self.reference_thrust = reference_thrust
        self.reference_thrust_effective = reference_thrust_effective
        self.y_ext_nacelle = y_ext_nacelle
        self.sfc_cruise_ref = sfc_cruise_ref
        self.sec_cruise_ref = sec_cruise_ref
        self.bli_effect = bli_effect
        self.bli_r_thrust_factor = bli_r_thrust_factor
        self.bli_thrust_factor = bli_thrust_factor
        self.rating_code = rating_code
        self.flight_data = flight_data
        self.mto_thrust_ref = mto_thrust_ref
        self.mcn_thrust_ref = mcn_thrust_ref
        self.mcl_thrust_ref = mcl_thrust_ref
        self.mcr_thrust_ref = mcr_thrust_ref
        self.fid_thrust_ref = fid_thrust_ref
        self.mass = mass
        self.c_g = c_g

#--------------------------------------------------------------------------------------------------------------------------------
class CharacteristicWeight(object):
    """
    Aircraft characteristic weights
    """
    INFO = {\
    "mwe":{"unit":"kg", "om":1.e4, "txt":"Manufacturer Weight Empty"},
    "owe":{"unit":"kg", "om":1.e4, "txt":"Operating Weight Empty (= mwe + m_op_item + m_cont_pallet)"},
    "mzfw":{"unit":"kg", "om":1.e4, "txt":"Maximum Zero Fuel Weight (= owe + n_pax_ref.m_pax_max)"},
    "mass_constraint_1":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on MZFW, must be kept positive"},
    "mlw":{"unit":"kg", "om":1.e4, "txt":"Maximum Landing Weight (close or equal to 1,07mzfw except for small aircraft where mlw = mtow)"},
    "mass_constraint_2":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on MLW, must be kept positive"},
    "mtow":{"unit":"kg", "om":1.e4, "txt":"Maximum Take Off Weight"},
    "mass_constraint_3":{"unit":"no_dim", "om":1.e0, "txt":"Constraint on MTOW, must be kept positive"},
    "mfw":{"unit":"kg", "om":1.e4, "txt":"Maximum Fuel Weight"},
    "battery_in_owe":{"unit":"kg", "om":1.e4, "txt":"Battery mass in OWE"}
    }
    def __init__(self, mwe = None,
                       owe = None,
                       mzfw = None,
                       mass_constraint_1 = None,
                       mlw = None,
                       mass_constraint_2 = None,
                       mtow = None,
                       mass_constraint_3 = None,
                       mfw = None,
                       battery_in_owe = None):
        self.mwe = mwe
        self.owe = owe
        self.mzfw = mzfw
        self.mass_constraint_1 = mass_constraint_1
        self.mlw = mlw
        self.mass_constraint_2 = mass_constraint_2
        self.mtow = mtow
        self.mass_constraint_3 = mass_constraint_3
        self.mfw = mfw
        self.battery_in_owe = battery_in_owe

#--------------------------------------------------------------------------------------------------------------------------------
class CenterOfGravity(object):
    """
    Operational positions of the center of gravity
    """
    INFO = {\
    "cg_range_optimization":{"unit":"int", "om":1.e0, "txt":"Wing position, HTP area and VTP area optimized according to HQ criteria, 0: no, 1:yes"},
    "mwe":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of MWE CG"},
    "owe":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of OWE CG"},
    "battery":{"unit":"m", "om":1.e1, "txt":"Longitudinal position of battery CG"},
    "max_fwd_mass":{"unit":"kg", "om":1.e4, "txt":"Aircraft mass at maximum forward CG"},
    "max_fwd_req_cg":{"unit":"m", "om":1.e1, "txt":"Required maximum forward aircraft CG"},
    "max_fwd_trim_cg":{"unit":"m", "om":1.e1, "txt":"Maximum trim-able forward CG"},
    "cg_constraint_1":{"unit":"no_dim", "om":1.e0, "txt":"Forward CG constraint n°1, must be kept positive"},
    "max_bwd_mass":{"unit":"kg", "om":1.e4, "txt":"Aircraft mass at maximum backward payload CG"},
    "max_bwd_req_cg":{"unit":"m", "om":1.e1, "txt":"Required maximum backward aircraft CG"},
    "max_bwd_stab_cg":{"unit":"m", "om":1.e1, "txt":"Maximum backward CG"},
    "cg_constraint_2":{"unit":"no_dim", "om":1.e0, "txt":"Backward CG constraint n°1, must be kept positive"},
    "max_bwd_oei_mass":{"unit":"kg", "om":1.e4, "txt":"Aircraft mass for OEI control criterion"},
    "max_bwd_oei_req_cg":{"unit":"m", "om":1.e1, "txt":"Required backward CG at max_bwd_oei_mass"},
    "max_bwd_oei_cg":{"unit":"m", "om":1.e1, "txt":"Maximum backward CG according to OEI control"},
    "cg_constraint_3":{"unit":"no_dim", "om":1.e0, "txt":"Backward CG constraint n°2, must be kept positive"}
    }
    def __init__(self, cg_range_optimization = None,
                       mwe = None,
                       owe = None,
                       battery = None,
                       max_fwd_mass = None,
                       max_fwd_req_cg = None,
					   max_fwd_trim_cg = None,
                       cg_constraint_1 = None,
                       max_bwd_mass = None,
                       max_bwd_req_cg = None,
                       max_bwd_stab_cg = None,
                       cg_constraint_2 = None,
                       max_bwd_oei_mass = None,
                       max_bwd_oei_req_cg = None,
                       max_bwd_oei_cg = None,
                       cg_constraint_3 = None):
        self.cg_range_optimization = cg_range_optimization
        self.mwe = mwe
        self.owe = owe
        self.battery = battery
        self.max_fwd_mass = max_fwd_mass
        self.max_fwd_req_cg = max_fwd_req_cg
        self.max_fwd_trim_cg = max_fwd_trim_cg
        self.cg_constraint_1 = cg_constraint_1
        self.max_bwd_mass = max_bwd_mass
        self.max_bwd_req_cg = max_bwd_req_cg
        self.max_bwd_stab_cg = max_bwd_stab_cg
        self.cg_constraint_2 = cg_constraint_2
        self.max_bwd_oei_mass = max_bwd_oei_mass
        self.max_bwd_oei_req_cg = max_bwd_oei_req_cg
        self.max_bwd_oei_cg = max_bwd_oei_cg
        self.cg_constraint_3 = cg_constraint_3

