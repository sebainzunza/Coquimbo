#!/usr/bin/env python3
"""
@author: sebainzunza
"""
#
#
#
import datetime
import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
# Defining the parameters of Larva loco
# Lagrangian array creates an array of size n
# where n is the number of particles 

class LarvaLocoElement(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for larva loco
    """
    variables = Lagrangian3DArray.add_variables([
        ('growth_rate', {'dtype': np.float32,
                         'units': 'um/day',
                         'default': 0,
                         'description': 'Growth rate of virtual loco larva'}),
        ('length', {'dtype': np.float32,
                    'units': 'um',
                    'default': 250,
                    'description': 'Length of virtual loco larva'})])

class LarvaLoco(OceanDrift):
    """Larva Loco (Conchalepas conchalepas) life cycle model based on the OpenDrift framework.
    """
    
    ElementType = LarvaLocoElement
    
    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'ocean_vertical_diffusivity': {'fallback': 0, 'profiles': True},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'land_binary_mask': {'fallback': None},
        'sea_water_temperature': {'fallback': 11},
        'sea_water_salinity': {'fallback': 34}}
        
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [0, -50]

    def __init__(self, *args, **kwargs):
        # Calling general constructor of parent class
        super(LarvaLoco, self).__init__(*args, **kwargs)
        
    def larva_growth(self,T):
        """Growth rate in um/dt as a function of sea temperature 
           for each virtual loco larva (Garavelli et al. 2016)
        """
    
        beta= 4.587
        dt = self.time_step.total_seconds()
    #    T = self.environment.sea_water_temperature
        PLD = np.exp(beta -1.34*np.log(T/15) - (0.28*(np.log(T/15))**2))
        GR = ((1900 - 250)/PLD) 
        self.elements.growth_rate = GR
        return GR*(dt/86400) 
    def update_larvaloco(self):
        length_larva = self.larva_growth(self.environment.sea_water_temperature)
        self.elements.length += length_larva
        
    def update(self):
        self.update_larvaloco()
        self.update_positions(self.environment.x_sea_water_velocity, self.environment.y_sea_water_velocity)
