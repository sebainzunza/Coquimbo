#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:55:48 2019
@author: joe
"""

# This file has been independently developed and is not part of
# the Open Drift Package.
#
# The model is intended to replicate the swimming patterns of Pelagia
# noctiluca, the most abundant jellyfish in the Western Mediterranean.
#
# Prerequisites and Dependencies:
# opendrift_p2 // Freely available on https://github.com/OpenDrift/opendrift
# pyephem project // Python package for performing high-precision astronomy
# computations. https://pypi.org/project/pyephem
# Model is python 3 compatible, but the developers provide no warranty or
# insurance
#
# Developers: Joe El Rahi, Deltares/UNIBO - info: joe.elrahi@studio.unibo.it
#                                                 joe.elrahi@deltares.nl
#             Marc Weeber, Deltares       - info: marc.weeber@deltares.nl
#
# Additional information on the behavioural modelling and references
# can be obtained by contacting one of the developers.
#
# Copyright 2019, Joe El Rahi & Marc Weeber, Deltares, Delft, The Netherlands


#The following lines should be checked before every simulation

#lines 219 & 250: change to equivalent latitude and longitude

#Important remark: This model has been developed for a Delta time step
# of 60 min. Check time step with o.list_configspec()
#For a different time step function should be adapted



from __future__ import division
import pandas
import logging
import numpy as np
import ephem
from datetime import datetime, timedelta
from opendrift.models.opendrift3D import \
     OpenDrift3DSimulation, Lagrangian3DArray
from opendrift.elements import LagrangianArray
from opendrift.elements.passivetracer import PassiveTracer
# Defining the parameters of Pelagia noctiluca
# Lagrangian array creates an array of size n
# where n is the number of particles 

class Pelagianoctiluca(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs
    """

    variables = LagrangianArray.add_variables([

    #add the test variables here
        ('time_stage', {'dtype': np.float32,
                        'units': 'min',
                        'default':1440.}),#Planulae stage development time
        ('stage',     {'dtype': np.float32,
                       'units':  '[]',
                       'default':1.}),#By default the initial stage is planulae
        ('time_stage_tot',{'dtype': np.float32,
                           'units':'min',
                           'default':1440.}),#total time of first stage
       
        ('vzact',{'dtype':np.float32,
                  'units':'[]',
                  'default':0}),
        ('bouy',{'dtype':np.float32,
                 'units':'[]',
                 'default':0})]) 
    
    
class PelagianoctilucaDrift(OpenDrift3DSimulation):
    """Pelagia noctiluca life cycle model based on the OpenDrift framework.
    
       Developed at Deltares, the Netherlands
    """
    
    ElementType = Pelagianoctiluca

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity',
                          'surface_downward_x_stress',
                          'surface_downward_y_stress',
                          'turbulent_kinetic_energy',
                          'turbulent_generic_length_scale',
                          'upward_sea_water_velocity'
                          ]

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)
    
    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-500, 0]

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.,
                       'surface_downward_x_stress': 0,
                       'surface_downward_y_stress': 0,
                       'turbulent_kinetic_energy': 0,
                       'turbulent_generic_length_scale': 0,
                       'upward_sea_water_velocity': 0
                       }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'stranded': 'red', 'eaten': 'yellow', 'died': 'magenta'}


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(PelagianoctilucaDrift, self).__init__(*args, **kwargs)

        # By default, jellyfish arriving to shore are "stranding"
        self.set_config('general:coastline_action', 'stranding')

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Information and references are available upon request by email
        """
        
        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                weight_upper + \
                Tprofiles[lower, range(Tprofiles.shape[1])] * \
                (1-weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                weight_upper + \
                Sprofiles[lower, range(Sprofiles.shape[1])] * \
                (1-weight_upper)

        ##Script developments starts here
        #Function 1 is determining the development stage, velocity and bouy
        #of the particles
        
        stage_growth_factor=np.interp(x=self.environment.sea_water_temperature,xp=[4.5,13.5,17.0,19.0,20.0,25.0],fp=[0.0,0.4,1.0,0.8,0.75,0.0])
        time_stage=self.elements.time_stage-(1)*stage_growth_factor
        self.elements.time_stage=time_stage
               
        for n in range(0,len(self.elements.time_stage)):
            if self.elements.time_stage[n] < 0:
                if self.elements.stage[n] ==1:
                    self.elements.time_stage_tot[n]=4080
                    self.elements.time_stage[n]=self.elements.time_stage_tot[n]
                    self.elements.stage[n]=self.elements.stage[n]+1
                elif self.elements.stage[n] ==2:
                   self.elements.time_stage_tot[n]=89280
                   self.elements.time_stage[n]=self.elements.time_stage_tot[n]
                   self.elements.stage[n]=self.elements.stage[n]+1
                elif self.elements.stage[n] ==3:
                   self.elements.time_stage_tot[n]=535680
                   self.elements.time_stage[n]=self.elements.time_stage_tot[n]
                   self.elements.stage[n]=self.elements.stage[n]+1
            #here you can add another elif for mortality which is not included
            #because this model is not expected to run for more than a year
        stage_development=1-(self.elements.time_stage / self.elements.time_stage_tot)
        
        for n in range(0,len(self.elements.stage)):
            if self.elements.stage[n]==1:
                self.elements.vzact[n]=0
                self.elements.bouy[n]=10
            if self.elements.stage[n]==2:
                self.elements.vzact[n]=np.interp(x=stage_development[n], xp=[0,1], fp=[0, (120./(24*60*60))])
                self.elements.bouy[n]=np.interp(x=stage_development[n], xp=[0,1],fp=[(24./(24*60*60)),(-24./(24*60*60))])
            elif self.elements.stage[n]==3:
                self.elements.vzact[n]=np.interp(x=stage_development[n], xp=[0,1], fp=[(120./(24*60*60)),(1800./(24*60*60))])
                self.elements.bouy[n]=np.interp(x=stage_development[n], xp=[0,1],fp=[(-24./(24*60*60)),(-80./(24*60*60))])
            elif self.elements.stage[n]==4:
                self.elements.vzact[n]=np.interp(x=stage_development[n], xp=[0,1], fp=[(1800./(24*60*60)),(2400./(24*60*60))])
                self.elements.bouy[n]=np.interp(x=stage_development[n], xp=[0,1],fp=[(-80./(24*60*60)),(-259./(24*60*60))])
        
        
        #Adjust the calculated swimming speed to the temperature
        motility_factor=np.interp(x=self.environment.sea_water_temperature,xp=[4.5,13.5,17,19,20,25],fp=[0,0.8,1,1,0.8,0])
        self.elements.vzact=self.elements.vzact*motility_factor
        #Adjusting the script for Diurnal Behaviour using Ephem
        test=self.get_time_array()[0]   
        momentime=test[-1] 
        spain=ephem.Observer()
        spain.lat,spain.lon = '39.02', '1.48'
        spain.date=momentime
        sr=(spain.previous_rising(ephem.Sun()))
        ss=(spain.next_setting(ephem.Sun()))
        srs=sr.datetime()
        sss=ss.datetime()
        a=momentime.hour
        b=srs.hour
        c=sss.hour
        
        if a > b and a < c:
            vz=self.elements.bouy-self.elements.vzact
        else:
            vz=self.elements.bouy+self.elements.vzact
        
            
        
        #W1 is the final speed used for the vertical behavior
        W1=vz
        #setting the maximum dive limit of the particle
        max_dive=np.interp(x=self.environment.sea_floor_depth_below_sea_level,xp=[0,50,100,300,900,1000,1600,5000],fp=[0,12,25,100,300,320,330,350])
        max_dive=-max_dive
        for n in range(0,len(self.elements.z)):
            if self.elements.z[n]<max_dive[n]:
                self.elements.z[n]=max_dive[n]
            if self.elements.z[n]>=0:
               self.elements.z[n]=0
                
        
        ##Script development ends here
        
        self.elements.terminal_velocity = W1
    def update(self):
        """Update positions and properties of buoyant particles."""

        # Turbulent Mixing
        self.update_terminal_velocity()
        
        self.vertical_mixing()
        # Horizontal advection
        self.advect_ocean_current()

        # Vertical advection
        if self.get_config('processes:verticaladvection') is True:
            self.vertical_advection()
