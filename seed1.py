#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import sys
sys.path.append('/home/loco/opendrift/')

import os
from datetime import datetime, timedelta
import numpy as np
#import pandas as pd
#colnames=['LON', 'LAT']
#lons= np.loadtxt("PuntosCosta.txt")[:, 0]
#lats= np.loadtxt("PuntosCosta.txt")[:, 1]
#data = pd.read_csv('PuntosCosta.txt')
#print(lons)
#exit()
#from opendrift.readers import reader_global_landmask
#from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)
#reader_landmask = reader_global_landmask.Reader('cobertura.txt')
filename = 'coquimbo_avg.nc';
#o.set_config('general:basemap_resolution','i')
#gridfile = '/home/seba/opendrift/Tutorial_opendrift/croco_grd.nc';
mosa_native = reader_ROMS_native.Reader(filename)
mosa_native.interpolation='linearND'
o.add_reader(mosa_native)
print(mosa_native)
#exit()
# Agregando elementos
time=mosa_native.start_time #tiempo de seeding
outname='seed1';
#z = np.linspace(0, -50, 100)  # Linearly increasing depth
#lons = np.linspace(-71.55, -71.33, 100) # lon de sedding elementos entre dos puntos
#lats = np.linspace(-29.7, -29.88, 100) # lat seeding elementos entre dos puntos
# o.seed_cone(lon=[-74.5, -74.1], lat=[-42.1, -41.6], number=1000, radius=[0, 5000], time=time)
# o.seed_element(lon=[-71.55, -71.33], lat=[-29.88, -27.7], number=1000, radius=[0, 5000], time=time)
o.seed_elements(lon=-71.33, lat=-29.87, time=time)
#o.set_config('environment:fallback:y_sea_water_velocity', 3)  # Adding some current to be able to visualise depth as color of trajectories
o.run(outfile = outname + '.nc')
print(o)
o.plot(linecolor='z', filename = outname + '.png')
#o.set_config('environment:fallback:y_sea_water_velocity', 0)
o.animation(filename = outname + '.mp4')
