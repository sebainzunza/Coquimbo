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
outname='seed1000cada1hr_15dias';
for hour in range(15*24):
       # print("=> Releasing {} particles within a radius of {} m on {} for each lat/lon location".format(
       #     confobj.releaseParticles,
       #     confobj.releaseRadius,
       #     confobj.startdate + timedelta(days=day)))
        o.seed_elements(lon=-71.33,
                        lat=-29.87,
                        number=100,
                        radius=1000,
                        time=time + timedelta(hours=hour),
                        z=-10) # "seafloor+0.5")  # z=-30 + (1.0 * random.randint(0, 10)))
time += timedelta(hours=1)
lons_start = o.elements_scheduled.lon
lats_start = o.elements_scheduled.lat

o.run()

index_of_first, index_of_last = o.index_of_activation_and_deactivation()
lons = o.get_property('lon')[0]
lats = o.get_property('lat')[0]
status = o.get_property('status')[0]
lons_end = lons[index_of_last, range(lons.shape[1])]
lats_end = lats[index_of_last, range(lons.shape[1])]
status_end = status[index_of_last, range(lons.shape[1])]

all = np.column_stack([lons_start, lats_start, lons_end, lats_end, status_end])
np.savetxt('test_inicialfinal.txt', all, fmt='%.4f')
print("----------")
#o.run(outfile = outname + '.nc')
#o.plot(linecolor='z', filename = outname + '.png')
#o.set_config('environment:fallback:y_sea_water_velocity', 0)
#o.animation(filename = outname + '.mp4')
