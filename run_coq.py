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
#from opendrift.models.larvaloco import SeaLice
from opendrift.models.oceandrift import OceanDrift
o = OceanDrift(loglevel=0)
#reader_landmask = reader_global_landmask.Reader('cobertura.txt')
filename = 'coquimbo_avg.nc';
filecoord= 'Hab_Rocoso_Coquimbo.csv'
#outname='figura1_runcoq';
#o.set_config('general:basemap_resolution','i')
#gridfile = '/home/seba/opendrift/Tutorial_opendrift/croco_grd.nc';
mosa_native = reader_ROMS_native.Reader(filename)
mosa_native.interpolation='linearND'
o.add_reader(mosa_native)
print(mosa_native)
#exit()
# Agregando elementos
start_time = mosa_native.start_time #tiempo de seeding
pasos_total = 24*30 # pasos total de simulacion 
total_days_seed = 15
num_particles = 2000   # 7200 =1 per hour  for 1 month
lons = []
lats = []
with open(filecoord, 'r') as f:
    lines = f.readlines()
    for line in lines:
        line_split = line.split()
        lons.append(float(line_split[0]))
        lats.append(float(line_split[1]))
    f.close()
aux_latlon=np.arange(0,num_particles)*0
#radiuses = aux_latlon + self.confobj.release_radius
all = np.array([0, 0, 0, 0, 0, 0, 0])
print(len(lats))
paso = 1
for lats, lons in zip(lats, lons):
#    avance = "Punto_{}_de_{}".format(str(paso),str(len(lats)))
#    print(avance)
    o.seed_elements(lon=lons,
       lat=lats,
       number=num_particles,
       time=[start_time, start_time + timedelta(days=total_days_seed)],
       z=-20)

    #self.logger.info('Elements scheduled for {} : {}'.format(self.confobj.species, o.elements_scheduled))

    lons_start = o.elements_scheduled.lon
    lats_start = o.elements_scheduled.lat
    time_start = o.elements_scheduled_time

    o.run(steps=pasos_total,
        time_step=timedelta(minutes=60),
        export_variables=['status','time'])

    index_of_first, index_of_last = o.index_of_activation_and_deactivation()
    lons = o.get_property('lon')[0]
    lats = o.get_property('lat')[0]
    status = o.get_property('status')[0]
    tiempo = o.get_time_array()[0]

    lons_end = lons[index_of_last, range(lons.shape[1])]
    lats_end = lats[index_of_last, range(lons.shape[1])]
    status_end = status[index_of_last, range(lons.shape[1])]
    time_end   = [tiempo[i] for i in index_of_last]

    uno = np.column_stack([time_start, lons_start, lats_start, time_end, lons_end, lats_end, status_end])

    all = np.row_stack([all, uno])
    paso = paso + 1

np.savetxt('Loco_Coq_Uniforme_IF_30dias_1000part.txt', all, fmt='%s %.4f %.4f %s %.4f %.4f %i')
#o.plot(linecolor='z', filename = outname + '.png')
#o.set_config('environment:fallback:y_sea_water_velocity', 0)
#o.animation(filename = outname + '.mp4')
print("Fin simulacion")

