import sys
sys.path.append('/home/loco/opendrift/')

import os
from datetime import datetime, timedelta
import numpy as np
#import pandas as pd
#from opendrift.readers import reader_global_landmask
#from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
#from opendrift.models.larvaloco import SeaLice
from opendrift.models.oceandrift import OceanDrift
o = OceanDrift(loglevel=0)

filename = 'coquimbo_avg.nc';
filecoord= 'Hab_Rocoso_Coquimbo.csv'
mosa_native = reader_ROMS_native.Reader(filename)
mosa_native.interpolation='linearND'
o.add_reader(mosa_native)
print(mosa_native)

start_date = mosa_native.start_time #tiempo de seeding
end_date = mosa_native.end_time
dias_deriva = 15
num_particles = 1000   # 7200 =1 per hour  for 1 month
mindep = -40
maxdep = -60
st_lons = []
st_lats = []
with open(filecoord, 'r') as f:
    lines = f.readlines()
    for line in lines:
        line_split = line.split()
        st_lons.append(float(line_split[0]))
        st_lats.append(float(line_split[1]))
    f.close()

zd = np.random.rand(num_particles)
zd = zd*(maxdep - mindep) + mindep

for i in range(0,1):
    ini_date = start_date + timedelta(days=i*30)
    print(ini_date)
    all = all = np.array([0, 0, 0, 0, 0, 0, 0])
    paso = 1

    for lat, lon in zip(st_lats, st_lons):
        print(lat)
        print(lon)

        o.seed_elements(lon=lon,
           lat=lat,
           radius = 1000,
           number=num_particles,
           time=ini_date,
           z=zd)

        lons_start = o.elements_scheduled.lon
        lats_start = o.elements_scheduled.lat
        time_start = o.elements_scheduled_time

        o.run(end_time = ini_date + timedelta(days=dias_deriva),
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

np.savetxt('Loco_Coq_Uniforme_IF_'+str(dias_deriva)+'dias_'+str(num_particles)+'part.txt', all, fmt='%s %.4f %.4f %s %.4f %.4f %i')
#o.plot(linecolor='z', filename = outname + '.png')
#o.set_config('environment:fallback:y_sea_water_velocity', 0)
#o.animation(filename = outname + '.mp4')
print("Fin simulacion")

