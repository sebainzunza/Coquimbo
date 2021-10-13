#
# > conda activate opendrift
# > python3
#
import opendrift
o = opendrift.open('seed1000cada1hr_15dias.nc')
o.plot(filename='figura1.png',fast=True)
