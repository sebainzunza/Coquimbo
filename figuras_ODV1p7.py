#
# > conda activate opendrift
# > python3
#
import opendrift
import numpy as np

o = opendrift.open('seed1000cada1hr_15dias.nc', elements=np.arange(0,392000,100)
o.plot(filename='figura1.png',fast=True)
