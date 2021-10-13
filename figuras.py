#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import sys
sys.path.append('/home/loco/opendrift/')
import opendrift
import os
from datetime import datetime, timedelta
import numpy as np
o = opendrift.open('seed1000cada1hr_15dias.nc')
o.plot(filename='figura1.png',fast=true)
