import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d

my_dpi = 200

filename = '/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/output.txt'

ctr = 0

#important variables
size_t = 0
size_x = 0
space_step = 0
time_step = 0
evolution = [[]]

with open(filename) as f:

  for line in f:
    ##first line of .txt has width, length, nvar and tilesize
    if(ctr == 0):
      size_t = float(line.split()[0])
      size_x = float(line.split()[1])
      space_step = float(line.split()[2])
      time_step = float(line.split()[3])
    else:
      dumb = []
      for x in line.split():
        dumb.append(float(x))
      evolution.append(dumb)
    ctr += 1

#now we have everything stored