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
evolution = []

with open(filename) as f:

  for line in f:
    ##first line of .txt has width, length, nvar and tilesize
    if(ctr == 0):
      size_t = int(line.split()[0])
      size_x = int(line.split()[1])
      space_step = float(line.split()[2])
      time_step = float(line.split()[3])
      ctr += 1
    else:
      dumb = []
      for x in line.split():
        dumb.append(float(x))
      evolution.append(dumb)

#now we have everything from the output.txt stored

#here we shall create the array for the discretization in x
x = []

for i in range(size_x):
  x.append(i*space_step)

for i in range(size_t):
  plt.scatter(x, evolution[i], marker = "o", s=10)
  plt.xlabel('Space (x)')
  plt.xlabel('Value (f(x,t))')
  plt.title(f"Time (t={i*time_step})")
  plt.savefig('/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/time_' + "{0:.5f}".format(i*time_step, 5) + '.png', dpi = my_dpi)
  plt.close()
