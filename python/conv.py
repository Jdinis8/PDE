import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d

my_dpi = 200

filename = '/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/conv.txt'

ctr = 0

#important variables
size_t = 0
size_x = 0
space_step = 0
time_step = 0
x0 = 0
conv_data = []

with open(filename) as f:

  for line in f:
    ##first line of .txt has width, length, nvar and tilesize
    if(ctr == 0):
      size_t = int(line.split()[0])
      size_x = int(line.split()[1])
      space_step = float(line.split()[2])
      time_step = float(line.split()[3])
      x0 = float(line.split()[4])
      ctr += 1
    else:
      for x in line.split():
        conv_data.append(float(x))

#now we have everything from the conv.txt stored

#here we shall create the array for the discretization in x
x = []

for i in range(size_x):
  x.append(x0 + i*space_step)

plt.scatter(x, conv_data, marker="o", s=10)
plt.xlabel('Space (x)')
plt.ylabel('Convergence high-medium/medium-low')
plt.xlim([x0, x0 + size_x*space_step])
plt.ylim([min(conv_data), max(conv_data)])
plt.title(f"Time (t=" + "{0:.5f}".format(size_t*time_step, 5) + "s)")
plt.savefig('/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/convergencetest_time_' +
                "{0:.5f}".format(size_t*time_step, 5) + '.png', dpi=my_dpi)
plt.close()
