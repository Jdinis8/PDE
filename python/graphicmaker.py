import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d

my_dpi = 200

filename = os.getcwd() + '/graphics/output.txt'

ctr = 0

#important variables
size_t = 0
size_x = 0
space_step = 0
time_step = 0
x0 = 0
evolution = []

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
      dumb = []
      for x in line.split():
        dumb.append(float(x))
      evolution.append(dumb)

#now we have everything from the output.txt stored

#here we shall create the array for the discretization in x
x = []

for i in range(size_x):
  x.append(x0 + i*space_step)

maxvalue = evolution[0][0]
minvalue = evolution[0][0]

for i in range(len(evolution)):
  for j in range(len(evolution[i])):
    if(evolution[i][j] > maxvalue):
      maxvalue = evolution[i][j]
    if(evolution[i][j] < minvalue):
      minvalue = evolution[i][j]

fctr = 0
for i in range(size_t):
  if fctr%1000 == 0:
    plt.scatter(x, evolution[i], marker = "o", s=10)
    plt.xlabel('Space (x)')
    plt.ylabel('Value (f(x,t))')
    plt.xlim([x0, x0 + size_x*space_step])
    plt.ylim([minvalue, maxvalue])
    plt.title(f"Time (t=" + "{0:.5f}".format(i*time_step, 5) + "s)")
    plt.savefig(os.getcwd() + '/graphics/time_' + "{0:.5f}".format(i*time_step, 5) + '.png', dpi=my_dpi)
    plt.close()
  fctr = fctr +1
