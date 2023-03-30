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
time_step = 0
xl = []
evolution = []

with open(filename) as f:

  for line in f:
    ##first line of .txt has width, length, nvar and tilesize
    if(ctr == 0):
      time_step = float(line.split()[0])
      ctr += 1
    elif(ctr ==1):
      for x in line.split():
        xl.append(float(x))
      ctr += 1
    else:
      dumb = []
      for x in line.split():
        dumb.append(float(x))
        ctr += 1
      evolution.append(dumb)

#now we have everything from the output.txt stored
fctr = 0
for i in range(len(evolution)):
  if fctr%100 == 0:
    plt.plot(xl, evolution[i], marker = "o")
    plt.xlabel('Space (x)')
    plt.ylabel('Value (u(x,t))')
    plt.xlim([min(xl), max(xl)])
    plt.ylim([-1, 1])
    plt.title(f"Time (t=" + "{0:.5f}".format(i*time_step, 5) + "s)")
    plt.savefig('/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/time_' + "{0:.5f}".format(i*time_step, 5) + '.png', dpi = my_dpi)
    plt.close()
  fctr = fctr +1
