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
tot = []

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
      conv_data = []
      for x in line.split():
        conv_data.append(float(x))
      tot.append(conv_data)

#now we have everything from the conv.txt stored

#here we shall create the array for the discretization in x
x = []

for i in range(size_x):
  x.append(x0 + i*space_step)

plt.plot(x, tot[0], marker="o", c="hotpink", linewidth = 1, markersize = 4, label = r'$f^p(u_h - u_{fh})$')
plt.plot(x, tot[1], marker="o", c="#88c999", linewidth = 1, markersize = 4, label = r'$u_{fh} - u_{f^2h}$')
plt.xlabel('Space (x)')
plt.xlim([x0, x0 + size_x*space_step])
plt.ylim([min(min(tot))-0.1*abs(min(min(tot))),
          max(max(tot))+0.1*abs(max(max(tot)))])
plt.title(f"Time (t=" + "{0:.5f}".format(size_t*time_step, 5) + "s)")
plt.legend()
plt.savefig('/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/convergencetest_time_' +
                "{0:.5f}".format(size_t*time_step, 5) + '.png', dpi=my_dpi)
plt.close()


##now for L2 norm

filename = '/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/l2.txt'
ctr = 0
x = []
y = []

with open(filename) as f:

  for line in f:
    ##first line of .txt has width, length, nvar and tilesize
    if(ctr == 0):
      for item in line.split():
        x.append(float(item))
      ctr = ctr + 1
    elif(ctr == 1):
      for item in line.split():
        y.append(float(item))
    else:
      break

plt.plot(x, y, marker="o", c="hotpink", linewidth = 1, markersize = 2)
plt.ylim([-0.5,
          max(y) + 0.1*abs(max(y))])
plt.xlabel('Time (t)')
plt.ylabel('L2 Norm')
plt.axhline(y=2, color="#88c999")
#plt.gca().invert_xaxis() #in case you want to change the order of the graphic
plt.savefig('/home/machado/Desktop/IST/4ano_2semestre/TAFC/code/graphics/l2norm.png', dpi=my_dpi)
plt.close()
