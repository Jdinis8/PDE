import numpy as np
import sys
import os
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d

my_dpi = 200

files = os.listdir('/home/machado/Desktop/FireTracer/code/Generated/MapStack');

filename = '';

f_ctr = 0;
cond_ter = True;
clim_max = 1;

for file in files:

  #cleaning stuff
  ctr = 0;
  y_ctr = 0;
  x_ctr = 0;
  width = 0;
  length = 0;
  nvar = 0;
  tilesize = 0;

  xp = [];
  yp = [];
  x2d = [];
  y2d = [];
  mapinfo = [];
  terrain = [];
  terrain2d = [];
  veg = [];
  veg2d = [];
  fire = [];
  fire2d = [];

  filename = "/home/machado/Desktop/FireTracer/code/Generated/MapStack/" + file;

  with open(filename) as f:

    for line in f:

      ##first line of .txt has width, length, nvar and tilesize
      if(ctr == 0):
        for x in line.split():
          if(x.isdigit()):
            mapinfo.append(float(x));
        if(len(mapinfo) == 4):
          width = mapinfo[0];
          length = mapinfo[1];
          nvar = mapinfo[2];
          tilesize = mapinfo[3];
      
      ##then, terrain. The terrain information is contained in
      ##the first "length" lines
      if(ctr > 0 and ctr < length+1):
        dumb = [];
        y_ctr = 0;
        x2d.append((ctr-1)*tilesize);
        y2d.append((ctr-1)*tilesize);
        for x in line.split():
          dumb.append(float(x));
          terrain.append(float(x));
          xp.append(x_ctr*tilesize);
          yp.append(y_ctr*tilesize);
          y_ctr = y_ctr + 1;
        x_ctr = x_ctr + 1;
        terrain2d.append(dumb);
      
      ##vegetation. There is a two line gap between terrain info
      ##and vegetation info. Thus, the diff in the counter. It ends in (length+2) + length + 1 = 2*length + 3
      elif(ctr > length+2 and ctr < 2*length+3):
        dumb = [];
        for x in line.split():
          dumb.append(float(x));
          veg.append(float(x));
        veg2d.append(dumb);
      
      ##fire. The same two line separation happens here
      elif(ctr > 2*length+4 and ctr < 3*length+5):
        dumb = [];
        for x in line.split():
          dumb.append(float(x));
          fire.append(float(x));
        fire2d.append(dumb);
        
      ctr = ctr + 1;
    
    ##printing a contour plot for the terrain and another one for vegetation
    if(cond_ter):
      clim_max = np.max(veg2d); ##we retrieve the maximum to define the limits of our color map below
      ticks = np.linspace(0, clim_max, 10, endpoint=True) ##defining the ticks for the colorbar
      
      cond_ter = False; ##we only do this once
      plt.contourf(x2d, y2d, terrain2d, 20, cmap='viridis');
      plt.colorbar();
      plt.savefig('/home/machado/Desktop/FireTracer/code/Generated/Terrain/MapContourTerrain.png', dpi = my_dpi);
      plt.close();

      plt.contourf(x2d, y2d, veg2d, 20, cmap='viridis');
      plt.colorbar();
      plt.savefig('/home/machado/Desktop/FireTracer/code/Generated/Vegetation/MapContourVegetation.png', dpi = my_dpi);
      plt.close();

    plt.contourf(x2d, y2d, fire2d, 20, cmap='inferno');
    plt.contourf(x2d, y2d, veg2d, 20, cmap='viridis', alpha=0.5);
    plt.colorbar(ticks=ticks, label='Vegetation Density')
    plt.clim(0, clim_max);
    plt.savefig('/home/machado/Desktop/FireTracer/code/Generated/Python/MapContourFire_' + file.rsplit('.',1)[0] + '.png', dpi = my_dpi);
    plt.close();

    f_ctr = f_ctr + 1;



#if you want 3d things
#fig = plt.figure()
#ax = plt.axes(projection='3d')
##ax.plot3D(xp,yp,terrain, 'red')
#ax.scatter(xp,yp,terrain,c=terrain,cmap='viridis', linewidth = 0.01)
#ax.set_zlabel('Elevation(m)')
#ax.set_xlabel('X Pos(m)')
#ax.set_ylabel('Y Pos(m)');
#plt.savefig('3dmap.png', dpi = my_dpi)
##plt.show();
#plt.close();
