#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.animation as animation
import os, sys
import time

# draw the figure so the animations will work
fig = plt.gcf()
fig.show()
fig.canvas.draw()


def plt_dynamic(X_dat,Y_dat,Z_dat):
   # Convert from pandas dataframes to numpy arrays
   X, Y, Z, = np.array([]), np.array([]), np.array([])
   for i in range(len(X_dat)):
        X = np.append(X,X_dat[i])
        Y = np.append(Y,Y_dat[i])
        Z = np.append(Z,Z_dat[i])

   # create x-y points to be used in heatmap
   xi = np.linspace(X.min(),X.max(),1000)
   yi = np.linspace(Y.min(),Y.max(),1000)

   # Z is a matrix of x-y values
   zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')

   # I control the range of my colorbar by removing data 
   # outside of my range of interest
   zmin = 293
   zmax = 380
   zi[(zi<zmin) | (zi>zmax)] = None

   # Create the contour plot
   CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
                  vmax=zmax, vmin=zmin)
   plt.colorbar() 
   plt.draw()
      


# Open a file
path = "/home/versaci/project/solver/plot/data/"
dirs = os.listdir( path )
dirs.sort()
# This would print all the files and directories
for file in dirs:
   print file
   # Load data from CSV
   dat = np.genfromtxt(path+file, delimiter='\t',skip_header=0)
   X_dat = dat[:,0]
   Y_dat = dat[:,1]
   Z_dat = dat[:,2]
   plt.show()
   plt_dynamic(X_dat,Y_dat,Z_dat) 
   plt.pause(0.05)
   fig.canvas.draw()
