"""
=================
An animated image
=================

This example demonstrates how to animate an image.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import os.path

fig = plt.figure()

def updatefig(file):
    global xi, yi
    # Load data from CSV
    dat = np.genfromtxt(path+file, delimiter=' ',skip_header=0)
    X_dat = dat[:,0]
    Y_dat = dat[:,1]
    Z_dat = dat[:,2]

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
    zi = griddata((X, Y),Z,(xi[None,:], yi[:,None]), method='cubic')

    im.set_array(zi)
    return im,

path = "../test/data/"
dirs = os.listdir(path)
for file in dirs:
   print(file)

   



# I control the range of my colorbar by removing data 
# outside of my range of interest
#zmin = 293
#zmax = 380
#zi[(zi<zmin) | (zi>zmax)] = None

#im = plt.imshow(zi,animated=True)

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
plt.show()
