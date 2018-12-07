import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
#import matplotlib.animation as animation


# Load data from CSV
dat = np.genfromtxt('GenericName002.ext', delimiter=' ',skip_header=0)
X_dat = dat[:,0]
Y_dat = dat[:,1]
Z_dat = dat[:,2]

# Convert from pandas dataframes to numpy arrays
X, Y, Z, = np.array([]), np.array([]), np.array([])
for i in range(len(X_dat)):                  
        X = np.append(X,X_dat[i])
        Y = np.append(Y,Y_dat[i])
        Z = np.append(Z,Z_dat[i])

print X,Y,Z

# create x-y points to be used in heatmap
xi = np.linspace(X.min(),X.max(),10)
yi = np.linspace(Y.min(),Y.max(),10)

# Z is a matrix of x-y values
zi = griddata((X, Y),Z,(xi[None,:], yi[:,None]), method='cubic')

# I control the range of my colorbar by removing data 
# outside of my range of interest
zmin = 293
zmax = 380
zi[(zi<zmin) | (zi>zmax)] = None

# Create the contour plot
CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
                  vmax=zmax, vmin=zmin)
plt.colorbar()  
plt.show()
