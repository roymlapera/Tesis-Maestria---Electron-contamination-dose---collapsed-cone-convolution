from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def spherical2cartesian(r, th, phi, grid, x, y, z, order=3):

    # Build relationship between Cartesian and spherical coordinates.
    X, Y, Z = np.meshgrid(x,y,z)
    new_r = np.sqrt(X*X+Y*Y+Z*Z)
    new_th = np.arccos(Z/new_r)
    new_phi = np.arctan2(Y, X)

    # Find these values for the input grid
    ir = interp1d(r, np.arange(len(r)), bounds_error=False)
    ith = interp1d(th, np.arange(len(th)))
    iphi = interp1d(phi, np.arange(len(phi)))

    new_ir = ir(new_r.ravel())
    new_ith = ith(new_th.ravel())
    new_iphi = iphi(new_phi.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    # Interpolate to Cartesian coordinates.
    return map_coordinates(grid, np.array([new_ir, new_ith, new_iphi]),
                            order=order).reshape(new_r.shape)


# Build 3D arrays for spherical coordinates.
r, th, phi = np.mgrid[0:201,0:201,0:201]

r = r/20.0               # r goes from 0 to 10.
th = th/200.0*np.pi         # Theta goes from 0 to pi
phi = phi/200.0*2*np.pi     # Phi goes from 0 to 2pi

# Density is spherically symmetric.  Only depends on r.
density = np.exp(-(r*np.cos(th))**2/0.5)

# Build ranges for function
r = np.linspace(0,200,200)
th = np.linspace(0,np.pi,200)
phi = np.linspace(-np.pi,np.pi,200)
x = np.linspace(-50,50,50)
y = np.linspace(-50,50,50)
z = np.linspace(-50,50,50)

# print (density)

# Map to Cartesian coordinates.
density = spherical2cartesian(r, th, phi, density, x, y, z, order=3)

# Plot density, now in spherical coordinates.
# not in cartesian coordinates.
plt.figure()
plt.imshow(density[:,:,10], cmap=cm.hot)
plt.show()


