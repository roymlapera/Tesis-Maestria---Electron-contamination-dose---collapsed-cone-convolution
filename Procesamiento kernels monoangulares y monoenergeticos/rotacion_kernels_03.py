import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy import interpolate

from norm_check import norm_check
from plot_polar import plotPolarContour


data = np.load('kernel_polienergetico_roy_24.npz')

kernel = data['kernel']
rAxis = data['rAxis']
thetaAxis = np.radians(data['thetaAxis'])

# print(kernel.sum())

# norm, kernel_norm = norm_check(kernel, rAxis, np.degrees(thetaAxis))

# print(norm)
# print(kernel.shape)

# plotPolarContour(kernel, rAxis, np.degrees(thetaAxis), 4.0, drawType= 'hemi')

r_axis = (rAxis + np.append(0.0, rAxis[:-1]))/2 # en cm
theta_axis = (thetaAxis + np.append(0.0, thetaAxis[:-1]))/2 # en cm

f = interpolate.RectBivariateSpline(r_axis, theta_axis, kernel)

u = np.array([np.sin(theta_axis), np.zeros(len(theta_axis)), np.cos(theta_axis)]).T

# print(u.shape)

# plt.figure()
# plt.plot(u[:,0], u[:,2],'.')
# plt.axis('equal')
# plt.show()

A = np.zeros_like(kernel)

theta_p = np.radians(45)

for p in range(360):

    phi_p = np.radians(p)

    v = np.array([np.sin(theta_p)*np.cos(phi_p), np.sin(theta_p)*np.sin(phi_p), np.cos(theta_p)]).T
    # print(v.shape)

    alfa = np.arccos(u.dot(v))

    Alfa, RAxis = np.meshgrid(alfa, r_axis)
    A += f.ev(RAxis, Alfa)

A = A/360



norm, kernel_norm = norm_check(A, rAxis, np.degrees(thetaAxis))

print(norm)

print(kernel_norm.shape)

plotPolarContour(A, rAxis, np.degrees(thetaAxis), 4.0, drawType= 'hemi')