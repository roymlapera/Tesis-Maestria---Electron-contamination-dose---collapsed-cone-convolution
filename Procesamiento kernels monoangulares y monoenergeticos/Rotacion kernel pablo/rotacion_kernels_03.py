import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

from scipy import interpolate

from norm_check import norm_check
from plot_polar import plotPolarContour

# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'
directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_histogramas as histo

def CreaKernelAsterisco(kernel,rAxis,thetaAxis,theta_p):
	r_axis = (rAxis + np.append(0.0, rAxis[:-1]))/2 # en cm
	theta_axis = (thetaAxis + np.append(0.0, thetaAxis[:-1]))/2 # en cm

	f = interpolate.RectBivariateSpline(r_axis, theta_axis, kernel)

	u = np.array([np.sin(theta_axis), np.zeros(len(theta_axis)), np.cos(theta_axis)]).T    # ESTE ES ELL VECTOR TK

	A = np.zeros_like(kernel)

	for p in range(360):

	    phi_p = np.radians(p)

	    v = np.array([np.sin(theta_p)*np.cos(phi_p), np.sin(theta_p)*np.sin(phi_p), np.cos(theta_p)]).T  #ESTE ES EL VECTOR S1 PRIMA
	    # print(v.shape)

	    alfa = np.arccos(u.dot(v))   #alfa es el angulo phi sub t del paper

	    # print(alfa.shape)

	    Alfa, RAxis = np.meshgrid(alfa, r_axis)
	    A += f.ev(RAxis, Alfa)

	A = A/360

	return A

def CreaKernelPoliangular(kernel,rAxis,thetaAxis,delta_alpha,angulo,peso_angular):

	kernel_poliangular = np.zeros(kernel.shape)

	for i in range(len(peso_angular)):
		kernel_rot = CreaKernelAsterisco(kernel,rAxis,thetaAxis,angulo[i])
		kernel_poliangular += peso_angular[i] * kernel_rot
		if (i%10==0):
			plotPolarContour(kernel_rot, rAxis, np.degrees(thetaAxis), 4, drawType= 'hemi')
		print(str(100*(i+1)/len(peso_angular)) + ' %')

	delta_beta = 2*np.pi/360

	return kernel_poliangular*delta_alpha*delta_beta/(2*np.pi)

# ----------------------------------------------------------------------------------------------------------------------------

data = np.load('kernel_polienergetico_roy_final.npz')

kernel = data['kernel']
rAxis = data['rAxis']
thetaAxis = np.radians(data['thetaAxis'])

archivo_peso_angular = 'pesos_angulares_FE_reg_central_90bines_c_corregir_ecinetica.txt'
peso_angular = histo.Abre_archivo_pesos_angulares(archivo_peso_angular)

# peso_angular = np.ones(len(peso_angular))

delta_alpha = (np.pi/2)/90

alpha = np.arange(delta_alpha/2,np.pi/2,delta_alpha)
print(alpha)

peso_angular = peso_angular[:len(alpha)]

peso_angular /= peso_angular.sum()*delta_alpha
# print(peso_angular.sum()*delta_alpha)

# plt.figure()
# plt.plot(alpha,peso_angular,'b.')
# plt.show()

# norm, kernel_norm = norm_check(kernel, rAxis, np.degrees(thetaAxis))
# print(norm)

# kernel /= norm

# norm, kernel_norm = norm_check(kernel, rAxis, np.degrees(thetaAxis))
# print(norm)

# A = CreaKernelPoliangular(kernel,rAxis,thetaAxis,delta_alpha,alpha[:88],peso_angular[:88])

# A = CreaKernelAsterisco(kernel,rAxis,thetaAxis,np.radians(45))

# norm, kernel_norm = norm_check(A, rAxis, np.degrees(thetaAxis))
# print(norm*360)

# plotPolarContour(A, rAxis, np.degrees(thetaAxis), 4, drawType= 'hemi')

# np.savez("kernel_poliangular_roy_final.npz", kernel = A, rAxis = rAxis, thetaAxis = np.degrees(thetaAxis))

# ------------------------------------------------------------------------------------------------------------------------------------

# data1 = np.load('kernel_polienergetico_roy_final.npz')
# data2 = np.load('kernel_poliangular_roy_final_c_corr.npz')
# data3 = np.load('kernel_poliangular_roy_final_s_corr.npz')

# kernel1 = data1['kernel']
# kernel2 = data2['kernel']
# kernel3 = data3['kernel']

# rAxis = data1['rAxis']
# thetaAxis = data1['thetaAxis']

# norm1, kernel1_norm = norm_check(kernel1, rAxis, thetaAxis)
# norm2, kernel2_norm = norm_check(kernel2, rAxis, thetaAxis)
# norm3, kernel3_norm = norm_check(kernel3, rAxis, thetaAxis)

# plt.figure()
# plt.imshow((abs(kernel2-kernel1)/kernel1)[:20,:90], vmin=0, vmax=1)
# plt.figure()
# plt.imshow((abs(kernel3-kernel1)/kernel1)[:20,:90], vmin=0, vmax=1)
# plt.show()

# plt.figure()
# for i in range(18):
# 	plt.semilogy(rAxis,kernel[:,5*i],'r.-')
# plt.ylabel('Kernel')
# plt.xlabel('Radio [cm]')
# plt.show()

# ------------------------------------------------------------------------------------------------------------------------------------

data = np.load('kernel_poliangular_roy_final_v2_7.npz')

kernel = data['kernel']
rAxis = data['rAxis']
thetaAxis = data['thetaAxis']

plotPolarContour(kernel, rAxis, thetaAxis, 4.0, drawType= 'hemi')

