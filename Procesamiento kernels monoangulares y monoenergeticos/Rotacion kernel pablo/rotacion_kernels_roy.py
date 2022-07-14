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

def CreaKernelPoliangular(kernels,rAxis,thetaAxis,angulo,bins_A):

	kernel_poliangular = np.zeros(kernels[0].shape)

	for i in range(len(kernels)):
		kernel_rot = CreaKernelAsterisco(kernels[i],rAxis,thetaAxis,angulo[i])
		# plotPolarContour(kernel_rot, rAxis, np.degrees(thetaAxis), 4, drawType= 'hemi')
		kernel_poliangular += kernel_rot
		print(str(100*(i+1)/bins_A) + ' %')

	return kernel_poliangular/(np.pi**2)

# ----------------------------------------------------------------------------------------------------------------------------

data0 = np.load('kernel_polienergetico_roy_final0.npz')
data1 = np.load('kernel_polienergetico_roy_final1.npz')
data2 = np.load('kernel_polienergetico_roy_final2.npz')
data3 = np.load('kernel_polienergetico_roy_final3.npz')
data4 = np.load('kernel_polienergetico_roy_final4.npz')
data5 = np.load('kernel_polienergetico_roy_final5.npz')
data6 = np.load('kernel_polienergetico_roy_final6.npz')
data7 = np.load('kernel_polienergetico_roy_final7.npz')
data8 = np.load('kernel_polienergetico_roy_final8.npz')

rAxis = data0['rAxis']
thetaAxis = np.radians(data0['thetaAxis'])


kernel6 = data6['kernel']
kernel7 = data7['kernel']
kernel8 = data8['kernel']
kernel0 = kernel6*0
kernel1 = kernel6*0
kernel2 = kernel6*0
kernel3 = kernel6*0
kernel4 = kernel6*0
kernel5 = kernel6*0



kernels = [kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,kernel8]

norm, kernel_norm = norm_check(kernel0, rAxis, np.degrees(thetaAxis))
print(norm)

# ---------------------------------------------------------------------------------------------------------------------------------------

bins_A=9
delta_alpha = np.pi/2/bins_A

alpha = np.arange(delta_alpha/2,np.pi/2,delta_alpha)

print(alpha*180/np.pi)
kernels = kernels[::-1]

# bins_A = 3
# alpha = alpha[:bins_A]
# kernels = kernels[:bins_A]


# archivo_peso_angular = 'pesos_angulares_FE_reg_central_90bines_c_corregir_ecinetica.txt'
# peso_angular = histo.Abre_archivo_pesos_angulares(archivo_peso_angular)

# plt.figure()
# plt.plot(peso_angular,'b.')

# peso_angular = peso_angular.reshape( (10,9) ).sum(axis=0)
# print(peso_angular)

# plt.figure()
# plt.plot(peso_angular,'b.')
# plt.show()



# peso_angular = peso_angular[:len(alpha)]

# peso_angular /= peso_angular.sum()*delta_alpha
# # print(peso_angular.sum()*delta_alpha)

# plt.figure()
# plt.plot(alpha,peso_angular,'b.')
# plt.show()

A = CreaKernelPoliangular(kernels,rAxis,thetaAxis,alpha,bins_A)

norm, kernel_norm = norm_check(A, rAxis, np.degrees(thetaAxis))
print(norm)

plotPolarContour(A, rAxis, np.degrees(thetaAxis), 4, drawType= 'hemi')

np.savez("kernel_poliangular_roy_final_v30.npz", kernel = A, rAxis = rAxis, thetaAxis = np.degrees(thetaAxis))

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

# data = np.load('kernel_poliangular_roy_final_prueba.npz')

# kernel = data['kernel']
# rAxis = data['rAxis']
# thetaAxis = data['thetaAxis']

# plotPolarContour(kernel, rAxis, thetaAxis, 4.0, drawType= 'hemi')

