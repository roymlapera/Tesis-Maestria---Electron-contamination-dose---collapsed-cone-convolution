# Pasa kernel a cartesianas sin interpolacion v4

from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

# -----------------------------------------------------------------------------------------------------------------------------------------------

path_kernels = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'
# # path_regiones = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'

tamano_voxel_pequeno = 0.0005
tamano_voxel_grande = 0.1

archivo = np.load(path_kernels + 'kernel_poliangular_reg_central_cart_01_001.npz')
kernel1 = archivo['kernel']

archivo = np.load(path_kernels + 'kernel_poliangular_reg_central_cart_01_0001.npz')
kernel2 = archivo['kernel']

archivo = np.load(path_kernels + 'kernel_poliangular_reg_central.npz')
kernel_completo = archivo['kernel']

th = np.arange(0,180,1)+0.5

kernel3 = np.zeros(kernel_completo.shape[0])

print(kernel3.shape[0])

for i in range(kernel_completo.shape[0]):
	kernel3[i] = np.interp(0,th,kernel_completo[i,:])

archivo = np.load(path_kernels + 'kernel_MC_1000M_01cm.npz')
kernel_MC = archivo['dosis']

l = int(np.cbrt(kernel_MC.shape[0]))
kernel_MC = kernel_MC.reshape([l,l,l])
kernel_MC = kernel_MC[kernel_MC.shape[0]//2+1:,:,:]

kernel1 /= kernel1.sum()*tamano_voxel_grande**3
kernel2 /= kernel2.sum()*tamano_voxel_grande**3

kernel_MC /= kernel_MC.sum()*tamano_voxel_grande**3

eje_z1 = np.arange(1,len(kernel1[kernel1.shape[0]//2,kernel1.shape[1]//2,:])+1)/10
eje_z2 = np.arange(1,len(kernel2[kernel1.shape[0]//2,kernel1.shape[1]//2,:])+1)/10

eje_z3 = np.array( list(np.arange(0,2,0.02)-0.01) + list(np.arange(2,5,0.05)-0.025) + list(np.arange(5,10,0.1)-0.05) )

plt.figure()
plt.plot(eje_z1,kernel1[kernel1.shape[0]//2,kernel1.shape[1]//2,:], 'b.-', label='Cartesianas - tamano_vox_chico = 0.01 cm')
plt.plot(eje_z2,kernel2[kernel2.shape[0]//2,kernel2.shape[1]//2,:], 'm.-', label='Cartesianas - tamano_vox_chico = 0.001 cm')
plt.plot(eje_z3,kernel3, 'g.-', label='Esfericas')
plt.plot(eje_z2,kernel_MC[:,kernel_MC.shape[1]//2,kernel_MC.shape[1]//2], 'r.-', label='MC - Generado con DOSXYZnrc')

plt.xlabel('Profundidad [cm]')
plt.ylabel('Fraccion de energia depositada por unidad de volumen [1/cm3]')
plt.yscale('log')
plt.legend()
plt.show()