from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels as ker

import funciones_con_histogramas as histo


file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'

# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel1 = "kernel_poliangular_corregido.npz"
archivo1 = np.load(file_path + file_kernel1)
kernel1 = archivo1['kernel']

rAxis     = archivo1['rAxis']
thetaAxis = archivo1['thetaAxis']

file_kernel2 = "kernel_poliangular.npz"
archivo2 = np.load(file_path + file_kernel2)
kernel2 = archivo2['kernel']

# ker.plotPolarContour_full(kernel1, rAxis, thetaAxis, 4)

# ker.plotPolarContour_full(kernel2, rAxis, thetaAxis, 4)

error_rel = abs(kernel1-kernel2)/kernel2

plt.figure()
plt.imshow(error_rel)
plt.show()

print(error_rel.min(),error_rel.max())
