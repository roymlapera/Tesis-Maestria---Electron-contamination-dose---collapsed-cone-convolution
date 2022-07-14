import sys
import numpy as np

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels as ker

archivo = np.load("kernel_poliangular_final.npz")

kernel_poli = archivo['kernel']
rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis']

ker.plotPolarContour_full(kernel_poli, rAxis, thetaAxis, 0.5)


