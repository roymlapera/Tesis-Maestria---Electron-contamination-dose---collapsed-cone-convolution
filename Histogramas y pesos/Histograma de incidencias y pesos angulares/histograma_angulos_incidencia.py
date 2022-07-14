import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'
# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Genero los pesos angulares

# path_EF = 'C:/Users/OWner/Documents/Facu/Bariloche 2/'
path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion de EFs/'

# path_angulos = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Histogramas y pesos/Histograma de incidencias y pesos angulares/Angulos/'
path_angulos = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts\Histogramas y pesos/Histograma de incidencias y pesos angulares/Angulos/'

filename_F       = 'angulos_vs_fluencia_energetica_reg_central.txt'
filename_FE      = 'angulos_vs_fluencia_energetica_reg_central_ecinetica.txt'

bins = 90
maxEnergy = 6.0

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Guardo los datos del EF en un archivo txt

tamano_pack_particulas = -1

xinf = -5 #cm
xsup = 5 #cm
yinf = -5 #cm
ysup = 5 #cm

nro_EF = 0
# angulos, weights = ef.DevuelveAngulosDeArchivo(path_EF, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup,nro_EF)
# ef.Guarda_data(path_angulos + filename_FE,angulos,weights)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

angulos, weights       = ef.Abre_energias(path_angulos + filename_FE)
angulos_v, weights_v   = ef.Abre_energias(path_angulos + filename_F)

#la energia del EF es la energia total, le resto la energia en reposo

Area_EF = 30*30 #cm2
weights /= Area_EF
# weights_v  /= Area_EF

plt.figure()
histograma, bin_edges, patches = plt.hist(angulos, bins=bins, histtype='step', normed=True, color='g', weights=weights, label='Sin corrección')
histograma2, bin_edges, patches = plt.hist(angulos, bins=bins, histtype='step', normed=True, color='orange', weights=weights/(2*np.pi*np.sin(angulos)), label='Con corrección')
# histograma_v, bin_edges_v, patches_v = plt.hist(angulos_v, bins=bins, histtype='step', normed=True, color='m', weights=weights_v/(2*np.pi*np.sin(angulos_v)))
plt.ylabel('dΨ/dα normalizado [1/rad]', fontsize=14)
plt.xlabel('Ángulo cenital α [rad]', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()

# print(histograma.sum()*np.pi/2/90)
# print(histograma_v.sum()*np.pi/2/90)

np.savetxt(path_angulos + 'pesos_angulares_FE_reg_central_90bines_s_corregir_ecinetica.txt', histograma)
np.savetxt(path_angulos + 'pesos_angulares_FE_reg_central_90bines_c_corregir_ecinetica.txt', histograma2)

# ------------------------------------------------------------------------------------------------------------------------