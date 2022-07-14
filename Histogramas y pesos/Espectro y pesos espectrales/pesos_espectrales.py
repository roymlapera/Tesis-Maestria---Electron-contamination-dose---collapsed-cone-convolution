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

path_EF = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Histogramas por region/'
# path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion de EFs/'

path_energia = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Histogramas por region/Espectro de fluencia y de fluencia energetica/'
# path_energia = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Histogramas y pesos/Espectro y pesos espectrales/'

filename_F_reg_central      = 'pesos_espectrales.txt'

bins = 24
maxEnergy = 6.0

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Guardo los datos del EF en un archivo txt

tamano_pack_particulas = -1

xinf = -5 #cm
xsup = 5 #cm
yinf = -5 #cm
ysup = 5 #cm

nro_EF = 0
# energias, weights = ef.DevuelveEnergiasDeArchivo(path_EF, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup,nro_EF)
# ef.Guarda_data(path_energia + filename_F_reg_central,energias,weights)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

archivo_pesos_viejos = 'pesos_espectro_fluencia_energ_24bines.txt'
pesos_viejos = np.loadtxt(archivo_pesos_viejos)

# print(pesos_viejos.sum()*6/24)

energias, weights       = ef.Abre_energias(filename_F_reg_central)

# energias -= 0.511 #MeV

# #la energia del EF es la energia total, le resto la energia en reposo

Area_EF = 30*30 #cm2
weights /= Area_EF

histograma, bin_edges, patches = plt.hist(energias, range=[0,6] , bins=bins, histtype='step', normed=False, color='b', weights=weights)
histograma2, bin_edges, patches = plt.hist(energias, range=[0,6] , bins=72, histtype='step', normed=False, color='b', weights=weights)

histograma /= histograma.sum()*6/24
histograma2 /= histograma2.sum()*6/72

plt.figure()
plt.bar(np.arange(6/24/2,6,6/24), histograma, 6/24, color = 'None', edgecolor = 'g',label='24 bines', linewidth=3)
plt.bar(np.arange(6/72/2,6,6/72), histograma2, 6/72, color = 'None', edgecolor = 'orange',label='72 bines')

plt.ylabel('dΨ/dE normalizado [1/MeV]', fontsize=14)
plt.xlabel('E cinética [MeV]', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()

print(histograma.shape)

# np.savetxt(path_energia + 'pesos_FE_24bines.txt', histograma)

