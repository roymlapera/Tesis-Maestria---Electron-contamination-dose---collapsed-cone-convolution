import sys
import matplotlib.pyplot as plt
import numpy as np

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Genero los pesos angulares

path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Espacios de fases/'

path_angulos = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'

filename =  'angulos_region_central.txt'

# filename_pablo  =  'espectro_incidencia.txt'

bins = 45

maxEnergy = 6.0

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Guardo los datos del EF en un archivo txt

# particulas_salteadas = 10000000 

# xinf = -5 #cm
# xsup = 5 #cm
# yinf = -5 #cm
# ysup = 5 #cm

# angulos, weights = ef.DevuelveAngulosDeArchivo(path_EF, maxEnergy, particulas_salteadas, xinf, xsup, yinf, ysup)

# histo.Guarda_datos_angulos(path_angulos + filename,angulos,weights)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

angulos, weights = histo.Abre_datos_angulos(path_angulos + filename)

#Correccion por angulo solido (disminuimos la contribucion de cada angulo por un factor 2*pi*sin theta)

weights = weights/(2*np.pi*np.sin(angulos))

fig, ax = plt.subplots()

histograma, bin_edges, patches = plt.hist(angulos, bins=bins, histtype='step', normed=True, color='b', weights=weights)

histograma /= histograma.sum()*np.pi/90      #Nose por que es necesario normalizar de nuevo, capaz tiene que ver el hecho de que le pasos pesos (weights)

plt.xlabel('Angulo director respecto al eje Z', fontsize=14)

plt.grid(True)

plt.show()

histo.Guarda_histograma(path_angulos + 'histograma_angulos_region_central.txt', path_angulos + 'histograma_angulos_region_central_bin_edges.txt', histograma, bin_edges)



