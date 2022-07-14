import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.special
from scipy import optimize
import scipy.ndimage as ndimage

# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'
directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def funcion_ajuste_FE(sigma_x,sigma_y,l,A):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = a
    d = b

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    funcion_modela_endurecimiento = A*(x**2+y**2)+1
    return lambda x,y: cte*funcion_modela_endurecimiento*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )

def funcion_ajuste_FE_simple(sigma_x,sigma_y,l):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = a
    d = b

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return lambda x,y: cte*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )

def fit_function(funcion_ajuste, data, params):
    """Devuelve los parametros de ajuste bidimiensional de la funcion_ajuste"""
    xy = (np.indices(data.shape) - 60)/2
    errorfunction = lambda p: np.ravel(funcion_ajuste(*p)(*xy) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# ----------------------------------------------- OBTENGO EL HISTOGRAMA BIDIMENSIONAL DE LA FENERGETICA ----------------------------------------------------------------

bins = 121

maxEnergy = 6.0

# path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion del EF para el kernel de comparacion/'
path_EF = 'C:/Users/OWner/Documents/Facu/Bariloche 2/'

#Guardo los datos del EF en un archivo txt

fenergetica_archivo = 'fenergetica.txt'

#Analizo todas las part

# fenergetica, PosX, PosY = ef.PuntosXY_Fluenciaenergetica_particulas(path_EF, maxEnergy,0)

# histo.Guarda_datos_fenergetica(fenergetica_archivo, fenergetica, PosX, PosY)

# FE, x, y = histo.Abre_datos_fenergetica(fenergetica_archivo)

# Fenergetica, xedges, yedges = histo.Grafica_PosXY_fenergetica(x, y, FE, bins)

# np.savez("histo_fluencia_energetica_todos_los_electrones.npz", Fenergetica = Fenergetica, xedges = xedges, yedges = yedges)

# ---------------------------------------------------------------- PASO A UNIDADES DE MEV/CM2 -------------------------------------------------------------------------

limite_fantoma = 30 #cm   #el fantoma va de -limite a +limite

archivo_FE = np.load("histo_fluencia_energetica_todos_los_electrones.npz")

Fenergetica = archivo_FE['Fenergetica']

area_pixel = (limite_fantoma*2/bins)**2

Fenergetica /= area_pixel   #ahora la fluencia enrgetica esta en MeV/cm2

# ---------------------------------------------------------------- FILTRADO -------------------------------------------------------------------------

sigma_de_filtrado = 3
Fenergetica_filtrada = ndimage.gaussian_filter(Fenergetica, sigma=sigma_de_filtrado)

# plt.plot(np.arange(bins),Fenergetica_filtrada[bins//2])
# plt.plot(np.arange(bins),Fenergetica[bins//2], 'r.')
# plt.show()

# plt.imshow(Fenergetica)
# plt.show()
# plt.imshow(Fenergetica_filtrada, extent=[-30,30,-30,30])
# plt.xlabel('Posicion X [cm]')
# plt.ylabel('Posicion Y [cm]')
# plt.show()

Fenergetica_filtrada /= Fenergetica_filtrada.sum()*area_pixel

print(Fenergetica_filtrada.sum()*area_pixel)

# ---------------------------------------------------------------- AJUSTE -------------------------------------------------------------------------

x = np.linspace(-limite_fantoma,limite_fantoma,bins)
y = np.linspace(-limite_fantoma,limite_fantoma,bins)

X, Y = np.meshgrid(x,y)



sigma_x = 8

sigma_y = 9

l = 10

# alpha = 2.8E-2

seeds = [sigma_x,sigma_y,l]

params = fit_function(funcion_ajuste_FE_simple,Fenergetica_filtrada, seeds)

print(params)

fig = plt.figure()
ax = plt.subplot(111)
# plt.imshow(Fenergetica, cmap=cm.gist_earth_r)
ax.plot(x,Fenergetica_filtrada[:,61], 'b.', label='FE filtrada X')
ax.plot(x,funcion_ajuste_FE_simple(params[0],params[1],params[2])(x,0),'k-', label='ajuste')
ax.plot(x,funcion_ajuste_FE_simple(seeds[0],seeds[1],seeds[2])(x,0),'g', label='ajuste a ojo')
ax.legend()
plt.show()

# fig = plt.figure()
# ax = plt.subplot(111)
# # plt.imshow(Fenergetica, cmap=cm.gist_earth_r)
# ax.plot(x,Fenergetica_filtrada[61,:], 'r.', label='FE filtrada Y')
# ax.plot(x,funcion_ajuste_FE_simple(params[0],params[1],params[2])(0,x),'k-', label='ajuste')
# ax.plot(x,funcion_ajuste_FE_simple(seeds[0],seeds[1],seeds[2])(x,0),'m', label='ajuste a ojo')
# ax.legend()
# plt.show()

# # fig = plt.figure()
# # ax = plt.subplot(111)

# # # plt.imshow(Fenergetica, cmap=cm.gist_earth_r)
# # ax.plot(x,Fenergetica_filtrada[:,61]/5.28E-5, 'r.', label='FE filtrada')
# # ax.plot(x,funcion_ajuste_FE3(params[0],params[1])(x,0), label='ajuste')
# # ax.plot(x,funcion_ajuste_FE3(seeds[0],seeds[1])(x,0), label='ajuste a ojo')
# # ax.legend()
# # plt.show()

# ------------------------------------------------------------------------------------------------

#Plots de perfiles de ajuste de FE

# fig = plt.figure()
# ax = plt.subplot(111)

# # plt.imshow(Fenergetica, cmap=cm.gist_earth_r)
# ax.plot(x,Fenergetica_filtrada[:,61]/5.28E-5, 'r.', label='FE filtrada')
# ax.plot(x,funcion_ajuste_FE3(params[0],params[1])(x,0), label='ajuste')
# ax.plot(x,funcion_ajuste_FE3(seeds[0],seeds[1])(x,0), label='ajuste a ojo')
# ax.legend()
# plt.show()
