import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.special
from scipy import optimize
import scipy.ndimage as ndimage
import matplotlib.ticker as mtick

# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'
directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def FE(alpha,sigma,C):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""

    l = 10 #cm

    a = -l
    b = l
    c = a
    d = b

    sigma_x = sigma
    sigma_y = sigma

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return lambda x,y: cte*alpha*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )+C


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

fenergetica_archivo = 'fenergetica2.txt'

#Analizo todas las part

# fenergetica, PosX, PosY = ef.PuntosXY_Fluenciaenergetica_particulas(path_EF, maxEnergy,0)

# histo.Guarda_datos_fenergetica(fenergetica_archivo, fenergetica, PosX, PosY)

# FE, x, y = histo.Abre_datos_fenergetica(fenergetica_archivo)

# Fenergetica, xedges, yedges = histo.Grafica_PosXY_fenergetica(x, y, FE, bins)

# np.savez("histo_fluencia_energetica_todos_los_electrones2.npz", Fenergetica = Fenergetica, xedges = xedges, yedges = yedges)

# ---------------------------------------------------------------- PASO A UNIDADES DE MEV/CM2 -------------------------------------------------------------------------

limite_fantoma = 30 #cm   #el fantoma va de -limite a +limite

archivo_FE = np.load("histo_fluencia_energetica_todos_los_electrones.npz")

Fenergetica = archivo_FE['Fenergetica']

area = (limite_fantoma/bins)**2

Fenergetica /= area   #ahora la fluencia enrgetica esta en MeV/cm2

# ---------------------------------------------------------------- FILTRADO -------------------------------------------------------------------------

sigma_de_filtrado = 2
Fenergetica_filtrada = ndimage.gaussian_filter(Fenergetica, sigma=sigma_de_filtrado)

# plt.plot(np.arange(bins),Fenergetica_filtrada[bins//2])
# plt.plot(np.arange(bins),Fenergetica[bins//2], 'r.')
# plt.show()

# plt.imshow(Fenergetica)
# plt.show()
plt.imshow(Fenergetica_filtrada, extent=[-30,30,-30,30])
plt.xlabel('Posición X [cm]')
plt.ylabel('Posición Y [cm]')
plt.show()

print(Fenergetica_filtrada.sum()*area**2)

# Fenergetica_filtrada /= Fenergetica_filtrada.sum()*area**2



# print(Fenergetica_filtrada.sum()*area_pixel)

# ---------------------------------------------------------------- AJUSTE -------------------------------------------------------------------------

x = np.linspace(-limite_fantoma,limite_fantoma,bins)
y = np.linspace(-limite_fantoma,limite_fantoma,bins)

X, Y = np.meshgrid(x,y)

sigma = 7.439
alpha = 0.556*0.025
C     = 8.72E-5*0.025

seeds = [alpha,sigma,C]

params = fit_function(FE,Fenergetica_filtrada, seeds)
params[0] /=1.015
print(params)

# fig = plt.figure()
# ax = plt.subplot(111)
# ax.plot(x,Fenergetica_filtrada[:,bins//2], 'b.', label='FE filtrada X')
# ax.plot(x,FE(params[0],params[1],params[2])(x,0),'k-', label='Ajuste')
# ax.legend()
# ax.set_ylabel('Fluencia energética [MeV/$\mathregular{cm^2}$]')
# ax.set_xlabel('Posición X [cm]')
# ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
# ax.grid()

# fig = plt.figure()
# ax = plt.subplot(111)
# ax.plot(y,Fenergetica_filtrada[bins//2,:], 'r.', label='FE filtrada Y')
# ax.plot(y,FE(params[0],params[1],params[2])(0,y),'k-', label='Ajuste')
# ax.set_ylabel('Fluencia energética [MeV/$\mathregular{cm^2}$]')
# ax.set_xlabel('Posición Y [cm]')
# ax.legend()
# ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
# ax.grid()

# plt.show()

F_energ = np.zeros((bins,bins))

for i in range(bins):
	for j in range(bins):
		F_energ[i,j] = FE(params[0],params[1],params[2])(x[i],y[j])

plt.figure()
plt.imshow(F_energ)
plt.show()


