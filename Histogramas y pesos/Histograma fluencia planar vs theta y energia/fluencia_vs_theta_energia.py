import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.special
from scipy import optimize
import scipy.ndimage as ndimage

# directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'
directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def Grafica_fenergetica(energia, angulo, FE, bins):    #ACA poner PosX PosY y FE
    fig = plt.figure()
 
    hist, xedges, yedges = np.histogram2d(energia-0.511, angulos, bins=bins, range=[[0, 6], [0, np.pi/2]], normed=False, weights=FE/(2*np.pi*np.sin(angulos)))

    hist = np.rot90(hist)

    # print(hist.shape)

    for i in range(hist.shape[0]):
    	norm = hist[i,:].sum()*6/bins
    	if (norm!=0):
    		hist[i,:] /= norm

    norm = hist.sum()*(6/bins)*(np.pi/2/bins)

    hist /= norm
 
    plt.imshow(hist, extent=[0,6,0,np.pi/2], interpolation='None', cmap=cm.GnBu_r, aspect='auto')

    plt.xlabel('Energía cinética [MeV]')
    plt.ylabel('Ángulo de incidencia [rad]')

    plt.tight_layout()
    plt.show()

    return hist, xedges, yedges

# ----------------------------------------------- EXTRAIGO LA INFO DEL EF ----------------------------------------------------------------

# path_EF = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Histogramas por region/'
path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion de EFs/'

fenergetica_archivo = 'fenergetica_no_planar.txt'

maxEnergy = 6.0

nro_EF = 0

xinf = -5 #cm
xsup = 5 #cm
yinf = -5 #cm
ysup = 5 #cm

fenergetica, energia,angulo = ef.PuntosEnergiaAngulo_Fluenciaenergetica_particulas(path_EF, maxEnergy, xinf, xsup, yinf, ysup, nro_EF)

histo.Guarda_datos_fenergetica(fenergetica_archivo, fenergetica, energia, angulo)

# print(energia.min(),energia.max())
# print(angulo.min(),angulo.max())
# print(fenergetica.min(),fenergetica.max())

# ----------------------------------------------- OBTENGO EL HISTOGRAMA BIDIMENSIONAL DE LA FENERGETICA ----------------------------------------------------------------

# bins = 60

# FE, energia, angulos = histo.Abre_datos_fenergetica(fenergetica_archivo)

# Fenergetica, xedges, yedges = Grafica_fenergetica(energia, angulos, FE, bins)


# ----------------------------------------------- CHECK ----------------------------------------------------------------

# np.savetxt('hist_sin_corregir.txt',Fenergetica)

# FE_c_corr = np.loadtxt('hist_con_corregir.txt')
# FE_s_corr = np.loadtxt('hist_sin_corregir.txt')
# espectro  = np.loadtxt('pesos_fluencia_100bines.txt')

# espectro_FE_c_corr = np.zeros( FE_c_corr[0,:].shape )
# espectro_FE_s_corr = np.zeros( FE_s_corr[0,:].shape )

# for i in range(espectro_FE_c_corr.shape[0]):
#     espectro_FE_c_corr[i] = FE_c_corr[:,i].sum()*np.pi/2/bins
#     espectro_FE_s_corr[i] = FE_s_corr[:,i].sum()*np.pi/2/bins

# espectro_FE_c_corr /= espectro_FE_c_corr.sum()*6/bins
# espectro_FE_s_corr /= espectro_FE_s_corr.sum()*6/bins

# energia1 = np.arange(6/bins/2,6,6/bins)
# energia2 = np.arange(6/100/2,6,6/100)

# espectro *=energia2
# espectro /= espectro.sum()*6/espectro.shape[0]



# # plt.figure()
# # plt.imshow(FE_c_corr)
# # plt.figure()
# # plt.imshow(FE_s_corr)
# plt.figure()
# plt.plot(energia2, espectro, 'b.')
# plt.plot(energia1, espectro_FE_c_corr, 'rx')
# plt.plot(energia1, espectro_FE_s_corr, 'g+')
# plt.xlabel('Energía [MeV]')
# plt.show()

