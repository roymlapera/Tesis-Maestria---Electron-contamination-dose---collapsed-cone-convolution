from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.cm as cm

directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'
# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def Grafica_fenergetica(energia, angulo, FE, bins_E, bins_A):    #ACA poner PosX PosY y FE
 
    hist, xedges, yedges = np.histogram2d(energia-0.511, angulos, bins=[bins_E,bins_A], range=[[0, 6], [0, np.pi/2]], normed=False, weights=FE/(2*np.pi*np.sin(angulos)))

    hist = np.rot90(hist)

    print(hist.shape)

    # for i in range(hist.shape[0]):
    # 	norm = hist[i,:].sum()*6/bins_E
    # 	if (norm!=0):
    # 		hist[i,:] /= norm
    # 	print(hist[i,:].sum()*6/bins_E)

    hist[:,0] = 0
    hist[:,1] = 0
    hist[:,2] = 0

    norm = hist.sum()*(6/bins_E)*(np.pi/2/bins_E)

    hist /= norm

    # fig = plt.figure()
 
    # plt.imshow(hist, extent=[0,6,0,np.pi/2], interpolation='None', cmap=cm.GnBu_r, aspect='auto')

    # plt.xlabel('Energía cinética [MeV]')
    # plt.ylabel('Ángulo de incidencia [rad]')

    # plt.tight_layout()
    # plt.show()

    return hist, xedges, yedges

# ----------------------------------------------- OBTENGO EL HISTOGRAMA BIDIMENSIONAL DE LA FENERGETICA ----------------------------------------------------------------

fenergetica_archivo = 'fenergetica.txt'
# fenergetica_archivo = 'fenergetica_reg_central.txt'

bins_E = 24
bins_A = 9

FE, energia, angulos = histo.Abre_datos_fenergetica(fenergetica_archivo)

Fenergetica, xedges, yedges = Grafica_fenergetica(energia, angulos, FE, bins_E, bins_A)

print(Fenergetica.shape)

# ----------------------------------------------------------------------------------------------------------------------

Es = (xedges[:-1] + xedges[1:])/2

delta_A = np.pi/2/bins_A

angulos = np.arange(delta_A/2,np.pi/2,delta_A)
angulos = angulos[::-1]

print(angulos)

# Fenergetica = Fenergetica[::-1,:]

i = 0

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for c, z in zip(['b', 'g', 'y', 'orange','r', 'orange', 'y', 'g','b'], range(9)):
    ax.bar(Es, Fenergetica[z,:], width=6/bins_E ,zs=angulos[z], zdir='y', color=c, ec=c, alpha=0.5)
    np.savetxt('espectro' + str(i) + '.txt', Fenergetica[z,:])
    i += 1

ax.set_xlabel('Energía cinética [MeV]')
ax.set_ylabel('Ángulo cenital α [rad]')
ax.set_zlabel('dΨ/dE normalizada [1/MeV]')

plt.show()