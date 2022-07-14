# ABRE PDDS MONOENERGETICOS

import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt


archivo_pesos = 'pesos_espectro_fluencia_energ_24bines.txt'
pesos_espectro = np.genfromtxt(archivo_pesos)

energias = np.arange(0.125,6,0.25)
print(energias)
print(pesos_espectro.sum()*6/24)

plt.figure()
plt.plot(energias,pesos_espectro,'r.')
plt.show()

# # ------------------------------------------------------- PDD ------------------------------------------------------- 

tam_vox_cart = 0.05#cm

xmin = -tam_vox_cart/2
xmax = tam_vox_cart/2
ymin = -tam_vox_cart/2
ymax = tam_vox_cart/2
zmin = 0
zmax = 1.3

nro_vox_gra_x = round(abs(xmax-xmin)/tam_vox_cart)
nro_vox_gra_y = round(abs(ymax-ymin)/tam_vox_cart)
nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)

# #Inicializo los puntos que quiero en el sistema cartesiano

X = np.linspace(xmin+tam_vox_cart/2,xmax-tam_vox_cart/2,nro_vox_gra_x)
Y = np.linspace(ymin+tam_vox_cart/2,ymax-tam_vox_cart/2,nro_vox_gra_y)
Z = np.linspace(zmin+tam_vox_cart/2,zmax-tam_vox_cart/2,nro_vox_gra_z)


prof_norm = 0.1466 #cm      #profundidad a la que uno quiere normalizar los PDDs

PDDs = []

for i in range(24):
	archivo_dosis_CC = np.load('dosis_conv_conocolapsado_v2_fullWater_50M_' + str(energias[i]).replace('.', '') + 'MeV.npz')
	dosis = archivo_dosis_CC['dosis']
	dosis = dosis[dosis.shape[0]//2,dosis.shape[1]//2,:]*pesos_espectro[i]*6/24
	PDDs.append(dosis)

# PDDs = np.array(PDDs)
PDDs = np.array(PDDs)

PDDs_acum = PDDs.cumsum(axis=0)


eje_z_CC = np.arange(zmin,zmax,tam_vox_cart) + tam_vox_cart/2

# dosis_CC_norm = np.interp(prof_norm,eje_z_CC,dosis_CC[dosis_CC.shape[0]//2,dosis_CC.shape[1]//2,:])

# # # # ----------------------------------------

# # PDD calculado en DOSXYZnrc de electrones de 2.125 MeV cuya incidencia es normal a la sup

# archivo_dosis_monoAmonoE = np.load("Varian-20-IAEA_monoangular_y_monoenergetica_electrones_1000M_paraPDD.npz")
# dosis_monoAmonoE = archivo_dosis_monoAmonoE['dosis']
# error_monoAmonoE = archivo_dosis_monoAmonoE['errs']

# dosis_monoAmonoE = dosis_monoAmonoE.reshape((34,11,11))

# error_monoAmonoE = error_monoAmonoE.reshape((34,11,11))

# dosis_monoAmonoE = dosis_monoAmonoE[:,dosis_monoAmonoE.shape[1]//2,dosis_monoAmonoE.shape[2]//2]
# error_monoAmonoE = error_monoAmonoE[:,error_monoAmonoE.shape[1]//2,error_monoAmonoE.shape[2]//2]

# eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
# dosis_monoAmonoE_norm    = np.interp(prof_norm,eje_z,dosis_monoAmonoE)

# # # # -------------------------------------------------------------------------------------------------------------- 



for i in range(24):
	plt.figure()
	plt.plot(eje_z_CC, 100*PDDs_acum[i,:]/PDDs_acum.max(), linewidth=0.5)

plt.xlabel('Profundidad [cm]')
plt.legend()
plt.show()