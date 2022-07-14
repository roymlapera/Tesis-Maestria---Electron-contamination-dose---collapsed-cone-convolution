# Pasa kernel a cartesianas sin interpolacion

from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1],
             new_shape[2], arr.shape[2] // new_shape[2])
    return arr.reshape(shape).sum(axis=1).sum(axis=2).sum(axis=3)

def Integral_vol_del_kernel(kernel,rAxis,thetaAxis,phiAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion imprime por pantalla el resultado de la integral
    y devuelve la matriz de volumenes con ejes segun rAxis thetaAxis y phiAxis.'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < nro_radios-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel3D_esf.shape)

    for j in range(len(phiAxis)):
        for k in range(len(thetaAxis)):
            volume[:,k,j] = 2*np.pi/len(phiAxis) * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3) 

    Frac_energ_esf = kernel3D_esf*volume

    print('La integral del kernel3D_esf en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    print('La normalizo para que su integral valga 1')
    kernel3D_esf /= Frac_energ_esf.sum()
    Frac_energ_esf = kernel3D_esf*volume
    print('La integral del kernel3D_esf en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume

def esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, nro_voxels_pequenos, nro_voxels_grandes, limite):

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    #Chequeo normalizacion del kernel en coordenadas esfericas, de paso me queda hecha la matriz de volumenes

    volume = Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    # tamano_voxel_pequeno = limite/(nro_voxels_pequenos/2)
    # tamano_voxel_grande = limite/(nro_voxels_grandes/2)

    # kernel3D_cart_HD, a, b = mgrid[0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

    # # #Inicializo los puntos que quiero en el sistema cartesiano HD

    # x = linspace(-limite,limite,nro_voxels_pequenos)
    # y = linspace(-limite,limite,nro_voxels_pequenos)
    # z = linspace(-limite,limite,nro_voxels_pequenos)

    # #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

    # X, Y, Z = np.meshgrid(x,y,z)

    # #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

    # new_r = np.sqrt(X*X+Y*Y+Z*Z)
    # new_th = np.arccos(-Z/new_r)    
    # new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

    # r = [0] + rAxis
    # th = [0] + thetaAxis
    # phi = phiAxis

    # #Obtengo funciones a trozos ir, ith y iphi obtenidas de r, th y phi

    # ir = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value='extrapolate')
    # ith = interp1d(th, np.arange(len(th)), bounds_error=False, fill_value='extrapolate')
    # iphi = interp1d(phi, np.arange(len(phi)), bounds_error=False, fill_value='extrapolate')

    # #Aplico las funciones obtenidas para encontrar qué indices de la grilla del kernel le corresponde a (X,Y,Z) 

    # new_ir = ir(new_r.ravel())
    # new_ith = ith(new_th.ravel())
    # new_iphi = iphi(new_phi.ravel())

    # #Corrijo por si hay algunos valores mal interpolados

    # new_ir[new_r.ravel() > r.max()] = len(r)-1
    # new_ir[new_r.ravel() < r.min()] = 0

    # new_ith[new_th.ravel() > th.max()] = len(th)-1
    # new_ith[new_th.ravel() < th.min()] = 0

    # new_iphi[new_phi.ravel() > phi.max()] = len(phi)-1
    # new_iphi[new_phi.ravel() < phi.min()] = 0

    # new_ir   = new_ir.round().astype(int32).reshape(kernel3D_cart_HD.shape)
    # new_ith  = new_ith.round().astype(int32).reshape(kernel3D_cart_HD.shape)
    # new_iphi = new_iphi.round().astype(int32).reshape(kernel3D_cart_HD.shape)

    # semi_nro_vox = nro_voxels_pequenos//2

    # for i in range(nro_voxels_pequenos):
    # 	for j in range(nro_voxels_pequenos):
    # 		for k in range(nro_voxels_pequenos):
    # 			if (  (i-semi_nro_vox)**2 + (j-semi_nro_vox)**2 + (k-semi_nro_vox)**2 > (semi_nro_vox)**2  ):
    # 				kernel3D_cart_HD[i,j,k] = 0
    # 				continue
    # 			kernel3D_cart_HD[i][j][k] = Frac_energ_esf[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ] * tamano_voxel_pequeno**3 / volume[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ]
    
    # print('la integral del kernel3D_cart_HD en la esfera de radio 10 da: ' + str(kernel3D_cart_HD.sum()))

    # kernel3D_cart = rebin(kernel3D_cart_HD,(nro_voxels_grandes,nro_voxels_grandes,nro_voxels_grandes))

    return kernel3D_cart/tamano_voxel_grande**3 


file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360
nro_voxels_grandes = 201
nro_voxels_pequenos = 201


limite = 10 #cm

tamano_voxel_pequeno = limite/(nro_voxels_pequenos/2)
tamano_voxel_grande = limite/(nro_voxels_grandes/2)

kernel = archivo['kernel']
kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis'] * np.pi/180
phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

kernel3D_esf, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    kernel3D_esf[:,:,k] = kernel

# -----------------------------------------------------------------------------------------------------------------------

# # Mapeo a coordenadas cartesianas.

#Inicializo kernel3D_cart asi para que sea float64

kernel3D_cart, a, b = mgrid[0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5]

kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, nro_voxels_pequenos, nro_voxels_grandes, limite)

print('la integral del kernel3D_cart en la esfera de radio 10 da: ' + str(kernel3D_cart.sum()*tamano_voxel_grande**3))

# # -----------------------------------------------------------------------------------------------------------------------

# # Grafico el kernel en cartesianas

# figure()
# imshow(squeeze(kernel3D_cart[:,nro_voxels_grandes//2,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**1))
# show()

# #Grafico varios cortes del kernel

# figure()
# for i in range(nro_voxels//2-10,nro_voxels//2+10,1):
# 	imshow(squeeze(kernel3D[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-7), vmax=10**4), extent=(-limite,limite,-limite,limite))
# 	print(i)
# 	show()

# # -----------------------------------------------------------------------------------------------------------------------

# #Guardo el kernel en cordenadas cartesianas. Tamano de la matriz: nro_voxels*nro_voxels*nro_voxels

# np.savez("kernel_cart_200a200.npz", kernel = kernel3D_cart)



