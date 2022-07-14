# Pasa kernel a cartesianas sin interpolacion

from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys

def Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi):
    '''Importa los valores de kernel y la informacion espacial a partir de un archivo NPZ. Los devuelve en forma de matriz y vectores, respectivamente.'''
    archivo = np.load(file_path + npz_file)

    kernel = archivo['kernel']
    kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

    rAxis     = archivo['rAxis']
    thetaAxis = archivo['thetaAxis'] * np.pi/180
    phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

    kernel3D_esf, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

    # Asigno el kernel para los angulos phi, por haber simetria de revolucion

    for k in range(nro_conos_phi):
        kernel3D_esf[:,:,k] = kernel

    return kernel3D_esf, rAxis, thetaAxis, phiAxis

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1],
             new_shape[2], arr.shape[2] // new_shape[2])
    return arr.reshape(shape).sum(axis=1).sum(axis=2).sum(axis=3)

def Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion muestra por pantalla el resultado de la integral
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

    return volume

def Devuelve_frac_energ_en_vox_grande(Frac_energ_esf, volume, ir, ith, iphi, nro_voxels_pequenos, nro_voxels_grandes, limite, x_vox_grande, y_vox_grande, z_vox_grande):
    '''Calcula la fraccion energetica de un voxel centrado en (x_vox_grande, y_vox_grande, z_vox_grande) subdividiendolo en voxeles pequenos. 
    El fantoma considerad tiene tamano 2*limite. La fraccion de energia en el voxel grande es la suma de la fraccion de enegia de los voxeles pequenos. +
    Un voxel pequeno esta dentro del voxel grande si su centro esta en el grande.'''

    # El kernel esferico solo tiene valores para un radi menor a limite:

    if( x_vox_grande**2+y_vox_grande**2+z_vox_grande**2 < limite**2):
        # #Inicializo los puntos que quiero en el sistema cartesiano HD

        tamano_voxel_grande  = limite*2/nro_voxels_grandes
        tamano_voxel_pequeno = tamano_voxel_grande/nro_voxels_pequenos

        pos_primer_vox_pequeno = -tamano_voxel_grande/2+tamano_voxel_pequeno/2
        pos_ultimo_vox_pequeno =  tamano_voxel_grande/2-tamano_voxel_pequeno/2

        x = linspace( pos_primer_vox_pequeno , pos_ultimo_vox_pequeno , nro_voxels_pequenos ) + x_vox_grande
        y = linspace( pos_primer_vox_pequeno , pos_ultimo_vox_pequeno , nro_voxels_pequenos ) + y_vox_grande
        z = linspace( pos_primer_vox_pequeno , pos_ultimo_vox_pequeno , nro_voxels_pequenos ) + z_vox_grande 

        #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

        X, Y, Z = np.meshgrid(x,y,z)

        #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

        new_r = np.sqrt(X*X+Y*Y+Z*Z)
        new_th = np.arccos(-Z/new_r)    
        new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

        new_th[np.argwhere(new_r == 0)] = 0

        #Aplico las funciones obtenidas para encontrar qué indices de la grilla del kernel le corresponde a (X,Y,Z) 

        new_ir = ir(new_r.ravel())
        new_ith = ith(new_th.ravel())
        new_iphi = iphi(new_phi.ravel())

        #Corrijo por si hay algunos valores mal interpolados

        new_ir[new_r.ravel() > 10] = 210-1
        new_ir[new_r.ravel() < 0] = 0

        new_ith[new_th.ravel() > np.pi] = 180-1
        new_ith[new_th.ravel() < np.pi/180] = 0

        new_iphi[new_phi.ravel() > -np.pi] = 360-1
        new_iphi[new_phi.ravel() < np.pi] = 0

        new_ir   = new_ir.round().astype(int32).reshape(X.shape)
        new_ith  = new_ith.round().astype(int32).reshape(X.shape)
        new_iphi = new_iphi.round().astype(int32).reshape(X.shape)

        #Sumo la fraccion de energia depositada en el voxel grande

        semi_nro_vox = nro_voxels_pequenos//2
        frac_energ_voxel = 0

        for i in range(nro_voxels_pequenos):
            for j in range(nro_voxels_pequenos):
                for k in range(nro_voxels_pequenos):
                    frac_energ_voxel += Frac_energ_esf[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ] * tamano_voxel_pequeno**3 / volume[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ]
    else:
        frac_energ_voxel = 0

    # print(frac_energ_voxel)

    return frac_energ_voxel

    

def esfericas_a_cartesianas(kernel3D_esf, rAxis, thetaAxis, phiAxis, nro_voxels_pequenos, nro_voxels_grandes, limite):

    #Chequeo normalizacion del kernel en coordenadas esfericas, de paso me queda hecha la matriz de volumenes para calcular la fraccion energetica

    volume = Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

    Frac_energ_esf = kernel3D_esf * volume

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    #Obtengo funciones a trozos ir, ith y iphi obtenidas de r, th y phi

    r   = [0] + rAxis
    th  = [0] + thetaAxis
    phi = phiAxis

    ir = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value='extrapolate')
    ith = interp1d(th, np.arange(len(th)), bounds_error=False, fill_value='extrapolate')
    iphi = interp1d(phi, np.arange(len(phi)), bounds_error=False, fill_value='extrapolate')

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    # #Inicializo los puntos que quiero en el sistema cartesiano de voxeles grandes

    tamano_voxel_grande  = limite*2/nro_voxels_grandes

    pos_primer_vox_grande = -limite+tamano_voxel_grande/2
    pos_ultimo_vox_grande =  limite-tamano_voxel_grande/2

    x = linspace( pos_primer_vox_grande , pos_ultimo_vox_grande , nro_voxels_grandes )
    y = linspace( pos_primer_vox_grande , pos_ultimo_vox_grande , nro_voxels_grandes )
    z = linspace( pos_primer_vox_grande , pos_ultimo_vox_grande , nro_voxels_grandes )

    #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

    X, Y, Z = np.meshgrid(x,y,z)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    kernel3D_cart = np.zeros((nro_voxels_grandes,nro_voxels_grandes,nro_voxels_grandes))

    for l in range(nro_voxels_grandes):
        for m in range(nro_voxels_grandes):
            for n in range(nro_voxels_grandes):
                kernel3D_cart[l][m][n] = Devuelve_frac_energ_en_vox_grande(Frac_energ_esf, volume, ir, ith, iphi, nro_voxels_pequenos, nro_voxels_grandes, limite, X[l][m][n], Y[l][m][n], Z[l][m][n])
        sys.stdout.write("\r{:.2f}".format(((l*nro_voxels_grandes**2)/nro_voxels_grandes**3)*100) + '%\n')
        sys.stdout.flush()   
                
    print('la integral del kernel3D_cart en la esfera de radio 10 da: ' + str(kernel3D_cart.sum()))

    return kernel3D_cart/tamano_voxel_grande**3 

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

npz_file = "kernel_poliangular.npz"

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

kernel3D_esf, rAxis, thetaAxis, phiAxis = Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# # Mapeo a coordenadas cartesianas.

limite = 10 #cm
nro_voxels_grandes  = 51  #debe ser impar para que el kernel tenga un voxel central
nro_voxels_pequenos = 11  #debe ser impar para centrar bien el voxel grande

#Inicializo kernel3D_cart asi para que sea float64

kernel3D_cart, a, b = mgrid[0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5]

kernel3D_cart = esfericas_a_cartesianas(kernel3D_esf, rAxis, thetaAxis, phiAxis, nro_voxels_pequenos, nro_voxels_grandes, limite)

# # -----------------------------------------------------------------------------------------------------------------------

# # Grafico un corrte del kernel en cartesianas

figure()
imshow(squeeze(kernel3D_cart[:,:,nro_voxels_grandes//2-1]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**1))
show()

# #Grafico varios cortes del kernel

# figure()
# for i in range(nro_voxels//2-10,nro_voxels//2+10,1):
# 	imshow(squeeze(kernel3D[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-7), vmax=10**4), extent=(-limite,limite,-limite,limite))
# 	print(i)
# 	show()

# # -----------------------------------------------------------------------------------------------------------------------

# #Guardo el kernel en cordenadas cartesianas. Tamano de la matriz: nro_voxels*nro_voxels*nro_voxels

# np.savez("kernel_cart_200a200.npz", kernel = kernel3D_cart)



