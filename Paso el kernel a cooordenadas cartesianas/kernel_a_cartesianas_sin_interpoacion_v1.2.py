# Pasa kernel a cartesianas sin interpolacion

from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math

def rebin(arr, new_shape):
    if(arr.shape[0] == new_shape[0] and arr.shape[1] == new_shape[1] and arr.shape[2] == new_shape[2]):
        return arr
    else:
        shape = (new_shape[0], arr.shape[0] // new_shape[0], new_shape[1], arr.shape[1] // new_shape[1], new_shape[2], arr.shape[2] // new_shape[2])
        return arr.reshape(shape).sum(axis=1).sum(axis=2).sum(axis=3)

def Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis):
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

    # print('La normalizo para que su integral valga 1')
    # kernel3D_esf /= Frac_energ_esf.sum()
    # Frac_energ_esf = kernel3D_esf*volume
    # print('La integral del kernel3D_esf en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume

def esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, xmin, xmax, ymin, ymax, zmin, zmax, tamano_voxel_grande):

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    nro_vox_gra_x = round(abs(xmax-xmin)/tamano_voxel_grande)
    nro_vox_gra_y = round(abs(ymax-ymin)/tamano_voxel_grande)
    nro_vox_gra_z = round(abs(zmax-zmin)/tamano_voxel_grande)

    nro_vox_peq_x = round(abs(xmax-xmin)/tamano_voxel_pequeno)
    nro_vox_peq_y = round(abs(ymax-ymin)/tamano_voxel_pequeno)
    nro_vox_peq_z = round(abs(zmax-zmin)/tamano_voxel_pequeno)

    print(nro_vox_peq_x,nro_vox_peq_y,nro_vox_peq_z,nro_vox_gra_x,nro_vox_gra_y,nro_vox_gra_z)

    # #Inicializo los puntos que quiero en el sistema cartesiano HD

    x = linspace(xmin+tamano_voxel_pequeno/2,xmax-tamano_voxel_pequeno/2,nro_vox_peq_x)
    # print(x)
    y = linspace(ymin+tamano_voxel_pequeno/2,ymax-tamano_voxel_pequeno/2,nro_vox_peq_y)
    z = linspace(zmin+tamano_voxel_pequeno/2,zmax-tamano_voxel_pequeno/2,nro_vox_peq_z)

    #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

    X, Y, Z = np.meshgrid(x,y,z)

    #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

    new_r = np.sqrt(X*X+Y*Y+Z*Z)
    new_th = np.arccos(-Z/new_r)    
    new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

    #Aplico las funciones obtenidas para encontrar qué indices de la grilla del kernel le corresponde a (X,Y,Z) 

    new_ir = ir(new_r.ravel())
    new_ith = ith(new_th.ravel())
    new_iphi = iphi(new_phi.ravel())

    #Corrijo por si hay algunos valores mal interpolados

    new_ir[new_r.ravel() > r.max()] = 209
    new_ir[new_r.ravel() < r.min()] = 0

    new_ith[new_th.ravel() > th.max()] = 180
    new_ith[new_th.ravel() < th.min()] = 0

    new_iphi[new_phi.ravel() > phi.max()] = 360
    new_iphi[new_phi.ravel() < phi.min()] = 0

    frac_energ_3D_cart_HD, a, b = mgrid[0:nro_vox_peq_x//2:0.5,0:nro_vox_peq_y//2:0.5,0:nro_vox_peq_z//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

    new_ir   = new_ir.round().astype(int32).reshape(frac_energ_3D_cart_HD.shape)
    new_ith  = new_ith.round().astype(int32).reshape(frac_energ_3D_cart_HD.shape)
    new_iphi = new_iphi.round().astype(int32).reshape(frac_energ_3D_cart_HD.shape)

    for i in range(nro_vox_gra_x):
    	for j in range(nro_vox_gra_y):
    		for k in range(nro_vox_gra_z):
    			frac_energ_3D_cart_HD[i][j][k] = Frac_energ_esf[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ] * tamano_voxel_pequeno**3 / volume[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ]

    frac_energ_3D_cart = rebin(frac_energ_3D_cart_HD,(nro_vox_gra_x,nro_vox_gra_y,nro_vox_gra_z))

    # frac_energ_3D_cart = frac_energ_3D_cart_HD

    return frac_energ_3D_cart/tamano_voxel_grande**3

def Concatena_subregiones1(path_regiones):
    subregiones = []

    for i in range(3):
        for j in range(3):
            for k in range(3):
                archivo = np.load(path_regiones + 'kernel_cart_reg1_'+str(i)+str(j)+str(k)+'.npz')
                subregion = archivo['kernel']
                subregiones.append(subregion)
    
    aux = []

    for i in range(len(subregiones)//3):
        aux.append(np.concatenate((subregiones[3*i],subregiones[3*i+1],subregiones[3*i+1]), axis=2))

    aux2 = []

    for i in range(len(aux)//3):
        aux2.append(np.concatenate((aux[3*i],aux[3*i+1],aux[3*i+1]), axis=1))

    for i in range(len(aux2)//3):
        reg1 = np.concatenate((aux2[3*i],aux2[3*i+1],aux2[3*i+1]), axis=0)

    return reg1


file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

kernel = archivo['kernel']
kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis'] * np.pi/180
phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

Frac_energ_esf, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    Frac_energ_esf[:,:,k] = kernel

print(Frac_energ_esf.max())

volume = Integral_vol_del_kernel(Frac_energ_esf,rAxis,thetaAxis,phiAxis)

Frac_energ_esf = Frac_energ_esf * volume

# -----------------------------------------------------------------------------------------------------------------------

# Genero funciones para encontrar indices de ls puntos (X,Y,Z) en la grilla del kernel en esfericas

r = [0] + rAxis
th = [0] + thetaAxis
phi = phiAxis

#Obtengo funciones a trozos ir, ith y iphi obtenidas de r, th y phi

ir = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value='extrapolate')
ith = interp1d(th, np.arange(len(th)), bounds_error=False, fill_value='extrapolate')
iphi = interp1d(phi, np.arange(len(phi)), bounds_error=False, fill_value='extrapolate')

# -----------------------------------------------------------------------------------------------------------------------

# # Mapeo a coordenadas cartesianas.

nro_voxels_pequenos = 126
nro_voxels_grandes = 21

kernel3D_cart, a, b = mgrid[0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5,0:nro_voxels_pequenos//2:0.5]

path_regiones = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'
# path_regiones = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'

# Region 1

tamano_voxel_pequeno = 0.0025
tamano_voxel_grande = 0.05

# for i in range(3):
#     for j in range(3):
#         for k in range(3):
#             kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, -0.025+i*1.05, 1.025+i*1.05, -0.025+j*1.05, 1.025+j*1.05, 0.025-k*1.05, -1.025-k*1.05, tamano_voxel_grande)
#             np.savez(path_regiones + 'kernel_cart_reg1_'+str(i)+str(j)+str(k)+'.npz', kernel = kernel3D_cart)

i = 0
j = 0
k = 0

kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, -0.025+i*1.05, 1.025+i*1.05, -0.025+j*1.05, 1.025+j*1.05, 0.025-k*1.05, -1.025-k*1.05, tamano_voxel_grande)

figure()
for i in range(0,60,2):
  imshow(squeeze(kernel3D_cart[i,:,:]), cmap=cm.hot, norm=LogNorm())
  print(i)
  show()

# tamano_voxel_pequeno = 0.025
# tamano_voxel_grande = tamano_voxel_pequeno*2

# # Region 2

# kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, -0.025, 10.025, -0.025, 10.025, -3.125, -10.025, tamano_voxel_grande)
# np.savez(path_regiones + 'kernel_cart_reg2'+'.npz', kernel = kernel3D_cart)

# # Region 3

# kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, -0.025, 10.025, 3.125, 10.025, 0.025, -3.125, tamano_voxel_grande)
# np.savez(path_regiones + 'kernel_cart_reg3'+'.npz', kernel = kernel3D_cart)

# # Region 3

# kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, Frac_energ_esf, tamano_voxel_pequeno, 3.125, 10.025, -0.025, 3.125, 0.025, -3.125, tamano_voxel_grande)
# np.savez(path_regiones + 'kernel_cart_reg4'+'.npz', kernel = kernel3D_cart)

# # -----------------------------------------------------------------------------------------------------------------------

#Abro regiones y concateno

# path_regiones = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'
# # path_regiones = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'

# reg1 = Concatena_subregiones1(path_regiones)

# print(reg1.shape)

# figure()
# imshow(squeeze(reg1[:,15,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**6))

# figure()
# imshow(squeeze(reg1[:,30,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**6))

# figure()
# imshow(squeeze(reg1[:,45,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**6))
# show()

# print('la integral del kernel3D_cart en la esfera de radio 10 da: ' + str(kernel3D_cart.sum()*tamano_voxel_grande**3))

# # -----------------------------------------------------------------------------------------------------------------------

# # Grafico el kernel en cartesianas

# figure()
# imshow(squeeze(kernel3D_cart[:,nro_voxels_grandes//2,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**6))
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



