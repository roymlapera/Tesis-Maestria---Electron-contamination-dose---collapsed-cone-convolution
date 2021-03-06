# Pasa kernel a cartesianas sin interpolacion v4

from pylab import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels as ker

def Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi):
    '''Importa los valores de kernel y la informacion espacial a partir de un archivo NPZ. Los devuelve en forma de matriz y vectores, respectivamente.'''
    archivo = np.load(file_path + npz_file)

    kernel = archivo['kernel']

    rAxis     = archivo['rAxis']
    thetaAxis = archivo['thetaAxis'] * np.pi/180
    phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

    kernel3D_esf, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

    print(kernel3D_esf.shape)

    # Asigno el kernel para los angulos phi, por haber simetria de revolucion

    kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

    # kernel /= len(phiAxis)  # me quedo con la fraccion de energia de un gajo

    # ker.plotPolarContour_full(kernel, rAxis, thetaAxis, 4)

    for k in range(nro_conos_phi):
        kernel3D_esf[:,:,k] = kernel

    volume = Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

    I = (kernel3D_esf*volume).sum()
    kernel3D_esf /= I

    Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

    return kernel3D_esf, rAxis, thetaAxis, phiAxis

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

    print('La integral del kernel3D_esf en la esfera de radio ? da: ' + str(Frac_energ_esf.sum()))

    return volume

def Devuelve_frac_energ_en_vox_grande(Frac_energ_esf, volume, ir, ith, iphi, tamano_voxel_pequeno, tamano_voxel_grande, limite, x_vox_grande, y_vox_grande, z_vox_grande):
    '''Calcula la fraccion energetica de un voxel centrado en (x_vox_grande, y_vox_grande, z_vox_grande) subdividiendolo en voxeles pequenos. 
    El fantoma considerad tiene tamano 2*limite. La fraccion de energia en el voxel grande es la suma de la fraccion de enegia de los voxeles pequenos. +
    Un voxel pequeno esta dentro del voxel grande si su centro esta en el grande.'''

    # El kernel esferico solo tiene valores para un radi menor a limite:

    if( x_vox_grande**2+y_vox_grande**2+z_vox_grande**2 < limite**2):
        # #Inicializo los puntos que quiero en el sistema cartesiano 

        pos_primer_vox_pequeno = -tamano_voxel_grande/2+tamano_voxel_pequeno/2
        pos_ultimo_vox_pequeno =  tamano_voxel_grande/2-tamano_voxel_pequeno/2

        nro_voxels_pequenos = int(tamano_voxel_grande/tamano_voxel_pequeno)

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

        #Aplico las funciones obtenidas para encontrar qu?? indices de la grilla del kernel le corresponde a (X,Y,Z) 

        new_ir = ir(new_r.ravel())
        new_ith = ith(new_th.ravel())
        new_iphi = iphi(new_phi.ravel())

        #Corrijo por si hay algunos valores mal interpolados

        new_ir[new_r.ravel() > 10] = 209
        new_ir[new_r.ravel() < 0] = 0

        new_ith[new_th.ravel() > np.pi] = 180-1
        new_ith[new_th.ravel() < np.pi/180] = 0

        new_iphi[new_phi.ravel() > -np.pi] = 360-1
        new_iphi[new_phi.ravel() < np.pi] = 0

        new_ir   = new_ir.round().astype(int32).reshape(X.shape)
        new_ith  = new_ith.round().astype(int32).reshape(X.shape)
        new_iphi = new_iphi.round().astype(int32).reshape(X.shape)

        #Sumo la fraccion de energia depositada en el voxel grande

        frac_energ_voxel = 0

        for i in range(nro_voxels_pequenos):
            for j in range(nro_voxels_pequenos):
                for k in range(nro_voxels_pequenos):
                    frac_energ_voxel += Frac_energ_esf[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ] * tamano_voxel_pequeno**3 / volume[ new_ir[i][j][k] ][ new_ith[i][j][k] ][ new_iphi[i][j][k] ]
    else:
        frac_energ_voxel = 0

    # print(frac_energ_voxel)

    return frac_energ_voxel

def esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, tamano_voxel_pequeno, xmin, xmax, ymin, ymax, zmin, zmax, tamano_voxel_grande):

    # Genero funciones para encontrar indices de ls puntos (X,Y,Z) en la grilla del kernel en esfericas

    r = [0] + rAxis
    th = [0] + thetaAxis
    phi = phiAxis

    #Obtengo funciones a trozos ir, ith y iphi obtenidas de r, th y phi

    ir = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value='extrapolate')
    ith = interp1d(th, np.arange(len(th)), bounds_error=False, fill_value='extrapolate')
    iphi = interp1d(phi, np.arange(len(phi)), bounds_error=False, fill_value='extrapolate')

    # -------------------------------------------------------------------------------------------------------------------------------------------------------

    # #Inicializo los puntos que quiero en el sistema cartesiano

    nro_vox_gra_x = round(abs(xmax-xmin)/tamano_voxel_grande)
    nro_vox_gra_y = round(abs(ymax-ymin)/tamano_voxel_grande)
    nro_vox_gra_z = round(abs(zmax-zmin)/tamano_voxel_grande)

    # #Inicializo los puntos que quiero en el sistema cartesiano

    x = np.linspace(xmin+tamano_voxel_grande/2,xmax-tamano_voxel_grande/2,nro_vox_gra_x)
    y = np.linspace(ymin+tamano_voxel_grande/2,ymax-tamano_voxel_grande/2,nro_vox_gra_y)
    z = np.linspace(zmin-tamano_voxel_grande/2,zmax+tamano_voxel_grande/2,nro_vox_gra_z)

    print(x)
    print(y)
    print(z)

    print(x.shape,y.shape,z.shape)

    #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

    X, Y, Z = np.meshgrid(y,x,z)

    #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

    new_r = np.sqrt(X*X+Y*Y+Z*Z)
    new_th = np.arccos(-Z/new_r)    
    new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

    new_th[np.argwhere(new_r==0)]=0

    #Aplico las funciones obtenidas para encontrar qu?? indices de la grilla del kernel le corresponde a (X,Y,Z) 

    new_ir = ir(new_r.ravel())
    new_ith = ith(new_th.ravel())
    new_iphi = iphi(new_phi.ravel())

    #Corrijo por si hay algunos valores mal interpolados

    new_ir[new_r.ravel() >= r.max()] = 209
    new_ir[new_r.ravel() <= r.min()] = 0

    new_ith[new_th.ravel() >= th.max()] = 180-1
    new_ith[new_th.ravel() <= th.min()] = 0

    new_iphi[new_phi.ravel() >= phi.max()] = 360-1
    new_iphi[new_phi.ravel() <= phi.min()] = 0

    # frac_energ_3D_cart, a, b = mgrid[0:nro_vox_gra_x//2:0.5,0:nro_vox_gra_y//2:0.5,0:nro_vox_gra_z//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32
    frac_energ_3D_cart, a, b = mgrid[0:nro_vox_gra_x,0:nro_vox_gra_y,0:nro_vox_gra_z].astype(float64)

    # print(frac_energ_3D_cart.shape)

    new_ir   = new_ir.round().astype(int32).reshape(frac_energ_3D_cart.shape)
    new_ith  = new_ith.round().astype(int32).reshape(frac_energ_3D_cart.shape)
    new_iphi = new_iphi.round().astype(int32).reshape(frac_energ_3D_cart.shape)

    # ------------------------------------------------------------------------------------------------------------------------------------

    volume = Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

    Frac_energ_esf = kernel3D_esf*volume

    for i in range(nro_vox_gra_x):
        for j in range(nro_vox_gra_y):
            for k in range(nro_vox_gra_z):
                frac_energ_3D_cart[i][j][k] = Devuelve_frac_energ_en_vox_grande(Frac_energ_esf, volume, ir, ith, iphi, tamano_voxel_pequeno, tamano_voxel_grande, limite, X[i][j][k], Y[i][j][k], Z[i][j][k])
                sys.stdout.write("\r{:.2f}".format(100*(k+j*nro_vox_gra_z+i*nro_vox_gra_z*nro_vox_gra_y)/(nro_vox_gra_x*nro_vox_gra_y*nro_vox_gra_z)) + ' %')
                sys.stdout.flush()  

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
        aux.append(np.concatenate((subregiones[3*i],subregiones[3*i+1],subregiones[3*i+2]), axis=2))

    aux2 = []

    for i in range(len(aux)//3):
        aux2.append(np.concatenate((aux[3*i],aux[3*i+1],aux[3*i+2]), axis=1))

    for i in range(len(aux2)//3):
        reg1 = np.concatenate((aux2[3*i],aux2[3*i+1],aux2[3*i+2]), axis=0)

    return reg1

# -----------------------------------------------------------------------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

npz_file = "kernel_poliangular_roy_final_final.npz"

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

kernel3D_esf, rAxis, thetaAxis, phiAxis = Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# # Mapeo a coordenadas cartesianas.

limite = 10 #cm

# tamano_voxel_grande  = 0.05
# tamano_voxel_pequeno = 0.002

# kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, tamano_voxel_pequeno, -0.775, 0.775, -0.075, 0.075, 0, -1, tamano_voxel_grande)

path_regiones = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'
# # path_regiones = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'

# # Region 1

tamano_voxel_pequeno = 0.001
tamano_voxel_grande = 0.1

kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, tamano_voxel_pequeno, -0.05, 3.05, -0.05, 3.05, 0, -3.1, tamano_voxel_grande)
np.savez(path_regiones + 'kernel_cart_roy_final_final_reg1.npz', kernel = kernel3D_cart)

# figure()
# imshow(kernel3D_cart[:,0,:], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))
# show()

# archivo = np.load(path_regiones + 'fullWater_50M_2125MeV_cart_reg1.npz')
# reg1 = archivo['kernel']

# for i in range(kernel3D_cart.shape[0]):
#     for j in range(kernel3D_cart.shape[1]):
#         for k in range(kernel3D_cart.shape[2]):
#             reg1[i][j][k] = kernel3D_cart[i][j][k]

# reg1 = np.concatenate((reg1[:,::-1,:],reg1[:,1:,:]), axis=1)
# reg1 = np.concatenate((reg1[::-1,:,:],reg1[1:,:,:]), axis=0)

# print(reg1.shape)

# print(reg1.sum()*tamano_voxel_grande**3)

# # reg1 /= reg1.sum()*tamano_voxel_grande**3

# # print(reg1.sum()*tamano_voxel_grande**3)

# np.savez(path_regiones + 'fullWater_50M_2125MeV_cart_corregido.npz', kernel = reg1)

# figure()
# imshow(squeeze(reg1[:,reg1.shape[1]//2,:]).T, cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))

# figure()
# imshow(squeeze(reg1[reg1.shape[1]//2,:,:]).T, cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))

# figure()
# imshow(squeeze(reg1[:,:,2]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))
# show()

# ----------------------------------------------------------------------------------------------------------------------------------------

# archivo = np.load(path_regiones + 'kernel_poliangular_reg_central_cart_corregido_a_05.npz')
# reg1 = archivo['kernel']

# archivo = np.load(path_regiones + 'kernel_MC_10M_1.npz')
# kernel_MC = archivo['dosis']
# l = int(np.cbrt(kernel_MC.shape[0]))
# kernel_MC = kernel_MC.reshape([l,l,l])

# kernel_MC = kernel_MC[kernel_MC.shape[0]//2:,:,:]

# reg1 /= reg1[reg1.shape[0]//2,reg1.shape[1]//2,30]
# kernel_MC /= kernel_MC[30,kernel_MC.shape[1]//2,kernel_MC.shape[1]//2]

# # kernel_MC /= kernel_MC.sum()*tamano_voxel_grande**3

# plt.figure()
# plt.plot(reg1[reg1.shape[0]//2,reg1.shape[1]//2,:], 'b.')
# plt.plot(kernel_MC[:,kernel_MC.shape[1]//2,kernel_MC.shape[1]//2], 'r.')
# plt.yscale('log')
# plt.show()

# ----------------------------------------------------------------------------------------------------------------------------------------

# figure()
# for i in range(0,reg1.shape[1],10):
#     imshow(squeeze(reg1[:,i,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-5, vmax=10**4))
#     print(i)
#     show()

# for i in range(3):
#     for j in range(3):
#         for k in range(3):
#             kernel3D_cart = esfericas_a_cartesianas_sin_interpolar(rAxis, thetaAxis, phiAxis, kernel3D_esf, tamano_voxel_pequeno, -0.025+i*1.05, 1.025+i*1.05, -0.025+j*1.05, 1.025+j*1.05, 0.025-k*1.05, -1.025-k*1.05, tamano_voxel_grande)
#             np.savez(path_regiones + 'kernel_cart_reg1_'+str(i)+str(j)+str(k)+'.npz', kernel = kernel3D_cart)

# i = 0
# j = 0
# k = 0

# figure()
# for i in range(0,kernel3D_cart.shape[1],5):
#     imshow(squeeze(kernel3D_cart[:,i,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-5, vmax=10**4))
#     print(i)
#     show()

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
# path_regiones = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/Por regiones/'

# reg1 = Concatena_subregiones1(path_regiones)

# print(reg1.min(),reg1.max(),reg1.shape)

# figure()
# for i in range(0,reg1.shape[0],5):
#   imshow(squeeze(reg1[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-5), vmax=10**4))
#   print(i)
#   show()

# print('la integral del kernel3D_cart en la esfera de radio 10 da: ' + str(kernel3D_cart.sum()*tamano_voxel_grande**3))

# # -----------------------------------------------------------------------------------------------------------------------

# # Grafico el kernel en cartesianas

# figure()
# imshow(squeeze(kernel3D_cart[:,nro_voxels_grandes//2,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**6))
# show()

# #Grafico varios cortes del kernel

# figure()
# for i in range(nro_voxels//2-10,nro_voxels//2+10,1):
#   imshow(squeeze(kernel3D[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-7), vmax=10**4), extent=(-limite,limite,-limite,limite))
#   print(i)
#   show()

# # -----------------------------------------------------------------------------------------------------------------------

# #Guardo el kernel en cordenadas cartesianas. Tamano de la matriz: nro_voxels*nro_voxels*nro_voxels

# np.savez("kernel_cart_200a200.npz", kernel = kernel3D_cart)
