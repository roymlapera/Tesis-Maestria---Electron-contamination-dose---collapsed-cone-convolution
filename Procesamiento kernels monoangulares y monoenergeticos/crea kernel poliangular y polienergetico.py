import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'
# directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels as ker

import funciones_con_histogramas as histo

def CreaFracEnerg(kernel,rAxis,thetaAxis):
    '''A partir del kernel polienergetico (fracciopn de E depo por unidad de volumen esferico de revolucion) obtendo su distribucion de fraccion de E dep'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < 209:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel.shape)

    for k in range(len(thetaAxis)):
        volume[:,k] = 2*np.pi * (np.cos(ThetaAxis[k]*np.pi/180)-np.cos(thetaAxis[k]*np.pi/180))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    # print(volume[0,0])

    Frac_energ = kernel*volume

    print('La integral del kernel polienergetico en la esfera de radio 10 da: ' + str(Frac_energ.sum()))

    return Frac_energ

def Integral_vol_del_kernel(kernel,rAxis,thetaAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion imprime por pantalla el resultado de la integral
    y devuelve la matriz de volumenes con ejes segun rAxis thetaAxis y phiAxis.'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < 209:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel.shape)

    for k in range(len(thetaAxis)):
        volume[:,k] = 2*np.pi/360 * (np.cos(ThetaAxis[k]*np.pi/180)-np.cos(thetaAxis[k]*np.pi/180))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    # Frac_energ_esf = kernel*volume

    # print('La integral del kernel en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume

#-------------------------------------------------------------------------------- MAIN ------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'


archivo_pesos = 'histograma_angulos_region_central.txt'

NC = 180        #Nro de conos del kernel
NR = 210        #Nro de radios del kernel

#Creo los kernels rotados 

nc_rotacion = 2      #Guardo los kernels rotados cada 2 grados

archivo = np.load("kernel_polienergetico_pablo.npz")
# archivo = np.load("kernel_poliangular_reg_central_v2.npz")

kernel_poli = archivo['kernel']
rAxis     = archivo['rAxis'] 
thetaAxis = archivo['thetaAxis'] 

frac_energ = CreaFracEnerg(kernel_poli,rAxis,thetaAxis) 

volume = Integral_vol_del_kernel(kernel_poli,rAxis,thetaAxis)

kernel_poli = frac_energ/360/volume

# kernels = ker.Crea_lista_kernels_rotados(kernel_poli, NR, NC, nc_rotacion)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

# #Grafico un kernel para mostrar como se ve no mas

# kernel_poli = ker.GeneraKernel360(kernel_poli)

# kernel_rotado = ker.RotaKernel(kernel_poli,0,NR,NC)

# frac_energ = kernel_poli * ker.GeneraMatrizVolumenDeVoxelEsf(rAxis,thetaAxis)
# print(frac_energ.sum()*180)
# print(kernel_poli.sum())
# print(kernel_rotado.sum())

# ker.plotPolarContour_full_360(kernel_rotado, rAxis, thetaAxis, 4)

# np.savez("kernel_polienergetico.npz", kernel = kernel_rotado, rAxis = rAxis, thetaAxis = thetaAxis)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Importacion del peso angular y calculo de kernel poliangular

# pesos_angulares = histo.Abre_archivo_pesos_angulares(file_path + archivo_pesos)

# plt.figure()
# plt.plot(pesos_angulares,'b.')
# plt.show()

# kernel_poliangular = ker.CreaKernelPoliangular(kernels, pesos_angulares, NC, NR, nc_rotacion)

# print(kernel_poliangular.shape)

# kernel_poliangular = kernel_poliangular[:,:180]

# kernel_poliangular /= (kernel_poliangular*volume).sum()*360   #Fuerzo a que la norma sea 1

# print((kernel_poliangular*volume).sum()*360)

# np.savez("kernel_poliangular_reg_central.npz", kernel = kernel_poliangular, rAxis = rAxis, thetaAxis = thetaAxis)

# plt.figure()

# plt.imshow(kernel_poli, cmap=cm.jet, norm=LogNorm(vmin=1E-7,vmax=1E7))
# # plt.plot(np.log(kernel_poliangular[:,0]))
# plt.show()

print(rAxis,thetaAxis)

kernel_poli = np.concatenate( (kernel_poli,kernel_poli[:,::-1]) , axis=1 )

ker.plotPolarContour_full_360(kernel_poli, rAxis, thetaAxis, 4)


