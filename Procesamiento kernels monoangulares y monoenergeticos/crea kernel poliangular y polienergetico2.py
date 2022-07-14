import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from scipy.integrate import quad



directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'
# directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels_v2 as ker2
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
        volume[:,k] = 2*np.pi * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    # print(volume[0,0])

    Frac_energ = kernel*volume

    # print('La integral del kernel polienergetico en la esfera de radio 10 da: ' + str(Frac_energ.sum()))

    return Frac_energ

def Integral_vol_del_kernel(kernel,rAxis,thetaAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion imprime por pantalla el resultado de la integral
    y devuelve la matriz de volumenes con ejes segun rAxis thetaAxis y phiAxis.'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < (len(rAxis)-1):
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel.shape)

    for k in range(len(thetaAxis)):
        volume[:,k] = 2*np.pi/360 * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    Frac_energ_esf = kernel*volume

    print('La integral del kernel en la esfera de radio 10 da: ' + str(360*Frac_energ_esf.sum()))

    return volume

def Calcula_funcion_peso(file_path, archivo_pesos, delta_alpha, delta_beta,alpha_max):

    angulos = np.arange(delta_alpha/2,np.pi/2,delta_alpha)

    pesos_angulares = histo.Abre_archivo_pesos_angulares(file_path + archivo_pesos)

    angulos = angulos[angulos<alpha_max]

    pesos_angulares = pesos_angulares[:angulos.shape[0]]

    # normalizo los pesos para que esten normalizada la integral de p en alpha y beta

    pesos_angulares /= pesos_angulares.sum()*delta_alpha*2*np.pi

    # print(pesos_angulares.sum()*delta_alpha*2*np.pi)  #ahora esta normalizado en 2D, la integral que queremos ahcer para hacer el kernel poliangular es 2D 

    p = lambda th: interp1d(angulos,pesos_angulares, kind='cubic', bounds_error=False, fill_value=(pesos_angulares[0],0))(th)

    plt.figure()
    plt.plot(angulos*180/np.pi,pesos_angulares,'b.')
    plt.plot(np.linspace(0,90,1000),p(np.linspace(0,90,1000)*np.pi/180),'r')
    plt.show()

    print(pesos_angulares.sum()*delta_alpha*2*np.pi , 2*np.pi*quad(p, 0, alpha_max)[0])

    return p

#-------------------------------------------------------------------------------- MAIN ------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Procesamiento kernels monoangulares y monoenergeticos/'


archivo_pesos = 'pesos_angulares_FE_reg_central_45bines.txt'

NC_th = 180        #Nro de conos del kernel en direccion theta
NC_phi = 360        #Nro de conos del kernel en direccion phi
NR = 140#210        #Nro de radios del kernel

delta_alpha = np.pi/90  #para que nc_rotacion sea igual a 4 grados

nc_rotacion = int(delta_alpha / (np.pi/180))

archivo = np.load("kernel_polienergetico_roy.npz")

kernel_poli = archivo['kernel']
rAxis     = archivo['rAxis'] 
thetaAxis = archivo['thetaAxis'] * np.pi/180
phiAxis   = np.linspace(0,360,NC_phi) * np.pi/180

frac_energ = CreaFracEnerg(kernel_poli,rAxis,thetaAxis) 

volume = Integral_vol_del_kernel(kernel_poli,rAxis,thetaAxis)

kernel_poli = frac_energ/360/volume   # ahora el kernel corresponde a a la densidad de energia de voxeles esfericos y no de anillos esfericos

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Importacion del peso angular y calculo de la funcion peso

delta_beta = np.pi/9  #10 grados

alpha_max = 40 * np.pi/180  #40 grados

p = Calcula_funcion_peso(file_path, archivo_pesos, delta_alpha, delta_beta, alpha_max)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

kernel_poliangular = ker2.CreaKernelPoliangular(rAxis[:140], kernel_poli, p, NC_th, NR, delta_alpha, delta_beta,alpha_max)

# np.savez("kernel_poliangular_prueba.npz", kernel = kernel_poliangular, rAxis = rAxis, thetaAxis = thetaAxis)

archivo = np.load("kernel_poliangular_roy.npz")

kernel_poliangular = archivo['kernel']

volume = Integral_vol_del_kernel(kernel_poliangular,rAxis[:140],thetaAxis)

kernel_poliangular /= (kernel_poliangular*volume*360).sum()

print((kernel_poliangular*volume*360).sum())

# volume = Integral_vol_del_kernel(kernel_poliangular,rAxis,thetaAxis)

# plt.figure()
# plt.imshow(kernel_poliangular, cmap=cm.jet, norm=LogNorm(vmin=1E-7,vmax=1E7))
# plt.figure()
# # plt.imshow(kernel_poli, cmap=cm.jet, norm=LogNorm(vmin=1E-7,vmax=1E7))
# plt.plot(rAxis,np.log10(kernel_poli[:,60]), 'b.')
# # plt.plot(rAxis,np.log10(kernel_poliangular[:,60]), 'g.')
# plt.show()


kernel_poliangular = np.concatenate( (kernel_poliangular,kernel_poliangular[:,::-1]) , axis=1 )
# kernel_poli = np.concatenate( (kernel_poli,kernel_poli[:,::-1]) , axis=1 )

ker2.plotPolarContour_full_360(kernel_poliangular, rAxis[:140], thetaAxis*180/np.pi, 4)