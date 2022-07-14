import sys
import numpy as np
import matplotlib.pyplot as plt

# directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'
directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels_v2 as ker2

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
        volume[:,k] = 2*np.pi * (np.cos(ThetaAxis[k]*np.pi/180)-np.cos(thetaAxis[k]*np.pi/180))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    Frac_energ_esf = kernel*volume

    print('La integral del kernel en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume

# ------------------------------------------------------------------------------------------------

i = 8

# filename = 'espectro' + str(i)+'.txt'

filename = 'pesos_FE_24bines.txt'

pesos_espectro = np.genfromtxt(filename)

# pesos_espectro /= pesos_espectro.sum()*6/21

# pesos_espectro[0] = 0
# pesos_espectro[1] = 0
# pesos_espectro[2] = 0

print(pesos_espectro.sum()*6/24) 

energia_max = 6.0
bins = 24
delta_energia = energia_max/bins

# energias_espectro =  np.arange(3*6/24+delta_energia/2.0, 6.0, delta_energia)
energias_espectro =  np.arange(delta_energia/2.0, 6.0, delta_energia)

# plt.figure()
# plt.bar(energias_espectro,pesos_espectro, fill='False')
# plt.show()

# # ------------------------------------------------------------------------------------------------

kernel_pol_file = np.load('fullWater_50M_0375MeV.npz')
kernel_pol_0375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_1125MeV.npz')
kernel_pol_1125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_1875MeV.npz')
kernel_pol_1875MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_2625MeV.npz')
kernel_pol_2625MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_3375MeV.npz')
kernel_pol_3375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_4125MeV.npz')
kernel_pol_4125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_4875MeV.npz')
kernel_pol_4875MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_5625MeV.npz')
kernel_pol_5625MeV = kernel_pol_file['kernel']



kernel_pol_file = np.load('fullWater_50M_0125MeV.npz')
kernel_pol_0125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_0625MeV.npz')
kernel_pol_0625MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_0875MeV.npz')
kernel_pol_0875MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_1375MeV.npz')
kernel_pol_1375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_1625MeV.npz')
kernel_pol_1625MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_2125MeV.npz')
kernel_pol_2125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_2375MeV.npz')
kernel_pol_2375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_2875MeV.npz')
kernel_pol_2875MeV = kernel_pol_file['kernel']



kernel_pol_file = np.load('fullWater_50M_3125MeV.npz')
kernel_pol_3125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_3625MeV.npz')
kernel_pol_3625MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_3875MeV.npz')
kernel_pol_3875MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_4375MeV.npz')
kernel_pol_4375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_4625MeV.npz')
kernel_pol_4625MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_5125MeV.npz')
kernel_pol_5125MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_5375MeV.npz')
kernel_pol_5375MeV = kernel_pol_file['kernel']

kernel_pol_file = np.load('fullWater_50M_5875MeV.npz')
kernel_pol_5875MeV = kernel_pol_file['kernel']

# # ------------------------------------------------------------------------------------------------

kernel_polienergetico1   = kernel_pol_0125MeV*pesos_espectro[0] + kernel_pol_0375MeV*pesos_espectro[1] + kernel_pol_0625MeV*pesos_espectro[2] + kernel_pol_0875MeV*pesos_espectro[3] + kernel_pol_1125MeV*pesos_espectro[4] + kernel_pol_1375MeV*pesos_espectro[5] + kernel_pol_1625MeV*pesos_espectro[6] + kernel_pol_1875MeV*pesos_espectro[7]
kernel_polienergetico2   = kernel_pol_2125MeV*pesos_espectro[8] + kernel_pol_2375MeV*pesos_espectro[9] + kernel_pol_2625MeV*pesos_espectro[10] + kernel_pol_2875MeV*pesos_espectro[11] + kernel_pol_3125MeV*pesos_espectro[12] + kernel_pol_3375MeV*pesos_espectro[13] + kernel_pol_3625MeV*pesos_espectro[14] + kernel_pol_3875MeV*pesos_espectro[15]
kernel_polienergetico3   = kernel_pol_4125MeV*pesos_espectro[16] + kernel_pol_4375MeV*pesos_espectro[17] + kernel_pol_4625MeV*pesos_espectro[18] + kernel_pol_4875MeV*pesos_espectro[19] + kernel_pol_5125MeV*pesos_espectro[20] + kernel_pol_5375MeV*pesos_espectro[21] + kernel_pol_5625MeV*pesos_espectro[22] + kernel_pol_5875MeV*pesos_espectro[23]

kernel_polienergetico = kernel_polienergetico1 + kernel_polienergetico2 + kernel_polienergetico3

kernel_polienergetico = kernel_polienergetico*delta_energia

# # ------------------------------------------------------------------------------------------------

kernel_pol_file = np.load('fullWater_50M_0375MeV.npz')
rAxis = kernel_pol_file['rAxis']
thetaAxis = kernel_pol_file['thetaAxis']

Integral_vol_del_kernel(kernel_polienergetico,rAxis,thetaAxis)

# np.savez('kernel_polienergetico_roy_final'+str(i)+'.npz', kernel = kernel_polienergetico, rAxis = rAxis, thetaAxis = thetaAxis)
# np.savez('kernel_polienergetico_roy_final_final.npz', kernel = kernel_polienergetico, rAxis = rAxis, thetaAxis = thetaAxis)

# ------------------------------------------------------------------------------------------------

kernel_pol_file = np.load('kernel_polienergetico_roy_final_final.npz')

rAxis = kernel_pol_file['rAxis']
thetaAxis = kernel_pol_file['thetaAxis']
kernel_polienergetico_ff = kernel_pol_file['kernel']



# ker2.plotPolarContour_full(kernel_polienergetico, rAxis, thetaAxis, 4)

# plt.figure()
# plt.plot(rAxis, np.log10(kernel_polienergetico[:,0]))
# plt.show()

# ------------------------------------------------------------------------------------------------

# kernel_pol_file = np.load('kernel_polienergetico_roy_24.npz')

# rAxis = kernel_pol_file['rAxis']
# thetaAxis = kernel_pol_file['thetaAxis']
# kernel_polienergetico = kernel_pol_file['kernel']

# Integral_vol_del_kernel(kernel_polienergetico,rAxis,thetaAxis)

# ker2.plotPolarContour_full(kernel_polienergetico, rAxis, thetaAxis, 4)

deltaR = np.zeros(len(rAxis))
deltaR[0] = rAxis[0]

k=0
while k < len(rAxis)-1:
    deltaR[k+1] = rAxis[k+1]-rAxis[k]
    k+=1
r = [0] + rAxis
r = r - deltaR/2

plt.figure()
for i in range(9):
    plt.semilogy(r[r<0.2], kernel_polienergetico_ff[:,10*i][r<0.2], 'x-', label='Corte a ' + str(10*i) + '°')

plt.xlabel('Radio [cm]')
plt.ylabel('Fracción de energía depositada [1/cm$\mathregular{^3}$]')
plt.grid(True)
plt.legend()
plt.show()

