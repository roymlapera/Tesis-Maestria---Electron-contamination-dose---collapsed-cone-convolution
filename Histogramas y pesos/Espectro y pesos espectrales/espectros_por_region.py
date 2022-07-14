import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'
# directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def func_exp(x,A,B):
    # aqui las x son pv y devuelve dosis
    return A*np.exp(-B*x)

def PasaDe100A24Bines(energias_espectro_red,pesos_espectro_energ_red_norm):
	energ_24 = np.arange(6.0/24.0, 6.0 + 6.0/24.0 , 6.0/24.0)
	pesos_24 = np.zeros(24)

	pesos_24[0] = pesos_espectro_energ_red_norm[energias_espectro_red < energ_24[0]].sum()

	for i in range(1,24):
		pesos_24[i] = pesos_espectro_energ_red_norm[(energias_espectro_red >= energ_24[i-1]) & (energias_espectro_red < energ_24[i])].sum()

	norma_24 = pesos_24.sum()*6.0/24.0
	pesos_24_norm = pesos_24/norma_24

	return energ_24, pesos_24_norm


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Genero los pesos angulares

path_EF = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Creacion del EF para el kernel de comparacion/'
# path_EF = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion del EF para el kernel de comparacion/'

path_energia = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Histogramas y pesos/Espectro y pesos espectrales/'
# path_energia = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Histogramas y pesos/Espectro y pesos espectrales/'

filename_F_reg_central      = 'energia_vs_fluencia_reg_central.txt'

bins = 72
maxEnergy = 6.0

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Guardo los datos del EF en un archivo txt

# tamano_pack_particulas = -1

# xinf = -5 #cm
# xsup = 5 #cm
# yinf = -5 #cm
# ysup = 5 #cm

# nro_EF = 0
# energias, weights = ef.DevuelveEnergiasDeArchivo(path_EF, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup,nro_EF)
# ef.Guarda_data(path_energia + filename_F_reg_central,energias,weights)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

energias, weights       = ef.Abre_energias(path_energia + filename_F_reg_central)

#la energia del EF es la energia total, le resto la energia en reposo

Area_EF = 30*30 #cm2
weights /= Area_EF

energias -= 0.511

histograma, bin_edges, patches = plt.hist(energias, bins=bins, range=[0,6] , histtype='step', normed=False, color='b', weights=weights)

plt.figure()
plt.bar(np.arange(6/72/2,6,6/72), histograma, 6/72, color = 'None', edgecolor = 'b')

plt.ylabel('dΦ/dE normalizado [1/MeV]', fontsize=14)
plt.xlabel('E cinética [MeV]', fontsize=14)
plt.grid(True)
plt.show()

# np.savetxt(path_energia + 'pesos_fluencia_72bines.txt', histograma)

# ------------------------------------------------------------------------------------------------------------------------

# f = 'pesos_fluencia_72bines.txt'
# pesos_espectro = np.loadtxt(path_energia + f, dtype='double', comments="#")

# delta_energia = maxEnergy/bins
# energias_espectro =  np.arange(delta_energia/2.0, maxEnergy, delta_energia)

# pesos_espectro /= pesos_espectro.sum()*delta_energia

# # ------------------------------------------------------------------------------------------------------------------------

# energias_espectro = energias_espectro-0.511

# energias_espectro_red = energias_espectro[energias_espectro >= delta_energia/2.0]
# pesos_espectro_red = pesos_espectro[energias_espectro >= delta_energia/2.0]

# #Fiteo segun sugiere Sikora en Sikora et al (2007)

# popt, pcov = curve_fit(func_exp, energias_espectro_red[3:], pesos_espectro_red[3:])

# pE = func_exp(energias_espectro_red, popt[0], popt[1])

# # --------------------------------------------------------------------------------------

# # grafico el espectro y la funcion de ajuste

# fig01 = plt.figure(figsize = (8,8))
# plt.bar(energias_espectro_red-delta_energia/2.0, pesos_espectro_red, delta_energia, color = 'None', edgecolor = 'b')
# # plt.plot(energias_espectro_red, pE)
# # plt.plot(energias_espectro_red[3], pesos_espectro_red[3], '*b')
# plt.xlim([0.0,6.0])
# # plt.ylim([0.0,0.6])

# ax = plt.gca()
# # Set axes ticks and labels
# plt.setp(ax.get_xticklabels(), fontsize=16)
# plt.setp(ax.get_yticklabels(), fontsize=16)
# ax.set_xlabel('Energia cinetica [MeV]', fontsize = 16)
# # ax.set_ylabel('Dosis [Gy/part]', fontsize = 16)
# ax.xaxis.labelpad = 10
# ax.yaxis.labelpad = 10

# plt.tight_layout() 
# plt.grid(True)
# plt.show()
# fig01.savefig('espectro_part_ajust.png', dpi = 200, bbox_inches='tight')

# --------------------------------------------------------------------------------------------

# # # completo los valores faltantes con la funcion de ajuste 

# pesos_espectro_red[:3] = func_exp(energias_espectro_red[:3], popt[0], popt[1])

# # ---------------------------------------------------------------------------------------------

# fig02 = plt.figure(figsize = (8,8))
# plt.bar(energias_espectro_red-delta_energia/2.0, pesos_espectro_red, delta_energia, color = 'None', edgecolor = 'r')

# # plt.plot(energias_espectro_red, pE)

# # plt.plot(energias_espectro_red[3], pesos_espectro_red[3], '*b')

# plt.xlim([0.0,6.0])
# # plt.ylim([0.0,0.6])
# ax = plt.gca()

# # Set axes ticks and labels
# plt.setp(ax.get_xticklabels(), fontsize=16)
# plt.setp(ax.get_yticklabels(), fontsize=16)

# ax.set_xlabel('Energia cinetica [MeV]', fontsize = 16)
# # ax.set_ylabel('Dosis [Gy/part]', fontsize = 16)

# ax.xaxis.labelpad = 10
# ax.yaxis.labelpad = 10

# #ax.set_ylim(rango)

# # legend = ax.legend(prop={'size': 14}, loc=(0.1,0.65))

# plt.tight_layout() 

# plt.show()

# fig02.savefig('espectro_part_posta.png', dpi = 200, bbox_inches='tight')

# # ---------------------------------------------------------------------------------------------

# # encuentro la fluencia energetica

# pesos_espectro_energ_red = pesos_espectro_red*energias_espectro_red 

# pesos_espectro_energ_red_norm = pesos_espectro_energ_red/(pesos_espectro_energ_red.sum()*delta_energia)    # Normalizo la fluencia energetica

# # --------------------------------------------------------------------------------------------

# fig03 = plt.figure(figsize = (8,8))
# plt.bar(energias_espectro_red-delta_energia/2.0, pesos_espectro_energ_red_norm, delta_energia, color = 'None', edgecolor = 'r')
# plt.xlim([0.0,6.0])
# # plt.ylim([0.0,0.4])
# ax = plt.gca()
# # Set axes ticks and labels
# plt.setp(ax.get_xticklabels(), fontsize=16)
# plt.setp(ax.get_yticklabels(), fontsize=16)
# ax.set_xlabel('Energia cinetica [MeV]', fontsize = 16)
# ax.set_ylabel('Fluencia energetica [1/MeV]', fontsize = 16)
# ax.xaxis.labelpad = 10
# ax.yaxis.labelpad = 10
# #ax.set_ylim(rango)
# # legend = ax.legend(prop={'size': 14}, loc=(0.1,0.65))
# plt.tight_layout() 
# plt.show()
# fig03.savefig('espectro_part_posta_2.png', dpi = 200, bbox_inches='tight')

# # --------------------------------------------------------------------------------------------

# # Hago rebinning de 100 a 8 bines de energia para generar los pesos de los 8 kernels monoenergeticos

# energ_24, pesos_24_norm = PasaDe100A24Bines(energias_espectro_red,pesos_espectro_energ_red_norm)

# fig04 = plt.figure(figsize = (8,8))
# plt.bar(energ_24 - 6.0/24.0/2, pesos_24_norm, 6.0/24.0, color = 'None', edgecolor = 'b')
# plt.bar(energias_espectro_red-delta_energia/2.0, pesos_espectro_energ_red_norm, delta_energia, color = 'None', edgecolor = 'r')

# plt.xlim([0.0,6.0])
# # plt.ylim([0.0,0.4])
# ax = plt.gca()
# # Set axes ticks and labels
# plt.setp(ax.get_xticklabels(), fontsize=16)
# plt.setp(ax.get_yticklabels(), fontsize=16)
# ax.set_xlabel('Energia cinetica [MeV]', fontsize = 16)
# ax.set_ylabel('Fluencia energetica [MeV]', fontsize = 16)
# ax.xaxis.labelpad = 10
# ax.yaxis.labelpad = 10
# #ax.set_ylim(rango)
# # legend = ax.legend(prop={'size': 14}, loc=(0.1,0.65))
# plt.tight_layout() 
# plt.show()
# fig04.savefig('espectro_energ_24.png', dpi = 200, bbox_inches='tight')

# # --------------------------------------------------------------------------------------------

# np.savetxt('pesos_espectro_fluencia_energ_24bines.txt', pesos_24_norm)

# # --------------------------------------------------------------------------------------------
