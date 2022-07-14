import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import convolve
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.special


directorio_scripts = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Scrips Python/'
# directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_espacio_fases as ef
import funciones_con_histogramas as histo

def DevuelveDosisDeKernelDesplazado(m,n,kernel,dosis, FE):
	# v_range1 = slice(max(0, m), max(min(m + kernel.shape[0], dosis.shape[0]), 0))
	# h_range1 = slice(max(0, n), max(min(n + kernel.shape[1], dosis.shape[1]), 0))
	# z_range1 = slice(max(0, 0), max(min(0 + kernel.shape[2], dosis.shape[2]), 0))

	# v_range2 = slice(max(0, -m), min(-m + dosis.shape[0], kernel.shape[0]))
	# h_range2 = slice(max(0, -n), min(-n + dosis.shape[1], kernel.shape[1]))
	# z_range2 = slice(max(0, -0), min(-0 + dosis.shape[2], kernel.shape[2]))

	# print(dosis[v_range1, h_range1,z_range1].shape,kernel[v_range2, h_range2,z_range2].shape)

	# pad = ( ( m-kernel.shape[0]//2 , dosis.shape[0]-(m-kernel.shape[0]//2) ) , ( n-kernel.shape[1]//2 , dosis.shape[1]-(n-kernel.shape[1]//2) ) , (0,0) )

	# kernel_padded = np.pad(kernel,pad,mode='constant')

	kernel_grande = np.zeros(dosis.shape, dtype=float)

	if(m<0):
		kernel = kernel[abs(m):,:,:]
		m += abs(m)

	if(m>dosis.shape[0]-kernel.shape[0]):
		kernel = kernel[:abs(m-dosis.shape[0]),:,:]

	if(n<0):
		kernel = kernel[:,abs(n):,:]
		n += abs(n)

	if(n>dosis.shape[1]-kernel.shape[1]):
		kernel = kernel[:,:abs(n-dosis.shape[1]),:]

	kernel_grande[ m: m + kernel.shape[0] , n : n + kernel.shape[1] ] = kernel

	dosis += kernel_grande * FE

	return dosis

def FE(x,y,l,sigma,alpha):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = a
    d = b

    sigma_x = sigma
    sigma_y = sigma

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return cte*alpha*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )

# --------------------------------------------------------------------------------------------------------------------------------------------------

kernel_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Convolucion pencil beam/'
# kernel_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Convolucion pencil beam/'

kernel_name = "kernel_MC_1000M_01cm.npz"

archivo = np.load(kernel_path + kernel_name)

# kernel = archivo['kernel']
kernel = archivo['dosis']
print(kernel.shape)
l = int(np.cbrt(kernel.shape[0]))
kernel = kernel.reshape([l,l,l])
kernel = np.transpose(kernel,(2,1,0))
print(kernel.shape)

limite_fantoma = 10.5
limite_kernel = 5.05

# print(kernel.sum()*(limite_kernel/(kernel.shape[0]/2))**3)

bins  = int((kernel.shape[0]-1)*limite_fantoma//limite_kernel + 1)
print(bins)

# ------------------------------------------------------ RESULTADOS DE AJUSTE DE LA FENERGETICA ------------------------------------------------------------------

sigma = 7.439

alpha = 0.556

l = 10

C = 8.72E-5

# ------------------------------------------------------ CONVOLUCION ------------------------------------------------------------------

kernel = kernel[:,:,0:kernel.shape[2]//2]   #me quedo con la parte inferior del kernel para computar menos voxeles

# kernel = kernel[:,:,::-1]

dosis = np.zeros((bins,bins,kernel.shape[2]), dtype=float)

print(kernel.shape)

# ---------------------------------------------------------------------------------------------------------------------------------------

# Pruebo tirar un pencil beam con peso 1 

#m puede ir [-kernel.shape[0]//2,dosis.shape[0]-kernel.shape[0]//2]

#n puede ir [-kernel.shape[1]//2,dosis.shape[1]-kernel.shape[1]//2]

# m = 241-80
# n = 241-80

# # m  y n son las posiciones del voxel [0,0], no del centro del kernel

# dosis = DevuelveDosisDeKernelDesplazado(m,n,kernel,dosis,1)#Fenergetica[m+(bins-1)//2][n+(bins-1)//2])
# dosis = DevuelveDosisDeKernelDesplazado(0,0,kernel,dosis,1)#Fenergetica[m+(bins-1)//2][n+(bins-1)//2])

# # plt.figure()
# # for i in range(0,dosis.shape[0],10):
# # 	plt.imshow(np.squeeze(dosis[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-4, vmax=10**5))
# # 	print(i)
# # 	plt.show()

# plt.figure()
# plt.imshow(dosis[:,bins//2+1,:], cmap=cm.hot, norm=LogNorm(vmin=10**-4, vmax=10**5))
# plt.figure()
# plt.imshow(dosis[:,:,1], cmap=cm.hot, norm=LogNorm(vmin=10**-4, vmax=10**5))
# plt.show()
# ---------------------------------------------------------------------------------------------------------------------------------------

#La convolucion

# x = np.linspace(-limite_fantoma,limite_fantoma,bins)
# y = np.linspace(-limite_fantoma,limite_fantoma,bins)

# X, Y = np.meshgrid(x,y)

# for m in range(-kernel.shape[0]//2, dosis.shape[0]-kernel.shape[0]//2):
# 	for n in range(-kernel.shape[1]//2, dosis.shape[1]-kernel.shape[1]//2):
# 		m_centro_kernel = m + kernel.shape[0]//2
# 		n_centro_kernel = n + kernel.shape[1]//2
# 		# print(m_centro_kernel,n_centro_kernel)
# 		dosis = DevuelveDosisDeKernelDesplazado(m,n,kernel,dosis,FE(X[m_centro_kernel,n_centro_kernel],Y[m_centro_kernel,n_centro_kernel],l,sigma,alpha)+C)
# 		sys.stdout.write("\r{:.2f}".format((n_centro_kernel+1+(m_centro_kernel+1)*dosis.shape[0])/(dosis.shape[0]*dosis.shape[1])*100) + ' %')
# 		sys.stdout.flush() 

# np.savez("dosis_conv_pencilbeam_kernel_01mm_prueba_p_PDD.npz", dosis = dosis)

# archivo_dosis = np.load("dosis_conv_pencilbeam_kernel_01mm_prueba_p_PDD.npz")
# dosis = archivo_dosis['dosis']

# plt.figure()
# plt.plot(dosis[dosis.shape[0]//2,dosis.shape[1]//2,:],'b.')

# plt.figure()
# # plt.imshow(dosis[:,dosis.shape[0]//2,:], cmap=cm.hot)
# plt.figure()
# plt.imshow(kernel[:,kernel.shape[0]//2,:], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))

# plt.show()

# ---------------------------------------------------------- APERTURA DE ARCHIVO DE DOSIS GUARDADO: PDD ---------------------------------------------------------

# archivo_dosis = np.load("dosis_conv_pencilbeam_reg_central_241_vox.npz")

# dosis = archivo_dosis['dosis']

# # plt.figure()
# # plt.imshow(dosis_MC[:,1,:], cmap=cm.hot, vmin=0, vmax=10**-16)
# # plt.show()

# eje_z_CS =  np.arange(0,dosis.shape[2])  * limite_fantoma*2/(bins-1) + limite_fantoma*2/(bins-1)/2

# plt.figure()
# plt.plot(eje_z_CS,dosis[bins//2+1,bins//2+1,:]/dosis[bins//2+1,bins//2+1,3],'b.-', linewidth=0.5,  label='PDD_CS')
# plt.legend()
# plt.show()

# # ---------------------------------------------------------------------------------------------- ---------------------------------------------------------

# archivo_dosis41 = np.load("dosis_conv_pencilbeam_121_vox.npz")
# dosis41 = archivo_dosis41['dosis']
# eje_z_CS41 =  np.arange(0,dosis41.shape[2])  * limite_fantoma*2/(dosis41.shape[1]-1) + limite_fantoma*2/(dosis41.shape[1]-1)/2

# archivo_dosis81 = np.load("dosis_conv_pencilbeam_241_vox.npz")
# dosis81 = archivo_dosis81['dosis']
# eje_z_CS81 =  np.arange(0,dosis81.shape[2])  * limite_fantoma*2/(dosis81.shape[1]-1) + limite_fantoma*2/(dosis81.shape[1]-1)/2

# archivo_dosis161 = np.load("dosis_conv_pencilbeam_481_vox.npz")
# dosis161 = archivo_dosis161['dosis']
# eje_z_CS161 =  np.arange(0,dosis161.shape[2])  * limite_fantoma*2/(dosis161.shape[1]-1) + limite_fantoma*2/(dosis161.shape[1]-1)/2

# print(eje_z_CS41)
# print(eje_z_CS81)
# print(eje_z_CS161)

# # plt.figure()
# plt.plot(eje_z_CS41,dosis41[dosis41.shape[1]//2+1,dosis41.shape[1]//2+1,:]/dosis41[dosis41.shape[1]//2+1,dosis41.shape[1]//2+1,1],'b.-', linewidth=0.5, label='PDD_CS41')
# plt.plot(eje_z_CS81,dosis81[dosis81.shape[1]//2+1,dosis81.shape[1]//2+1,:]/dosis81[dosis81.shape[1]//2+1,dosis81.shape[1]//2+1,3],'g*-', linewidth=0.5, label='PDD_CS81')
# plt.plot(eje_z_CS161,dosis161[dosis161.shape[1]//2+1,dosis161.shape[1]//2+1,:]/dosis161[dosis161.shape[1]//2+1,dosis161.shape[1]//2+1,6],'m+-', linewidth=0.5, label='PDD_CS161')

# plt.legend()
# plt.show()

# # # ---------------------------------------------------------------------------------------------- ---------------------------------------------------------

# archivo_dosis_MC = np.load("dosis_MC_1000M.npz")

# dosis_MC = archivo_dosis_MC['dosis']
# dosis_MC = dosis_MC.reshape((17,3,101))


# eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10.5,5)+2.5 )).ravel()

# print(eje_z_MC[4])

# plt.plot(eje_z_MC[:-1], dosis_MC[:,1,51]/dosis_MC[4,1,51],'r^-', linewidth=0.5, label='PDD_MC')
# plt.legend()
# plt.show()



# # ---------------------------------------------------------- PERFIL DE DOSIS A 0.5 cm ---------------------------------------------------------

# # eje_x_CS =  np.arange(0,60)  * limite_fantoma*2/(bins-1) + limite_fantoma*2/(bins-1)/2
# # eje_x_CS = np.asarray( list(-1*eje_x_CS[::-1]) + [0] + list(eje_x_CS) )

# # plt.figure()
# # plt.plot(eje_x_CS, dosis[:,bins//2+1,1]/dosis[:,bins//2+1,1].max(),'b.', label='Perfil_dosis_0.75cm_CS')
# # plt.legend()

# # eje_x_MC = np.hstack(( np.arange(0.5,20,0.5),np.arange(2,5,0.5)+0.25)).ravel()

# # plt.plot(eje_z_MC[:-1], dosis_MC[:,1,51]/dosis_MC[:,1,51].max(),'r.', label='Perfil_dosis_0.75cm_MC')
# # plt.legend()
# # plt.show()

# # plt.show()