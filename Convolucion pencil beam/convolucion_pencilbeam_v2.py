import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import convolve
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import scipy.special

def FE(x,y,l,sigma,alpha,C):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = a
    d = b

    sigma_x = sigma
    sigma_y = sigma

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return cte*alpha*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )+C

def Calcula_Dosis_en_vox(FE,l,sigma,alpha,C,kernel,tam_vox_cart,i,j,k):
	dosis = 0
	for m in range(kernel.shape[0]):
		for n in range(kernel.shape[1]):
			if(k<kernel.shape[2]):
				dosis += kernel[m,n,k] * FE(tam_vox_cart*(i-m),tam_vox_cart*(j-n),l,sigma,alpha,C)
	return dosis*tam_vox_cart*tam_vox_cart

def Calcula_Dosis(dosis,FE,l,sigma,alpha,C,kernel,tam_vox_cart):
	for i in range(dosis.shape[0]):
		for j in range(dosis.shape[1]):
			for k in range(dosis.shape[2]):
				dosis[i][j][k] = Calcula_Dosis_en_vox(FE,l,sigma,alpha,C,kernel,tam_vox_cart,i,j,k)
	return dosis

# --------------------------------------------------------------------------------------------------------------------------------------------------

# kernel_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Convolucion pencil beam/'
kernel_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Convolucion pencil beam/'

# kernel_name = "Varian-20-IAEA_electrones_1000M_kernel_vactowater.npz"
# kernel_name = "kernel_cartesiano_roy_final_final_c_interp.npz"
kernel_name = "kernel_cartesiano_roy_final_final_s_interp.npz"

archivo = np.load(kernel_path + kernel_name)

kernel = archivo['kernel']

# kernel = archivo['dosis']
# print(kernel.shape)

# print(kernel.sum()*0.1*0.1*0.1)
# kernel /= kernel.sum()*0.1*0.1*0.1
# l = int(np.cbrt(kernel.shape[0]))
# kernel = kernel.reshape([50,101,101])
# kernel = np.transpose(kernel,(2,1,0))
# # kernel = kernel[:,:,kernel.shape[2]//2:]
# print(kernel.shape)

# kernel = kernel[20:81,20:81,:31]

# print(kernel.min(),kernel.max())

# plt.figure()
# plt.imshow(kernel[:,30,:], cmap=cm.hot, norm=LogNorm(vmin=10**-8, vmax=10**2), extent=[0,3.1,-3.05,3.05])
# plt.xlabel('Posición Z [cm]')
# plt.ylabel('Posición X [cm]')
# plt.colorbar()
# plt.show()


# ------------------------------------------------------- CONVOLUCION ------------------------------------------------------- 

tam_vox_cart = 0.1#cm

xmin = -0.05
xmax = 0.05
ymin = -0.05
ymax = 0.05
zmin = 0
zmax = 3

nro_vox_gra_x = round(abs(xmax-xmin)/tam_vox_cart)
nro_vox_gra_y = round(abs(ymax-ymin)/tam_vox_cart)
nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)

# ------------------------------ RESULTADOS DE AJUSTE DE LA FENERGETICA  -----------------------------------------------

sigma = 7.439
alpha = 0.556
l = 10
C = 8.72E-5

# ------------------------------------------------------- CONVOLUCION ------------------------------------------------------- 

X = np.linspace(xmin+tam_vox_cart/2,xmax-tam_vox_cart/2,nro_vox_gra_x)
Y = np.linspace(ymin+tam_vox_cart/2,ymax-tam_vox_cart/2,nro_vox_gra_y)
Z = np.linspace(zmin+tam_vox_cart/2,zmax-tam_vox_cart/2,nro_vox_gra_z)

dosis = np.zeros((nro_vox_gra_x,nro_vox_gra_y,nro_vox_gra_z))

print(dosis.shape)

dosis = Calcula_Dosis(dosis,FE,l,sigma,alpha,C,kernel,tam_vox_cart)

# # ------------------------------------------------------- GUARDO PDD ------------------------------------------------------- 

# np.savez("dosis_conv_cart__kernel_s_interp.npz", dosis = dosis)

# plt.figure()
# plt.imshow(dosis[:,dosis.shape[1]//2+1,:], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**-5))
# plt.figure()
# plt.imshow(dosis[:,:,1], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**-5))
# plt.show()


# ---------------------------------------------------------- APERTURA DE ARCHIVO DE DOSIS GUARDADO: PDD ---------------------------------------------------------

prof_norm = 0.5 #cm

# ----------------------------------------

archivo_dosis_CS = np.load("dosis_conv_cart_01_001_50M_fullwater_2125MeV.npz")
dosis_CS = dosis#archivo_dosis_CS['dosis']

print(dosis_CS.shape)

eje_z_CS = Z
dosis_CS_norm = np.interp(prof_norm,eje_z_CS,dosis_CS[dosis_CS.shape[0]//2,dosis_CS.shape[1]//2,:])

# ----------------------------------------

# archivo_dosis_CC = np.load("dosis_conv_conocolapsado_5_15_kernel_poliangular_reg_central.npz")
# dosis_CC = archivo_dosis_CC['dosis']
# print(dosis_CC.shape)

# eje_z_CC = Z
# dosis_CC_norm = np.interp(prof_norm,eje_z_CC,dosis_CC[dosis_CC.shape[0]//2,dosis_CC.shape[1]//2,:])

# # # ----------------------------------------

# PDD calculado en DOSXYZnrc de electrones contaminantes

archivo_dosis_econt = np.load("Varian-20-IAEA_electrones_1000M_paraPDD.npz")
dosis_econt = archivo_dosis_econt['dosis']
error_econt = archivo_dosis_econt['errs']

dosis_econt = dosis_econt.reshape((34,3,3))

error_econt = error_econt.reshape((34,3,3))

print(dosis_econt.shape)

dosis_econt = dosis_econt[:,dosis_econt.shape[1]//2,dosis_econt.shape[2]//2]
error_econt = error_econt[:,error_econt.shape[1]//2,error_econt.shape[2]//2]

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
dosis_econt_norm    = np.interp(prof_norm,eje_z,dosis_econt)

# # # ----------------------------------------

# PDD calculado en DOSXYZnrc de electrones de 2.125 MeV cuya incidencia es normal a la sup

archivo_dosis_monoAmonoE = np.load("Varian-20-IAEA_monoangular_y_monoenergetica_electrones_1000M_paraPDD.npz")
dosis_monoAmonoE = archivo_dosis_monoAmonoE['dosis']
error_monoAmonoE = archivo_dosis_monoAmonoE['errs']

dosis_monoAmonoE = dosis_monoAmonoE.reshape((34,11,11))

error_monoAmonoE = error_monoAmonoE.reshape((34,11,11))

dosis_monoAmonoE = dosis_monoAmonoE[:,dosis_monoAmonoE.shape[1]//2,dosis_monoAmonoE.shape[2]//2]
error_monoAmonoE = error_monoAmonoE[:,error_monoAmonoE.shape[1]//2,error_monoAmonoE.shape[2]//2]

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
dosis_monoAmonoE_norm    = np.interp(prof_norm,eje_z,dosis_monoAmonoE)

# ----------------------------------------



plt.figure()

plt.plot(eje_z_CS, 100*dosis_CS[dosis_CS.shape[0]//2,dosis_CS.shape[1]//2,:]/dosis_CS_norm,'gx-', linewidth=0.5, label='PDD_Conv cartesiana')

# plt.plot(eje_z_CC, 100*dosis_CC[dosis_CC.shape[0]//2,dosis_CC.shape[1]//2,:]/dosis_CC_norm,'m.-', linewidth=0.5, label='PDD_CC')

plt.errorbar(eje_z, 100*dosis_econt/dosis_econt_norm, yerr=100*error_econt*dosis_econt/dosis_econt_norm, linewidth=0.5, c='r', label='PDD_MC')

# print((abs(dosis_CS[dosis_CS.shape[0]//2,dosis_CS.shape[1]//2,:]/dosis_CS_norm-dosis_econt/dosis_econt_norm)/dosis_econt/dosis_econt_norm).max())

# plt.plot(eje_z, 100*dosis_monoAmonoE/dosis_monoAmonoE_norm,'r', linewidth=0.5, label='PDD_MC_monoAmonoE')

plt.ylabel('PDD [%]')
plt.xlabel('Profundidad [cm]')
plt.legend()
plt.grid(True)
plt.show()

# ---------------------------------------------------------- PDD COMPARA POLIANGULAR POLIENERGETICO ---------------------------------------------------------

# prof_norm = 1 #cm

# # ----------------------------------------

# archivo_dosis_poliA = np.load("dosis_conv_cart_kernel_poliangular_reg_central_cart_01_001_PDD.npz")
# dosis_poliA = archivo_dosis_poliA['dosis']

# eje_z_poliA = Z
# dosis_poliA_norm = np.interp(prof_norm,eje_z_poliA,dosis_poliA[dosis_poliA.shape[0]//2,dosis_poliA.shape[1]//2,:])

# # ----------------------------------------

# archivo_dosis_poliE = np.load("dosis_conv_cart_kernel_poliE_cart_01_001.npz")
# dosis_poliE = archivo_dosis_poliE['dosis']

# eje_z_poliE = Z
# dosis_poliE_norm = np.interp(prof_norm,eje_z_poliE,dosis_poliE[dosis_poliE.shape[0]//2,dosis_poliE.shape[1]//2,:])

# # ----------------------------------------

# archivo_dosis_MC = np.load("dosis_MC_1000M.npz")
# dosis_MC = archivo_dosis_MC['dosis']

# dosis_MC = dosis_MC.reshape((17,3,101))
# dosis_MC = np.transpose(dosis_MC,(2,1,0))

# eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10,5)+2.5 )).ravel()
# dosis_MC_norm = np.interp(prof_norm,eje_z_MC,dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:])

# # ----------------------------------------

# plt.figure()

# plt.plot(eje_z_poliA, 100*dosis_poliA[dosis_poliA.shape[0]//2,dosis_poliA.shape[1]//2,:]/dosis_poliA_norm,'b.-', linewidth=0.5, label='PDD_poliA')

# plt.plot(eje_z_poliE, 100*dosis_poliE[dosis_poliE.shape[0]//2,dosis_poliE.shape[1]//2,:]/dosis_poliE_norm,'m.-', linewidth=0.5, label='PDD_poliE')

# plt.plot(eje_z_MC, 100*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm,'r.-', linewidth=0.5, label='PDD_MC')

# plt.xlabel('Profundidad [cm]')
# plt.legend()
# plt.show()





