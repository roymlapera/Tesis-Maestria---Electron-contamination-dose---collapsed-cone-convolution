# ---------------------------------------------------------- PDD COMPARA POLIANGULAR POLIENERGETICO FE1 ---------------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math

# -----------------------------------------------------------------------------------------------------------------------------------------

prof_norm = 0.5 #cm

tam_vox_cart = 0.1#cm
zmin = 0
zmax = 4

nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)

Z = np.linspace(zmin+tam_vox_cart/2,zmax-tam_vox_cart/2,nro_vox_gra_z)

# ----------------------------------------

archivo_dosis_poliA_CS = np.load("dosis_conv_cart_01_001_kernel_poliangular_reg_central_FE1.npz")
dosis_poliA_CS = archivo_dosis_poliA_CS['dosis']

eje_z_poliA_CS = Z
dosis_poliA_CS_norm = np.interp(prof_norm,eje_z_poliA_CS,dosis_poliA_CS[dosis_poliA_CS.shape[0]//2,dosis_poliA_CS.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_poliE_CS = np.load("dosis_conv_cart_01_001_kernel_polienergetico_roy_FE1.npz")
dosis_poliE_CS = archivo_dosis_poliE_CS['dosis']

eje_z_poliE_CS = Z
dosis_poliE_CS_norm = np.interp(prof_norm,eje_z_poliE_CS,dosis_poliE_CS[dosis_poliE_CS.shape[0]//2,dosis_poliE_CS.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_poliA_CC = np.load("dosis_conv_conocolapsado_5_15_kernel_poliangular_reg_central_FE1.npz")
dosis_poliA_CC = archivo_dosis_poliA_CC['dosis']

eje_z_poliA_CC = Z
dosis_poliA_CC_norm = np.interp(prof_norm,eje_z_poliA_CC,dosis_poliA_CC[dosis_poliA_CC.shape[0]//2,dosis_poliA_CC.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_poliE_CC = np.load("dosis_conv_conocolapsado_5_15_kernel_polienergetico_roy_FE1.npz")
dosis_poliE_CC = archivo_dosis_poliE_CC['dosis']

eje_z_poliE_CC = Z
dosis_poliE_CC_norm = np.interp(prof_norm,eje_z_poliE_CC,dosis_poliE_CC[dosis_poliE_CC.shape[0]//2,dosis_poliE_CC.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_poliA_CC_EFinf = np.load("dosis_conv_conocolapsado_5_15_kernel_poliangular_reg_central_FE1_EFinf.npz")
dosis_poliA_CC_EFinf = archivo_dosis_poliA_CC_EFinf['dosis']

eje_z_poliA_CC_EFinf = Z
dosis_poliA_CC_EFinf_norm = np.interp(prof_norm,eje_z_poliA_CC_EFinf,dosis_poliA_CC_EFinf[dosis_poliA_CC_EFinf.shape[0]//2,dosis_poliA_CC_EFinf.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_poliE_CC_EFinf = np.load("dosis_conv_conocolapsado_5_15_kernel_polienergetico_roy_FE1.npz")
dosis_poliE_CC_EFinf = archivo_dosis_poliE_CC_EFinf['dosis']

eje_z_poliE_CC_EFinf = Z
dosis_poliE_CC_EFinf_norm = np.interp(prof_norm,eje_z_poliE_CC_EFinf,dosis_poliE_CC_EFinf[dosis_poliE_CC_EFinf.shape[0]//2,dosis_poliE_CC_EFinf.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_kernel_MC = np.load("dosis_conv_cart_01_001_kernel_MC_1000M_01.npz")
dosis_kernel_MC = archivo_dosis_kernel_MC['dosis']

eje_z_kernel_MC = Z
dosis_kernel_MC_norm = np.interp(prof_norm,eje_z_kernel_MC,dosis_kernel_MC[dosis_kernel_MC.shape[0]//2,dosis_kernel_MC.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_kernel_MC_vactowat = np.load("dosis_conv_cart_01_001_1000M_kernel_vactowater.npz")
dosis_kernel_MC_vactowat = archivo_dosis_kernel_MC_vactowat['dosis']

eje_z_kernel_MC_vactowat = Z
dosis_kernel_MC_vactowat_norm = np.interp(prof_norm,eje_z_kernel_MC_vactowat,dosis_kernel_MC_vactowat[dosis_kernel_MC_vactowat.shape[0]//2,dosis_kernel_MC_vactowat.shape[1]//2,:])

# ----------------------------------------

# archivo_dosis_MC = np.load("dosis_MC_1000M.npz")
# dosis_MC = archivo_dosis_MC['dosis']

# # dosis_MC = dosis_MC.reshape((17,3,101))
# # dosis_MC = np.transpose(dosis_MC,(2,1,0))

# eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10,5)+2.5 )).ravel()
# dosis_MC_norm = np.interp(prof_norm,eje_z_MC,dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:])

# ----------------------------------------

archivo_dosis_electrones = np.load("Varian-20-IAEA_electrones_1000M_paraPDD.npz")
dosis_electrones = archivo_dosis_electrones['dosis']
error_electrones = archivo_dosis_electrones['errs']

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
dosis_electrones_norm    = np.interp(prof_norm,eje_z,dosis_electrones)

# ----------------------------------------

archivo_dosis_electrones_monoangular = np.load("Varian-20-IAEA_monoangular_electrones_1000M_paraPDD.npz")
dosis_electrones_monoangular = archivo_dosis_electrones_monoangular['dosis']
error_electrones_monoangular = archivo_dosis_electrones_monoangular['errs']

dosis_electrones_monoangular = dosis_electrones_monoangular.reshape((34,11,11))
error_electrones_monoangular = error_electrones_monoangular.reshape((34,11,11))

dosis_electrones_monoangular = dosis_electrones_monoangular[:,dosis_electrones_monoangular.shape[1]//2,dosis_electrones_monoangular.shape[2]//2]
error_electrones_monoangular = error_electrones_monoangular[:,error_electrones_monoangular.shape[1]//2,error_electrones_monoangular.shape[2]//2]


eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
dosis_electrones_monoangular_norm    = np.interp(prof_norm,eje_z,dosis_electrones_monoangular)

# ----------------------------------------

# plt.figure()
# plt.imshow(dosis_electrones_monoangular.reshape((11,34,11)))
# plt.show()

# error_electrones_monoangular = error_electrones_monoangular.reshape((11,11,34))

# dosis_electrones_monoangular = dosis_electrones_monoangular[dosis_electrones_monoangular.shape[0]//2,:,dosis_electrones_monoangular.shape[2]//2]
# error_electrones_monoangular = error_electrones_monoangular[error_electrones_monoangular.shape[0]//2,error_electrones_monoangular.shape[1]//2,:]

# eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
# dosis_electrones_monoangular_norm    = np.interp(prof_norm,eje_z,dosis_electrones_monoangular)

# ----------------------------------------

plt.figure()

# plt.plot(eje_z_poliA_CS, 100*dosis_poliA_CS[dosis_poliA_CS.shape[0]//2,dosis_poliA_CS.shape[1]//2,:]/dosis_poliA_CS_norm,'b.-', linewidth=0.5, label='PDD_poliA_CS')

# plt.plot(eje_z_poliE_CS, 100*dosis_poliE_CS[dosis_poliE_CS.shape[0]//2,dosis_poliE_CS.shape[1]//2,:]/dosis_poliE_CS_norm,'m.-', linewidth=0.5, label='PDD_poliE_CS')

# plt.plot(eje_z_poliA_CC, 100*dosis_poliA_CC[dosis_poliA_CC.shape[0]//2,dosis_poliA_CC.shape[1]//2,:]/dosis_poliA_CC_norm,'g', linewidth=0.5, label='PDD_poliA_CC')

# plt.plot(eje_z_poliE_CC, 100*dosis_poliE_CC[dosis_poliE_CC.shape[0]//2,dosis_poliE_CC.shape[1]//2,:]/dosis_poliE_CC_norm,'c', linewidth=0.5, label='PDD_poliE_CC')

# plt.plot(eje_z_poliA_CC_EFinf, 100*dosis_poliA_CC_EFinf[dosis_poliA_CC_EFinf.shape[0]//2,dosis_poliA_CC_EFinf.shape[1]//2,:]/dosis_poliA_CC_EFinf_norm,'g.-', linewidth=0.5, label='PDD_poliA_CC_EFinf')

# plt.plot(eje_z_poliE_CC_EFinf, 100*dosis_poliE_CC_EFinf[dosis_poliE_CC_EFinf.shape[0]//2,dosis_poliE_CC_EFinf.shape[1]//2,:]/dosis_poliE_CC_EFinf_norm,'c.-', linewidth=0.5, label='PDD_poliE_CC_EFinf')

# plt.plot(eje_z_kernel_MC, 100*dosis_kernel_MC[dosis_kernel_MC.shape[0]//2,dosis_kernel_MC.shape[1]//2,:]/dosis_kernel_MC_norm,'r.-', linewidth=0.5, label='PDD_kernel_MC')

plt.errorbar(eje_z, 100*dosis_electrones/dosis_electrones_norm, yerr=100*error_electrones*dosis_electrones/dosis_electrones_norm, ls='None', marker='x', c='r', label='Electrones: solo eje del haz')

plt.errorbar(eje_z, 100*dosis_electrones_monoangular, yerr=100*error_electrones_monoangular*dosis_electrones_monoangular/dosis_electrones_monoangular_norm, ls='None', marker='x', c='r', label='Electrones_monoangular: solo eje del haz')

plt.plot(eje_z_kernel_MC_vactowat, 100*dosis_kernel_MC_vactowat[dosis_kernel_MC_vactowat.shape[0]//2,dosis_kernel_MC_vactowat.shape[1]//2,:]/dosis_kernel_MC_vactowat_norm,'c.-', linewidth=0.5, label='PDD_kernel_MC_vactowat')


plt.xlabel('Profundidad [cm]')
plt.legend()
plt.show()



