import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys
import scipy.special
import scipy.fftpack

# -----------------------------------------------------------------------------------------------------------------------------------------

archivo_dosis_electrones1 = np.load("Varian-20-IAEA_monoA_y_monoE0875_electrones_1000M_paraPDD_2.npz")
dosis_econt1 = archivo_dosis_electrones1['dosis']
error_econt1 = archivo_dosis_electrones1['errs']

print(dosis_econt1.shape)

dosis_econt1 = dosis_econt1.reshape((34,11,11))
error_econt1 = error_econt1.reshape((34,11,11))

dosis_econt1 = dosis_econt1[:,dosis_econt1.shape[1]//2,dosis_econt1.shape[2]//2]
error_econt1 = error_econt1[:,error_econt1.shape[1]//2,error_econt1.shape[2]//2]


archivo_dosis_electrones2 = np.load("Varian-20-IAEA_monoangular_y_monoenergetica_electrones_1000M_paraPDD.npz")
dosis_econt2 = archivo_dosis_electrones2['dosis']
error_econt2 = archivo_dosis_electrones2['errs']

dosis_econt2 = dosis_econt2.reshape((34,11,11))
error_econt2 = error_econt2.reshape((34,11,11))

dosis_econt2 = dosis_econt2[:,dosis_econt2.shape[1]//2,dosis_econt2.shape[2]//2]
error_econt2 = error_econt2[:,error_econt2.shape[1]//2,error_econt2.shape[2]//2]


archivo_dosis_electrones3 = np.load("Varian-20-IAEA_monoA_y_monoE5125_electrones_1000M_paraPDD.npz")
dosis_econt3 = archivo_dosis_electrones3['dosis']
error_econt3 = archivo_dosis_electrones3['errs']

dosis_econt3 = dosis_econt3.reshape((34,11,11))
error_econt3 = error_econt3.reshape((34,11,11))

dosis_econt3 = dosis_econt3[:,dosis_econt3.shape[1]//2,dosis_econt3.shape[2]//2]
error_econt3 = error_econt3[:,error_econt3.shape[1]//2,error_econt3.shape[2]//2]

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )

# -----------------------------------------------------------------------------------------------------------------------------------------

archivo_dosis_CC1 = np.load("dosis_conv_conocolapsado_kernel0875.npz")
dosis_CC1 = archivo_dosis_CC1['dosis']

archivo_dosis_CC2 = np.load("dosis_conv_conocolapsado_kernel2125.npz")
dosis_CC2 = archivo_dosis_CC2['dosis']

archivo_dosis_CC3 = np.load("dosis_conv_conocolapsado_kernel5125.npz")
dosis_CC3 = archivo_dosis_CC3['dosis']

zmin = 0.05
zmax = 3

tam_vox_cart = 0.05#cm

nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)
eje_z_CC = np.arange(zmin,zmax,tam_vox_cart) + tam_vox_cart/2  #-


# -----------------------------------------------------------------------------------------------------------------------------------------

prof_norm = 0.6 #cm

dosis_econt_norm1    = np.interp(prof_norm,eje_z,dosis_econt1)
dosis_CC1_norm = np.interp(prof_norm,eje_z_CC,dosis_CC1[dosis_CC1.shape[0]//2,dosis_CC1.shape[1]//2,:])

prof_norm = 0.45 #cm

dosis_econt_norm2    = np.interp(prof_norm,eje_z,dosis_econt2)
dosis_CC2_norm = np.interp(prof_norm,eje_z_CC-0.02,dosis_CC2[dosis_CC2.shape[0]//2,dosis_CC2.shape[1]//2,:])

prof_norm = 1.25 #cm

dosis_econt_norm3    = np.interp(prof_norm,eje_z,dosis_econt3)
dosis_CC3_norm = np.interp(prof_norm,eje_z_CC,dosis_CC3[dosis_CC3.shape[0]//2,dosis_CC3.shape[1]//2,:])


# -----------------------------------------------------------------------------------------------------------------------------------------

plt.figure()

# plt.plot(eje_z_CC, 100*dosis_CC1[dosis_CC1.shape[0]//2,dosis_CC1.shape[1]//2,:]/dosis_CC1_norm,'b.-', linewidth=0.5, label='CC: E cin = 0.875 MeV')

plt.plot(eje_z_CC-0.02, 100*dosis_CC2[dosis_CC2.shape[0]//2,dosis_CC2.shape[1]//2,:]/dosis_CC2_norm,'b.-', linewidth=0.5, label='CC: E cin = 2.125 MeV')

plt.plot(eje_z_CC, 100*dosis_CC3[dosis_CC3.shape[0]//2,dosis_CC3.shape[1]//2,:]/dosis_CC3_norm,'g.-', linewidth=0.5, label='CC: E cin = 5.125 MeV')

# plt.errorbar(eje_z, 100*dosis_econt1/dosis_econt_norm1, yerr=100*error_econt1*dosis_econt1/dosis_econt_norm1, linewidth=0.5, c='b', label='MC')

plt.errorbar(eje_z, 100*dosis_econt2/dosis_econt_norm2, yerr=100*error_econt2*dosis_econt2/dosis_econt_norm2, linewidth=0.5, c='b', label='MC: E cin = 2.125 MeV')

plt.errorbar(eje_z, 100*dosis_econt3/dosis_econt_norm3, yerr=100*error_econt3*dosis_econt3/dosis_econt_norm3, linewidth=0.5, c='g', label='MC: E cin = 5.125 MeV')

plt.ylabel('PDD [%]')
plt.xlabel('Profundidad [cm]')
plt.legend()
plt.grid(True)
plt.show()

