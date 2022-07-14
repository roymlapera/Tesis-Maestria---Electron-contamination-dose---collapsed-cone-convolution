import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys
import scipy.special
import scipy.fftpack

# -----------------------------------------------------------------------------------------------------------------------------------------

prof_norm = 1 #cm

archivo_dosis_electrones = np.load("Varian-20-IAEA_electrones_1000M_paraPDD.npz")
dosis_econt = archivo_dosis_electrones['dosis']
error_econt = archivo_dosis_electrones['errs']

# dosis_econt = dosis_econt.reshape((34,3,3))
# error_econt = error_econt.reshape((34,3,3))

# dosis_econt = dosis_econt[:,dosis_econt.shape[1]//2,dosis_econt.shape[2]//2]
# error_econt = error_econt[:,error_econt.shape[1]//2,error_econt.shape[2]//2]

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )

# -----------------------------------------------------------------------------------------------------------------------------------------

archivo_dosis_CS1 = np.load("dosis_conv_cart__kernel_c_interp.npz")
dosis_CS1 = archivo_dosis_CS1['dosis']

archivo_dosis_CS2 = np.load("dosis_conv_cart__kernel_s_interp.npz")
dosis_CS2 = archivo_dosis_CS2['dosis']

zmin = 0
zmax = 3

tam_vox_cart = 0.1

nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)
eje_z_CS = np.linspace(zmin+tam_vox_cart/2,zmax-tam_vox_cart/2,nro_vox_gra_z)

# -----------------------------------------------------------------------------------------------------------------------------------------

prof_norm = 0.5 #cm

dosis_econt_norm    = np.interp(prof_norm,eje_z,dosis_econt)
dosis_CS1_norm = np.interp(prof_norm,eje_z_CS,dosis_CS1[dosis_CS1.shape[0]//2,dosis_CS1.shape[1]//2,:])
dosis_CS2_norm = np.interp(prof_norm,eje_z_CS,dosis_CS2[dosis_CS2.shape[0]//2,dosis_CS2.shape[1]//2,:])

# -----------------------------------------------------------------------------------------------------------------------------------------

plt.figure()

plt.plot(eje_z_CS, 100*dosis_CS1[dosis_CS1.shape[0]//2,dosis_CS1.shape[1]//2,:]/dosis_CS1_norm,'bx-', linewidth=0.5, label='Kernel: c/ interpolación')

plt.plot(eje_z_CS, 100*dosis_CS2[dosis_CS2.shape[0]//2,dosis_CS2.shape[1]//2,:]/dosis_CS2_norm,'gx-', linewidth=0.5, label='Kernel: s/ interpolación')

plt.errorbar(eje_z, 100*dosis_econt/dosis_econt_norm, yerr=100*error_econt*dosis_econt/dosis_econt_norm, linewidth=0.5, c='r', label='MC')

plt.ylabel('PDD [%]')
plt.xlabel('Profundidad [cm]')
plt.legend()
plt.grid(True)
plt.show()

