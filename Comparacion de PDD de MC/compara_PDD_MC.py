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
dosis_electrones = archivo_dosis_electrones['dosis']
error_electrones = archivo_dosis_electrones['errs']

archivo_dosis_fotones = np.load("Varian-20-IAEA_fotones_1000M_paraPDD.npz")
dosis_fotones = archivo_dosis_fotones['dosis']
error_fotones = archivo_dosis_fotones['errs']

archivo_dosis_fot_elect = np.load("Varian-20-IAEA_fot-elec_1000M_paraPDD.npz")
dosis_fot_elect = archivo_dosis_fot_elect['dosis']
error_fot_elect = archivo_dosis_fot_elect['errs']

error_e_diferencia = error_fotones + error_fot_elect

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )

print(eje_z.shape,dosis_electrones.shape)

# -----------------------------------------------------------------------------------------------------------------------------------------

dosis_e_diferencia  = dosis_fot_elect - dosis_fotones
dosis_e_diferencia += dosis_electrones[-1] - dosis_e_diferencia[-1]

dosis_electrones_norm    = np.interp(prof_norm,eje_z,dosis_electrones)
dosis_e_diferencia_norm  = np.interp(prof_norm,eje_z,dosis_e_diferencia)

archivo_dosis_MC = np.load("dosis_MC_1000M.npz")
dosis_MC = archivo_dosis_MC['dosis']
error_MC = archivo_dosis_MC['errs']

dosis_MC = dosis_MC.reshape((17,3,101))
dosis_MC = np.transpose(dosis_MC,(2,1,0))
error_MC = error_MC.reshape((17,3,101))
error_MC = np.transpose(error_MC,(2,1,0))

eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10,5)+2.5 )).ravel()
dosis_MC_norm = np.interp(prof_norm,eje_z_MC,dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:])


# plt.figure()

# plt.errorbar(eje_z, 100*dosis_electrones/dosis_electrones_norm,'b.-', linewidth=0.5, label='Electrones: solo eje del haz')
# plt.errorbar(eje_z, 100*dosis_e_diferencia/dosis_e_diferencia_norm,'m.-', linewidth=0.5, label='Todas las part - fotones: solo eje del haz')

# # plt.plot(eje_z, 100*dosis_fotones/dosis_fotones_norm,'m.-', linewidth=0.5, label='Fotones: eje del haz')
# # plt.plot(eje_z, 100*(dosis_fot_elect - dosis_fotones)/dosis_fot_elect_norm,'g.-', linewidth=0.5, label='(Fotones+electrones) - Fotones: eje del haz')

# plt.errorbar(eje_z_MC, 100*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm,'r.-', linewidth=0.5, label='Electrones: dosis en todo el vol')

# plt.xlabel('Profundidad [cm]')
# plt.legend()
# plt.show()


plt.figure()

plt.errorbar(eje_z, 100*dosis_electrones/dosis_electrones_norm, yerr=100*error_electrones*dosis_electrones/dosis_electrones_norm, ls='None', marker='x', c='b', label='Electrones: solo eje del haz')

plt.errorbar(eje_z, 100*dosis_e_diferencia/dosis_e_diferencia_norm, yerr=100*error_e_diferencia*dosis_e_diferencia/dosis_e_diferencia_norm, ls='None', marker='x', c='orange', label='Todas las part - fotones: solo eje del haz')

# plt.errorbar(eje_z_MC, 100*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm , yerr=100*error_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm, ls='None', marker='o',c='r', label='Electrones: dosis en todo el vol' )

plt.xlabel('Profundidad [cm]')
plt.legend()
plt.show()
