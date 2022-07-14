# CONVOLUCION CON APROXIMACION DE CONO COLAPSADO

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import math
import sys
import scipy.special
import scipy.integrate

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

def Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi):
    '''Importa los valores de kernel y la informacion espacial a partir de un archivo NPZ. Los devuelve en forma de matriz y vectores, respectivamente.'''
    archivo = np.load(file_path + npz_file)

    kernel = archivo['kernel']
    kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

    rAxis     = archivo['rAxis']
    thetaAxis = archivo['thetaAxis'] * np.pi/180
    phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

    kernel3D_esf, a, b = np.mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

    # Asigno el kernel para los angulos phi, por haber simetria de revolucion

    for k in range(nro_conos_phi):
        kernel3D_esf[:,:,k] = kernel

    return kernel3D_esf, rAxis, thetaAxis, phiAxis

def Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion muestra por pantalla el resultado de la integral
    y devuelve la matriz de volumenes con ejes segun rAxis thetaAxis y phiAxis.'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < nro_radios-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel3D_esf.shape)
    ang_solido = np.zeros(kernel3D_esf.shape)

    for j in range(len(phiAxis)):
        for k in range(len(thetaAxis)):
            ang_solido[:,k,j] = 2*np.pi/len(phiAxis) * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))
            volume[:,k,j] = ang_solido[:,k,j] * (RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    Frac_energ_esf = kernel3D_esf*volume

    print('La integral del kernel3D_esf en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume,ang_solido

def rebin(arr, conos_sum_en_theta, conos_sum_en_phi):
    arr = np.array(np.split(arr, indices_or_sections=conos_sum_en_theta ,axis=1)).sum(axis=0)
    arr = np.array(np.split(arr, indices_or_sections=conos_sum_en_phi ,axis=2)).sum(axis=0)
    return arr

def Genera_kernels_cono_colapsado(kernel3D_esf,rAxis,thetaAxis,phiAxis,conos_sumados_en_theta,conos_sumados_en_phi):

    rad,th,phi = Devuelve_coordenadas_centro_voxeles(rAxis,thetaAxis,phiAxis)  #centro de los voxeles del kernel original

    print(len(rad),len(th),len(phi))

    for j,t in enumerate(th):
        kernel3D_esf[:,j,:] *= np.sin(t)

    kernels_cono_colapsado = rebin(kernel3D_esf, conos_sumados_en_theta, conos_sumados_en_phi)

    nro_conos_th_KCC = 90//conos_sumados_en_theta
    nro_conos_phi_KCC = 360//conos_sumados_en_phi

    kernels_cono_colapsado *= (np.pi/2/90)*(np.pi*2/360)

    for k,r in enumerate(rad):
        kernels_cono_colapsado[k,:,:] *= r**2

    return kernels_cono_colapsado

def Genera_kernel_continuo(m,n,r,kernels_cono_colapsado): 
    return  lambda x: interp1d(r, kernels_cono_colapsado[:,m,n], bounds_error=False, fill_value=0, kind='cubic')(x)

def Devuelve_coordenadas_centro_voxeles_KCC(rAxis,thetaAxis,phiAxis,conos_intgrados_th,conos_intgrados_phi):
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    deltaTh = np.zeros(len(thetaAxis))
    deltaTh[0] = thetaAxis[0]
    deltaPhi = np.zeros(len(phiAxis))
    deltaPhi[0] = phiAxis[0]

    k=0
    while k < len(rAxis)-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1
    r = [0] + rAxis
    r = r - deltaR/2

    k=0
    while k < len(thetaAxis)-1:
        deltaTh[k+1] = thetaAxis[k+1]-thetaAxis[k]
        k+=1
    th = [0] + thetaAxis
    th = th - deltaTh/2

    k=0
    while k < len(phiAxis)-1:
        deltaPhi[k+1] = phiAxis[k+1]-phiAxis[k]
        k+=1
    phi = [0] + phiAxis
    phi = phi - deltaPhi/2

    th = th[conos_intgrados_th//2+1::conos_intgrados_th]
    phi = phi[conos_intgrados_phi//2+1::conos_intgrados_phi]

    return r,th,phi

def Devuelve_coordenadas_centro_voxeles(rAxis,thetaAxis,phiAxis):
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    deltaTh = np.zeros(len(thetaAxis))
    deltaTh[0] = thetaAxis[0]
    deltaPhi = np.zeros(len(phiAxis))
    deltaPhi[0] = phiAxis[0]

    k=0
    while k < len(rAxis)-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1
    r = [0] + rAxis
    r = r - deltaR/2

    k=0
    while k < len(thetaAxis)-1:
        deltaTh[k+1] = thetaAxis[k+1]-thetaAxis[k]
        k+=1
    th = [0] + thetaAxis
    th = th - deltaTh/2

    k=0
    while k < len(phiAxis)-1:
        deltaPhi[k+1] = phiAxis[k+1]-phiAxis[k]
        k+=1
    phi = [0] + phiAxis
    phi = phi - deltaPhi/2

    return r,th,phi

def Calcula_dosis_por_conos_colapsados(kernel3D_esf, rAxis,thetaAxis,phiAxis,conos_intgrados_th,conos_intgrados_phi,
                                        kernels_cono_colapsado,
                                        l,sigma,alpha,C,
                                        x,y,z):

    # Encuentro las direcciones de los conos colapsados determinadas por los angulos th y phi

    r,th,phi = Devuelve_coordenadas_centro_voxeles_KCC(rAxis,thetaAxis,phiAxis,conos_intgrados_th,conos_intgrados_phi)

    # Calculo las coordenadas (xs[m,n],ys[m,n]) de los puntos interseccion entre la recta con direccion (th[m],phi[n]) que pasa por el punto (x,y,z) 
    # donde quiero calcular dosis y el plano de superficie del fantoma. Luego evaluo la FE en esos puntos y multiplico por el kernel colapsado evaluado
    # en la distancia entre el punto (xs[m,n],ys[m,n]) y (x,y,z) para cada direccion (th[m],phi[n]).

    z_fantoma = 0
    dosis = 0

    rho = np.arange(0.025,4,0.05)

    count = 0

    for m,t in enumerate(th):
        for n,p in enumerate(phi):
            xs = x-z*np.cos(p)*np.tan(t)
            ys = y-z*np.sin(p)*np.tan(t)
            d = z/np.cos(t) 
            if (abs(xs)<=limite_fantoma and abs(ys)<=limite_fantoma and d<rAxis[-1]):
                count +=1
                dosis += Genera_kernel_continuo(m,n,r,kernels_cono_colapsado)(d) * FE(xs,ys,l,sigma,alpha,C)
        # plt.figure()
        # plt.plot(rho,Genera_kernel_continuo(m,0,r,kernels_cono_colapsado)(rho),'r.', label='theta' + str(m))
        # plt.ylim(0,0.004)
        # plt.legend()
    # plt.plot(r,np.log10(kernels_cono_colapsado[:,m,0]),'b.')

    plt.show()

    print('Nro de puntos en superficie = ' + str(count))

    return dosis

def DevuelveNormaKCC(kernels_cono_colapsado,rAxis):
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < nro_radios-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    kernels_cono_colapsado = kernels_cono_colapsado.sum(axis=2).sum(axis=1)

    norma = 0

    for i,dr in enumerate(deltaR):
        norma += dr*kernels_cono_colapsado[i]

    print(norma)



# ------------------------------------------------------ RESULTADOS DE AJUSTE DE LA FENERGETICA ------------------------------------------------------

sigma = 7.439
alpha = 0.556
l     = 10       # l = TC/2
C     = 8.72E-5

# ------------------------------------------------------- CALCULO DEL KERNEL DE CONO COLAPSADO ------------------------------------------------------- 

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Convolucion cono colapsado/'
# file_path = 'C:/Users/OWner/Documents/Facu\Bariloche 2/Nuevos Scripts/Convolucion cono colapsado/'

npz_file = "fullWater_50M_2125MeV.npz"
# npz_file = "kernel_poliangular_reg_central.npz"

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

limite_fantoma = 30.05 #cm (30 voxeles de tamano tam_vox_cart + medio voxel)

kernel3D_esf, rAxis, thetaAxis, phiAxis = Importa_kernel_de_npz(file_path,npz_file,nro_radios,nro_conos_th,nro_conos_phi)

volume,ang_solido = Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis)

# ---------------------------------------------------------------------------------------------------------- 

nro_radios = 210
nro_conos_th = 90

kernel3D_esf = kernel3D_esf[:nro_radios,:nro_conos_th,:]
rAxis = rAxis[:nro_radios]
thetaAxis = thetaAxis[:nro_conos_th]

print((kernel3D_esf*volume[:nro_radios,:nro_conos_th,:]).sum())    #chequeo normalizacion del la parte del kernel que esta abajo de la superficie del fantoma

conos_intgrados_th  = 5
conos_intgrados_phi = 12

kernels_cono_colapsado = Genera_kernels_cono_colapsado(kernel3D_esf,rAxis,thetaAxis,phiAxis,conos_intgrados_th,conos_intgrados_phi)

DevuelveNormaKCC(kernels_cono_colapsado,rAxis)

# ------------------------------------------------------- CONVOLUCION ------------------------------------------------------- 

tam_vox_cart = 0.1#cm

xmin = -tam_vox_cart/2
xmax = tam_vox_cart/2
ymin = -tam_vox_cart/2
ymax = tam_vox_cart/2
zmin = 0
zmax = 4

nro_vox_gra_x = round(abs(xmax-xmin)/tam_vox_cart)
nro_vox_gra_y = round(abs(ymax-ymin)/tam_vox_cart)
nro_vox_gra_z = round(abs(zmax-zmin)/tam_vox_cart)

# #Inicializo los puntos que quiero en el sistema cartesiano

X = np.linspace(xmin+tam_vox_cart/2,xmax-tam_vox_cart/2,nro_vox_gra_x)
Y = np.linspace(ymin+tam_vox_cart/2,ymax-tam_vox_cart/2,nro_vox_gra_y)
Z = np.linspace(zmin+tam_vox_cart/2,zmax-tam_vox_cart/2,nro_vox_gra_z)

dosis = np.zeros((X.shape[0],Y.shape[0],Z.shape[0]))

# for i,x in enumerate(X):
#     for j,y in enumerate(Y):
#         for k,z in enumerate(Z):
#             dosis[i][j][k] = Calcula_dosis_por_conos_colapsados(kernel3D_esf, rAxis,thetaAxis,phiAxis,conos_intgrados_th,conos_intgrados_phi,kernels_cono_colapsado,l,sigma,alpha,C,x,y,z)
#             sys.stdout.write("\r{:.2f}".format((i*dosis.shape[1]*dosis.shape[2]+j*dosis.shape[2]+k)/(dosis.shape[0]*dosis.shape[1]*dosis.shape[2])*100) + ' %')
#             sys.stdout.flush() 

# # ------------------------------------------------------- GUARDO PDD ------------------------------------------------------- 

# np.savez("dosis_conv_conocolapsado_5_12_fullWater_50M_2125MeV.npz", dosis = dosis)
# np.savez("dosis_conv_conocolapsado_5_15_kernel_poliangular_reg_central_FE1_EFinf.npz", dosis = dosis)

# print(dosis.min(),dosis.max())


# archivo_dosis_CC = np.load("dosis_conv_conocolapsado_3_15_kernel_poliangular_roy_2D.npz")
# dosis_CC = archivo_dosis_CC['dosis']
# print(dosis.shape)

# plt.figure()
# plt.imshow(dosis_CC[:,0,:], cmap=cm.hot, norm=LogNorm(vmin=10**-10, vmax=10**-5))
# # plt.figure()
# # plt.imshow(dosis[:,:,1], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**-5))
# plt.show()


# ------------------------------------------------------- PDD ------------------------------------------------------- 

prof_norm = 0.5 #cm

archivo_dosis_CC = np.load("dosis_conv_conocolapsado_5_12_fullWater_50M_2125MeV.npz")
dosis_CC = archivo_dosis_CC['dosis']

eje_z_CC = np.arange(zmin,zmax,tam_vox_cart) + tam_vox_cart/2
dosis_CC_norm             = np.interp(prof_norm,eje_z_CC,dosis_CC[dosis_CC.shape[0]//2,dosis_CC.shape[1]//2,:])

# # # ------------------------------------------------------------------------------------------------------------ 

# archivo_dosis_MC = np.load("dosis_MC_1000M.npz")
# dosis_MC = archivo_dosis_MC['dosis']

# dosis_MC = dosis_MC.reshape((17,3,101))
# dosis_MC = np.transpose(dosis_MC,(2,1,0))

# eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10,5)+2.5 )).ravel()
# dosis_MC_norm = np.interp(prof_norm,eje_z_MC,dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:])

# # # ----------------------------------------

# archivo_dosis_electrones_monoangular = np.load("Varian-20-IAEA_monoangular_electrones_1000M_paraPDD.npz")
# dosis_electrones_monoangular = archivo_dosis_electrones_monoangular['dosis']
# error_electrones_monoangular = archivo_dosis_electrones_monoangular['errs']

# dosis_electrones_monoangular = dosis_electrones_monoangular.reshape((34,11,11))

# error_electrones_monoangular = error_electrones_monoangular.reshape((34,11,11))

# dosis_electrones_monoangular = dosis_electrones_monoangular[:,dosis_electrones_monoangular.shape[1]//2,dosis_electrones_monoangular.shape[2]//2]
# error_electrones_monoangular = error_electrones_monoangular[:,error_electrones_monoangular.shape[1]//2,error_electrones_monoangular.shape[2]//2]

# eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
# dosis_electrones_monoangular_norm    = np.interp(prof_norm,eje_z,dosis_electrones_monoangular)

# # # ----------------------------------------

archivo_dosis_monoAmonoE = np.load("Varian-20-IAEA_monoangular_electrones_1000M_paraPDD.npz")
dosis_monoAmonoE = archivo_dosis_monoAmonoE['dosis']
error_monoAmonoE = archivo_dosis_monoAmonoE['errs']

dosis_monoAmonoE = dosis_monoAmonoE.reshape((34,11,11))

error_monoAmonoE = error_monoAmonoE.reshape((34,11,11))

dosis_monoAmonoE = dosis_monoAmonoE[:,dosis_monoAmonoE.shape[1]//2,dosis_monoAmonoE.shape[2]//2]
error_monoAmonoE = error_monoAmonoE[:,error_monoAmonoE.shape[1]//2,error_monoAmonoE.shape[2]//2]

eje_z = np.concatenate( (np.arange(0,1,0.05)+0.05/2,np.arange(1,2,0.1)+0.1/2,np.arange(2,4,0.5)+0.5/2) )
dosis_monoAmonoE_norm    = np.interp(prof_norm,eje_z,dosis_monoAmonoE)

# # # -------------------------------------------------------------------------------------------------------------- 

plt.figure()

plt.plot(eje_z_CC, 100*dosis_CC[dosis_CC.shape[0]//2,dosis_CC.shape[1]//2,:]/dosis_CC_norm,'b.-', linewidth=0.5, label='PDD_cono_colapsado')

# # plt.plot(eje_z_CC, 100*dosis_CC_reg_central[dosis_CC_reg_central.shape[0]//2,dosis_CC_reg_central.shape[1]//2,:]/dosis_CC_norm_reg_central,'m.-', linewidth=0.5, label='PDD_cono_colapsado: reg_central')

# # plt.plot(eje_z_MC, 100*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm,'r.-', linewidth=0.5, label='PDD_MC')

# plt.plot(eje_z, 100*dosis_electrones_monoangular/dosis_electrones_monoangular_norm,'m', linewidth=0.5, label='PDD_MC')

plt.plot(eje_z, 100*dosis_monoAmonoE/dosis_monoAmonoE_norm,'m', linewidth=0.5, label='PDD_MC_monoAmonoE')

plt.xlabel('Profundidad [cm]')
plt.legend()
plt.show()

# # # ---------------------------------------------------------- PDD COMPARA POLIANGULAR POLIENERGETICO ---------------------------------------------------------

# # # prof_norm = 1 #cm

# # # # ----------------------------------------

# # # archivo_dosis_poliA = np.load(file_path + "Final/dosis_conv_conocolapsado_5_15_kernel_poliangular_reg_central.npz")
# # # dosis_poliA = archivo_dosis_poliA['dosis']

# # # eje_z_poliA = Z
# # # dosis_poliA_norm = np.interp(prof_norm,eje_z_poliA,dosis_poliA[dosis_poliA.shape[0]//2,dosis_poliA.shape[1]//2,:])

# # # # ----------------------------------------

# # # archivo_dosis_poliE = np.load(file_path + "Final/dosis_conv_conocolapsado_3_15_kernel_polienergetico_roy.npz")
# # # dosis_poliE = archivo_dosis_poliE['dosis']

# # # eje_z_poliE = Z
# # # dosis_poliE_norm = np.interp(prof_norm,eje_z_poliE,dosis_poliE[dosis_poliE.shape[0]//2,dosis_poliE.shape[1]//2,:])

# # # # ----------------------------------------

# # # archivo_dosis_MC = np.load(file_path + "Final/dosis_MC_1000M.npz")
# # # dosis_MC = archivo_dosis_MC['dosis']

# # # dosis_MC = dosis_MC.reshape((17,3,101))
# # # dosis_MC = np.transpose(dosis_MC,(2,1,0))

# # # eje_z_MC = np.hstack(( np.arange(0,2,0.2)+0.1,np.arange(2,5,0.5)+0.25,np.arange(5,10,5)+2.5 )).ravel()
# # # dosis_MC_norm = np.interp(prof_norm,eje_z_MC,dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:])

# # # # ----------------------------------------

# # # plt.figure()

# # # plt.plot(eje_z_poliA, 100*dosis_poliA[dosis_poliA.shape[0]//2,dosis_poliA.shape[1]//2,:]/dosis_poliA_norm,'b.-', linewidth=0.5, label='PDD_poliA')

# # # plt.plot(eje_z_poliE, 100*dosis_poliE[dosis_poliE.shape[0]//2,dosis_poliE.shape[1]//2,:]/dosis_poliE_norm,'m.-', linewidth=0.5, label='PDD_poliE')

# # # plt.plot(eje_z_MC, 100*dosis_MC[dosis_MC.shape[0]//2,dosis_MC.shape[1]//2,:]/dosis_MC_norm,'r.-', linewidth=0.5, label='PDD_MC')

# # # plt.xlabel('Profundidad [cm]')
# # # plt.legend()
# # # plt.show()








