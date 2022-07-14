from pylab import *
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt

def spherical2cartesian(r, th, phi, grid, x, y, z, order, nro_voxels_cartesianos, limite):

    # Build relationship between Cartesian and spherical coordinates.
    X, Y, Z = np.meshgrid(x,y,z)

    #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

    X *= limite/(nro_voxels_cartesianos//2)
    Y *= limite/(nro_voxels_cartesianos//2)
    Z *= limite/(nro_voxels_cartesianos//2)

    print('x max =',X.max(),' x min =',X.min())

    #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

    new_r = np.sqrt(X*X+Y*Y+Z*Z)
    new_th = np.arccos(-Z/new_r)    #aca me va a tirar un error cuando (X,Y,Z) = (0,0,0) 
    new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

    #Obtengo funciones a trozos ir, ith y iphi obtenidas de r, th y phi

    ir = interp1d(r, np.arange(len(r)), bounds_error=False, fill_value='extrapolate')
    ith = interp1d(th, np.arange(len(th)), bounds_error=False, fill_value='extrapolate')
    iphi = interp1d(phi, np.arange(len(phi)), bounds_error=False, fill_value='extrapolate')

    #Aplico las funciones obtenidas para encontrar qué indices de la grilla del kernel le corresponde a (X,Y,Z) 

    new_ir = ir(new_r.ravel())
    new_ith = ith(new_th.ravel())
    new_iphi = iphi(new_phi.ravel())

    #Corrijo por si hay algunos valores mal interpolados

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    new_ith[new_th.ravel() > th.max()] = len(th)-1
    new_ith[new_th.ravel() < th.min()] = 0

    new_iphi[new_phi.ravel() > phi.max()] = len(phi)-1
    new_iphi[new_phi.ravel() < phi.min()] = 0

    # Interpolo en las nuevas coordenadas encontradas para conocer su valor de kernel correspondiente

    kernel_cart = map_coordinates(grid, np.array([new_ir, new_ith, new_iphi]), order=order, mode='nearest').reshape(new_r.shape)

    for i in range(nro_voxels_cartesianos):
    	for j in range(nro_voxels_cartesianos):
    		for k in range(nro_voxels_cartesianos):
    			semi_nro_vox = nro_voxels_cartesianos//2
    			if (  (i-semi_nro_vox)**2 + (j-semi_nro_vox)**2 + (k-semi_nro_vox)**2 > (semi_nro_vox)**2):
    				kernel_cart[i,j,k] = 0

    return kernel_cart


# file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular_corregido.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360
nro_voxels = 486                       #siempre tiene que ser P A R, asi NO hay un voxel central

kernel = archivo['kernel']
kernel = kernel[:,:180]

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis']
phiAxis   = np.linspace(-180,180,nro_conos_phi)

#calculo coordenadas de los centros de los voxeles

deltaR = np.zeros(len(rAxis))
deltaR[0] = rAxis[0]
k=0
while k < nro_radios-1:
    deltaR[k+1] = rAxis[k+1]-rAxis[k]
    k+=1
r = rAxis - deltaR/2 

th = ( thetaAxis- 0.5 ) *np.pi/180

phi = ( phiAxis- 0.5 ) *np.pi/180 

#Inicializo kernel3D asi para que sea float64

kernel3D, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    kernel3D[:,:,k] = kernel

# print(kernel3D.sum())

#Inicializo los puntos que quiero en el sistema cartesiano

x = linspace(-nro_voxels//2,nro_voxels//2,nro_voxels)
y = linspace(-nro_voxels//2,nro_voxels//2,nro_voxels)
z = linspace(-nro_voxels//2,nro_voxels//2,nro_voxels)

# Los valores de x, y y z van a ir entre -limite y +limite

limite = 10

tamano_voxel = 4/(nro_voxels//2)

#Chequeo normalizacion del kernel en coordenadas esfericas

deltaR = np.zeros(len(rAxis))
deltaR[0] = rAxis[0]
k=0
while k < nro_radios-1:
    deltaR[k+1] = rAxis[k+1]-rAxis[k]
    k+=1

RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

volume = np.zeros(kernel.shape)
kernel_vol = 0

for k in range(nro_conos_th):
   volume[:,k] = 2*np.pi*(np.cos(ThetaAxis[k]*np.pi/180)-np.cos(thetaAxis[k]*np.pi/180))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

kernel_por_vol = kernel*volume

print('La integral del kernel en la esfera de radio 4, en coord esfericas, da: ' + str(kernel_por_vol.sum()))

# -----------------------------------------------------------------------------------------------------------------------

# Mapeo a coordenadas cartesianas.

order = 1            #orden de interpolacion cuando interpolo los valores en esfericas de (X,Y,Z) dados los valores de r,th y phi
kernel3D = spherical2cartesian(r, th, phi, kernel3D, x, y, z, order, nro_voxels, limite)

# #Correccion de kernel en 60 60 58

# print(kernel3D[60][60][58])
# kernel3D[60][60][58] = 1099873107.88

print('La integral del kernel en la esfera de radio 4, en coord cartesianas, da: ' + str(kernel3D.sum()*tamano_voxel**3))

# -----------------------------------------------------------------------------------------------------------------------

# Grafico el kernel en cartesianas

figure()
imshow(squeeze(kernel3D[:,nro_voxels//2,:]), cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**4))
show()

#Grafico varios cortes del kernel

# figure()
# for i in range(nro_voxels//2-10,nro_voxels//2+10,1):
# 	imshow(squeeze(kernel3D[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-7), vmax=10**4), extent=(-limite,limite,-limite,limite))
# 	print(i)
# 	show()

# -----------------------------------------------------------------------------------------------------------------------

#Guardo el kernel en cordenadas cartesianas. Tamano de la matriz: nro_voxels*nro_voxels*nro_voxels

# np.savez("kernel_poliangular_cartesianas_486_voxeles.npz", kernel = kernel3D)