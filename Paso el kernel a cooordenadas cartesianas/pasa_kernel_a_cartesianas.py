from pylab import *
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib import pyplot as plt

def spherical2cartesian(r, th, phi, grid, x, y, z, order, nro_voxels_cartesianos):

    #Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas
    X, Y, Z = np.meshgrid(x,y,z)

    print(x.shape,y.shape,z.shape)

    #Paso (X,Y,Z) a coordenadas esfericas: (new_r,new_th,new_phi)

    new_r = np.sqrt(X*X+Y*Y+Z*Z)
    new_th = np.arccos(-Z/new_r)    #aca me va a tirar un error cuando (X,Y,Z) = (0,0,0) 
    new_phi = np.arctan2(Y, X)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

    print(new_r.min(),new_r.max())
    print(new_th.min(),new_th.max())
    print(new_phi.min(),new_phi.max())

    #En este punto tengo singularidad por ser new_r nulo

    new_th[(nro_voxels_cartesianos-1)//2,(nro_voxels_cartesianos-1)//2,(nro_voxels_cartesianos-1)//2] = 0

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

    print(kernel_cart.shape)

    for i in range(nro_voxels_cartesianos):
    	for j in range(nro_voxels_cartesianos):
    		for k in range(nro_voxels_cartesianos):
    			if (  i**2 + j**2 + k**2 > (nro_voxels_cartesianos)**2):
    				kernel_cart[i,j,k] = 0

    return kernel_cart

# --------------------------------------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular_roy_final_final.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

kernel = archivo['kernel']
kernel = kernel[:,:180]

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis']
phiAxis   = np.linspace(-180,180,nro_conos_phi)

# --------------------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------------------

#Inicializo kernel3D asi para que sea float64

kernel3D, a, b = mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    kernel3D[:,:,k] = kernel

# print(kernel3D.sum())

# --------------------------------------------------------------------------------------------------------------

#Chequeo normalizacion

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
   volume[:,k] = 2*np.pi/360*(np.cos(ThetaAxis[k]*np.pi/180)-np.cos(thetaAxis[k]*np.pi/180))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

kernel_por_vol = kernel*volume

print('La integral del kernel esf en la esfera de radio 4, en coord esfericas, da: ' + str(kernel_por_vol.sum()*360))

kernel /= (kernel_por_vol.sum()*360)

kernel_por_vol = kernel*volume

print('La integral del kernel esf en la esfera de radio 4, en coord esfericas, da: ' + str(kernel_por_vol.sum()*360))

# --------------------------------------------------------------------------------------------------------------

#Inicializo los puntos que quiero en el sistema cartesiano

xmin = -0.05
xmax =  3.05
ymin = -0.05
ymax =  3.05
zmin =  0
zmax = -3.1

tamano_voxel = 0.1

# #Inicializo los puntos que quiero en el sistema cartesiano

nro_vox_gra_x = round(abs(xmax-xmin)/tamano_voxel)
nro_vox_gra_y = round(abs(ymax-ymin)/tamano_voxel)
nro_vox_gra_z = round(abs(zmax-zmin)/tamano_voxel)

print(nro_vox_gra_x)

# #Inicializo los puntos que quiero en el sistema cartesiano

x = np.linspace(xmin+tamano_voxel/2,xmax-tamano_voxel/2,nro_vox_gra_x)
y = np.linspace(ymin+tamano_voxel/2,ymax-tamano_voxel/2,nro_vox_gra_y)
z = np.linspace(zmin-tamano_voxel/2,zmax+tamano_voxel/2,nro_vox_gra_z)

print(x)
print(y)
print(z)

# -----------------------------------------------------------------------------------------------------------------------

# Mapeo a coordenadas cartesianas.

order = 1            #orden de interpolacion cuando interpolo los valores en esfericas de (X,Y,Z) dados los valores de r,th y phi
kernel_cart = spherical2cartesian(r, th, phi, kernel3D, x, y, z, order,nro_vox_gra_x)

# #Correccion de kernel en 60 60 58

# print(kernel3D[60][60][58])
# kernel3D[60][60][58] = 1099873107.88

# print('La integral del kernel cart en la esfera de radio 4, en coord cartesianas, da: ' + str((kernel_cart.sum())*tamano_voxel**3))

# -----------------------------------------------------------------------------------------------------------------------

# Grafico el kernel en cartesianas

figure()
imshow(kernel3D[:,:,0], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))
show()

#Grafico varios cortes del kernel

# figure()
# for i in range((nro_voxels-1)//2-10,(nro_voxels-1)//2+10,1):
# 	imshow(squeeze(kernel3D[i,:,:]), cmap=cm.hot, norm=LogNorm(vmin=10**(-7), vmax=10**4), extent=(-limite,limite,-limite,limite))
# 	print(i)
# 	show()

# print(kernel3D[i,:,:])
# -----------------------------------------------------------------------------------------------------------------------

#Guardo el kernel en cordenadas cartesianas. Tamano de la matriz: nro_voxels*nro_voxels*nro_voxels

# np.savez("kernel_cartesiano_roy_final_final.npz", kernel = kernel3D)