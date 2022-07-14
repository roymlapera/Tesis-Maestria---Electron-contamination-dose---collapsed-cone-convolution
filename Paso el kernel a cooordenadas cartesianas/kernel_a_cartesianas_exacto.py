# Pasa kernel a cartesianas de manera exacta

import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import sys

def Integral_vol_del_kernel(kernel3D_esf,rAxis,thetaAxis,phiAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion imprime por pantalla el resultado de la integral
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

    for j in range(len(phiAxis)):
        for k in range(len(thetaAxis)):
            volume[:,k,j] = 2*np.pi/len(phiAxis) * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3) 

    Frac_energ_esf = kernel3D_esf*volume

    print('La integral del kernel3D_esf en la esfera de radio 10 da: ' + str(Frac_energ_esf.sum()))

    return volume


def arccoseno(y,x):
	if(x==0): 
			return np.pi/2
	else: 
		return np.arccos(y/x)
	

def Calcula_limites_de_voxel_cart(rAxis,thetaAxis,phiAxis,x,y,z,tam_vox_cart):

	vert = np.array([[x+tam_vox_cart/2,y+tam_vox_cart/2,z+tam_vox_cart/2] , 
			[x+tam_vox_cart/2,y+tam_vox_cart/2,z-tam_vox_cart/2] , 
			[x+tam_vox_cart/2,y-tam_vox_cart/2,z+tam_vox_cart/2] , 
			[x+tam_vox_cart/2,y-tam_vox_cart/2,z-tam_vox_cart/2] ,
			[x-tam_vox_cart/2,y+tam_vox_cart/2,z+tam_vox_cart/2] , 
			[x-tam_vox_cart/2,y+tam_vox_cart/2,z-tam_vox_cart/2] ,
			[x-tam_vox_cart/2,y-tam_vox_cart/2,z+tam_vox_cart/2] , 
			[x-tam_vox_cart/2,y-tam_vox_cart/2,z-tam_vox_cart/2]])

	vert = vert.reshape((8,3))

	r = np.sqrt(vert[:,0]**2+vert[:,1]**2+vert[:,2]**2)
	th = np.arccos(-vert[:,2]/r)
	phi = np.arctan2(vert[:,1], vert[:,0])      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

	ir_min = int(max(np.argwhere(rAxis<r.min())))
	ir_max = int(min(np.argwhere(rAxis>r.max())))

	ith_min = int(max(np.argwhere(thetaAxis<th.min())))
	ith_max = int(min(np.argwhere(thetaAxis>th.max())))

	iphi_min = int(max(np.argwhere(phiAxis<phi.min())))
	iphi_max = int(min(np.argwhere(phiAxis>phi.max())))

	return ir_min,ir_max,ith_min,ith_max,iphi_min,iphi_max


def f(z,y,x,r_min,r_max,th_min,th_max,phi_min,phi_max):

	r = np.sqrt(x*x+y*y+z*z)
	th = arccoseno(-z,r)
	phi = np.arctan2(y, x)      #esta funcion le da un valor pi/2 a new_phi cuando X = 0

	return (r<r_max) * (r>r_min) * (th<th_max) * (th>th_min) * (phi<phi_max) * (phi>phi_min)


def Calcula_frac_energ_en_vox_cart(rAxis,thetaAxis,phiAxis,kernel_esf,x0,y0,z0,tam_vox_cart):

	ir_min,ir_max,ith_min,ith_max,iphi_min,iphi_max = Calcula_limites_de_voxel_cart(rAxis,thetaAxis,phiAxis,x0,y0,z0,tam_vox_cart)

	r_min = rAxis[ir_min]
	r_max = rAxis[ir_max]
	th_min = thetaAxis[ith_min]
	th_max = thetaAxis[ith_max]
	phi_min = phiAxis[iphi_min]
	phi_max = phiAxis[iphi_max]

	frac_energ = 0

	for i in range(ir_min,ir_max,1):
		for j in range(ith_min,ith_max,1):
			for k in range(iphi_min,iphi_max,1):
				frac_energ += integrate.tplquad(lambda z, y, x: f(z,y,x,r_min,r_max,th_min,th_max,phi_min,phi_max),x0-tam_vox_cart/2,x0+tam_vox_cart/2,lambda x: y0-tam_vox_cart/2,lambda x: y0+tam_vox_cart/2,lambda x,y: z0-tam_vox_cart/2,lambda x,y: z0+tam_vox_cart/2)[0] 

	return frac_energ


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'
# file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular_reg_central.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360

kernel = archivo['kernel']
# kernel = kernel[:,:180]  # me saco de encima la mitad del kernel que contiene la misma info

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis'] * np.pi/180
phiAxis   = np.linspace(-180,180,nro_conos_phi) * np.pi/180

kernel_esf, a, b = np.mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    kernel_esf[:,:,k] = kernel		# kernel_esf contiene la fraccion de energia depositada por volumen de voxel esferico

volume = Integral_vol_del_kernel(kernel_esf,rAxis,thetaAxis,phiAxis)

tamano_vox_cart = 0.1

xmin=-0.5
xmax=0.5
ymin=-0.5
ymax=0.5
zmin=-0.5
zmax=-1.5

nro_vox_gra_x = round(abs(xmax-xmin)/tamano_vox_cart)
nro_vox_gra_y = round(abs(ymax-ymin)/tamano_vox_cart)
nro_vox_gra_z = round(abs(zmax-zmin)/tamano_vox_cart)

# #Inicializo los puntos que quiero en el sistema cartesiano

x = np.linspace(xmin,xmax,nro_vox_gra_x)
y = np.linspace(ymin,ymax,nro_vox_gra_y)
z = np.linspace(zmin,zmax,nro_vox_gra_z)

print(x,y,z)

#Obtengo los valores (X,Y,Z) que quiero tener del kernel en las nuevas coordenadas

kernel_cart, a, b = np.meshgrid(y,x,z)

for i,x0 in enumerate(x):
	for j,y0 in enumerate(y):
		for k,z0 in enumerate(z):
			print(i,j,k)
			print(x0,y0,z0)
			kernel_cart[i][j][k] = Calcula_frac_energ_en_vox_cart(rAxis,thetaAxis,phiAxis,kernel_esf,x0,y0,z0,tamano_vox_cart)
			sys.stdout.write("\r{:.2f}".format(100*(k+j*nro_vox_gra_z+i*nro_vox_gra_z*nro_vox_gra_y)/(nro_vox_gra_x*nro_vox_gra_y*nro_vox_gra_z)) + ' %')
			sys.stdout.flush()

plt.figure()
plt.imshow(kernel_cart[:,kernel_cart.shape[1]//2,:], cmap=cm.hot, norm=LogNorm(vmin=10**-7, vmax=10**7))
plt.show()



