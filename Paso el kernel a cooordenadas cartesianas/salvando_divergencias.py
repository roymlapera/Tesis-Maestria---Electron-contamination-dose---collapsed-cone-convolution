import sys
import numpy as np
import matplotlib.pyplot as plt

directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_kernels as ker

nro_ptos = 50


def Devuelve_radio_inferior(rAxis,l):
    aux = rAxis
    aux = np.delete(aux, np.where(rAxis>l))
    return aux.max()

def Devuelve_radio_superior(rAxis,l):
    aux = rAxis
    aux = np.delete(aux, np.where(rAxis<l))
    return aux.min()

def Devuelve_pto_random_en_voxel_esferico(rAxis, thetaAxis, phiAxis, i, j, k):
    point = [ np.random.uniform(rAxis[i],rAxis[i+1]), np.random.uniform(thetaAxis[j],thetaAxis[j+1]), np.random.uniform(phiAxis[k],phiAxis[k+1]) ]
    return point



# -------------------------------------------------------------------------------------------------------------------------------------------------------

# Cargo el kernel poliangular y lo hago volumetrico

# file_path = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Scripts/Paso el kernel a cooordenadas cartesianas/'

file_path = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Nuevos Scripts/Paso el kernel a cooordenadas cartesianas/'

file_kernel = "kernel_poliangular_corregido.npz"

archivo = np.load(file_path + file_kernel)

nro_radios = 210
nro_conos_th = 180
nro_conos_phi = 360
nro_voxels = 81                      #siempre tiene que ser impar, asi hay un voxel central

limite = 4

kernel = archivo['kernel']
kernel = kernel[:,:180]

rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis'] * np.pi/180
phiAxis   = np.linspace(1,nro_conos_phi, nro_conos_phi) * np.pi/180

#Inicializo kernel3D asi para que sea float64

kernel3D, a, b = np.mgrid[0:nro_radios//2:0.5,0:nro_conos_th//2:0.5,0:nro_conos_phi//2:0.5]   #lo hago asi para que kernel3D sea de float64 y no int32

# Asigno el kernel para los angulos phi, por haber simetria cilindrica

for k in range(nro_conos_phi):
    kernel3D[:,:,k] = kernel

# Ahora uso kernel 3D:

forma_kernel3D = kernel3D.shape

# -------------------------------------------------------------------------------------------------------------------------------------------------------

l = limite / ((nro_voxels-1)//2)

print('Lado del voxel cartesiano: ' + str(l) + ' cm')

r_inf = Devuelve_radio_inferior(rAxis,l)
r_sup = Devuelve_radio_superior(rAxis,np.sqrt(2)*l)

print('Radio inf: ', r_inf, ', Radio sup: ', r_sup)

# Me quedo con la parte del kernel cuyos radios son menores que el radio superior (saco los que incluyen al cubo completamente)

kernel3D = kernel3D[:np.where(rAxis==r_sup)[0][0],:,:]

# MC

nro_ptos_dentro_del_cubo = 0

alpha = np.zeros(kernel3D.shape)             #Factor de correccion de la contribucion de energia en el cubo segun el area del voxel esferico dentro del cubo
count = 0

for i in range(np.where(rAxis==r_sup)[0][0]):
    for j in range(nro_conos_th):
        for k in range(nro_conos_phi):

            if rAxis[i]<r_inf:
                    alpha[i,j,k] = nro_ptos     # Luego voy a dividir por nro_ptos para que me de alpha = 1
            else:
                points = []
                for p in range(nro_ptos):
                    point = Devuelve_pto_random_en_voxel_esferico([0] + list(rAxis), [0] + list(thetaAxis), [0] + list(phiAxis), i, j, k)
                    points.append(point)

                points = np.transpose(np.array(points))

                # print(points.shape)
                x = points[0] * np.sin(points[1]) * np.cos(points[1]) 
                y = points[0] * np.sin(points[1]) * np.sin(points[2]) 
                z = points[0] * np.cos(points[1])

                # print(len(x))

                esta_adentro_del_cubo = (abs(x)<l) * (abs(y)<l) * (abs(z)<l)

                alpha[i,j,k] = esta_adentro_del_cubo.sum()

            count += 1
            sys.stdout.write("\r{0}".format(1.0*count/(np.where(rAxis==r_sup)[0][0]*nro_conos_th*(nro_conos_phi-1))*100))
            sys.stdout.flush()  

alpha /= nro_ptos

# --------------------------------------------------------------------------------------------------------------------------------------------------------------

# np.savetxt('alpha.txt', alpha[:,:,0])

# alpha_guardado = np.loadtxt('alpha.txt')

# rad, theta = np.meshgrid([0] + list(rAxis[:4]), [0] + list(thetaAxis)) #rectangular plot of polar data
# X = rad*np.cos(theta)
# Y = rad*np.sin(theta)

# fig = plt.figure()
# plt.pcolormesh(X, Y, alpha_guardado.T, cmap='PuBu_r', vmin=0, vmax=10) #X,Y & data2D must all be same dimensions
# plt.axis('equal')
# plt.show()

# --------------------------------------------------------------------------------------------------------------------------------------------------------------

producto = kernel3D.ravel() * alpha.ravel()
valor_kernel_pto_central = producto.sum()

print('\n',valor_kernel_pto_central)

