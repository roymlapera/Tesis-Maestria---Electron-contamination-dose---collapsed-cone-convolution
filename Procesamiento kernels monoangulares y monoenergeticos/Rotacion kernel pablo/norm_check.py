import numpy as np    


def norm_check(kernel, rAxis, thetaAxis):

    # para chequear la normalizacion de kernels discretizados cada un grado
    # aunque puede servir para otras discretizaciones
    
    kernel_pol = kernel
    kernel_pol = kernel_pol[:,:180] # por las dudas sea de 360 grados

    theta_axis_distal = np.radians(thetaAxis)
    theta_axis_proximal = np.append([0],theta_axis_distal[:-1])

    r_axis_distal = rAxis
    r_axis_proximal = np.append(0.0, r_axis_distal[:-1]) # en cm

    deltaR = r_axis_distal - r_axis_proximal
        
    
    # Chequeo de la normalizacion

    volume = np.zeros(kernel_pol.shape)

    for k in range(len(theta_axis_distal)):
        volume[:,k] = 2.0*np.pi*( np.cos(theta_axis_proximal[k]) - np.cos(theta_axis_distal[k]))*(r_axis_proximal*deltaR**2+deltaR*r_axis_proximal**2+(deltaR**3)/3 )

    # print kernel_pol.sum()

    kernel_pol_int = kernel_pol*volume

    return kernel_pol_int.sum(), kernel_pol_int
