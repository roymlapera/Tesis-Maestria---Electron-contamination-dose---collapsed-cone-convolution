import lmfit
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import scipy.special
import scipy.ndimage as ndimage

def funcion_ajuste_FE_simple(x,y,l,sigma,alpha):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = a
    d = b

    sigma_x = sigma
    sigma_y = sigma

    cte = 1/(np.pi*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return cte*alpha*( scipy.special.erf((b-x)/raiz_2/sigma_x) - scipy.special.erf((a-x)/raiz_2/sigma_x) )*( scipy.special.erf((d-y)/raiz_2/sigma_y) - scipy.special.erf((c-y)/raiz_2/sigma_y) )
# ---------------------------------------------------------------- DATA -------------------------------------------------------------------------

bins = 121

limite_fantoma = 30 #cm   #el fantoma va de -limite a +limite

archivo_FE = np.load("histo_fluencia_energetica_todos_los_electrones.npz")

Fenergetica = archivo_FE['Fenergetica']

area_pixel = (limite_fantoma*2/bins)**2

Fenergetica /= area_pixel   #ahora la fluencia enrgetica esta en MeV/cm2

x = np.linspace(-limite_fantoma,limite_fantoma,bins)
y = np.linspace(-limite_fantoma,limite_fantoma,bins)

# ---------------------------------------------------------------- FILTRADO -------------------------------------------------------------------------

sigma_de_filtrado = 5
Fenergetica_filtrada = ndimage.gaussian_filter(Fenergetica, sigma=sigma_de_filtrado)

# plt.plot(np.arange(bins),Fenergetica_filtrada[bins//2])
# plt.plot(np.arange(bins),Fenergetica[bins//2], 'r.')

# plt.show()

# plt.imshow(Fenergetica)
# plt.show()
plt.imshow(Fenergetica_filtrada)
plt.show()

#Normalizo FE:

# Fenergetica_filtrada /= Fenergetica_filtrada.sum()*area_pixel

# print(0.5*(Fenergetica_filtrada[:,bins//2+1].min() + Fenergetica_filtrada[bins//2+1,:].min()))

# Fenergetica_filtrada -= 0.5*(Fenergetica_filtrada[:,bins//2+1].min() + Fenergetica_filtrada[bins//2+1,:].min())

# # ---------------------------------------------------------------- AJUSTE -------------------------------------------------------------------------

# X, Y = np.meshgrid(x,y)

# l = 10

# p = lmfit.Parameters()
# p.add('sigma', value = 7, min=5, max=8)
# # p.add('sigma_x', value = 9, min=0, max=10)
# # p.add('sigma_y', value = 9, min=0, max=10)
# p.add('alpha', value = 0.75, min=0, max=1)
# # p.add('C', value = 3E-6, min=0, max=2.9E-5)

# def residual(p):
# 	v = p.valuesdict()
# 	alpha = v['alpha']
# 	# sigma_x = v['sigma_x']
# 	# sigma_y = v['sigma_y']
# 	sigma = v['sigma']
# 	# C = v['C']

# 	res_rel = []

# 	for i in range(bins):
# 		for j in range(bins):
# 			if(abs(X[i][j])<20 and abs(Y[i][j])<20):    #Solo ajusto los pixeles cuya distancia al centro sea menor a 15 cm
# 				res_rel.append( (funcion_ajuste_FE_simple(X[i][j],Y[i][j],l,sigma,alpha) - Fenergetica_filtrada[i][j])/Fenergetica_filtrada[i][j] )

# 	return res_rel



# mi = lmfit.minimize(residual, p, method='leastsq')
# lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)

# fig = plt.figure()

# plt.plot(x, Fenergetica_filtrada[:,61], 'b.', label='FE filtrada X')
# plt.plot(x, funcion_ajuste_FE_simple(x,0,l,mi.params.valuesdict()['sigma'],mi.params.valuesdict()['alpha']), 'b', label='ajuste')
# # plt.plot(x, funcion_ajuste_FE_simple(x,0,8,9,2.8E-2,10), 'g', label='ajuste a ojo')

# fig = plt.figure()

# plt.plot(x, Fenergetica_filtrada[61,:], 'r.', label='FE filtrada Y')
# plt.plot(x, funcion_ajuste_FE_simple(0,y,l,mi.params.valuesdict()['sigma'],mi.params.valuesdict()['alpha']), 'r', label='ajuste')
# # plt.plot(x, funcion_ajuste_FE_simple(0,y,8,9,2.8E-2,10), 'm', label='ajuste a ojo')

# plt.show()