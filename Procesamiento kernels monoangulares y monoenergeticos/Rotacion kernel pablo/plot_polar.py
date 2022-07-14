import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

def plotPolarContour(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax, drawType):

	if drawType == 'full':
		# flippedkernel2draw = np.fliplr(kernel2draw)
		# kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
		thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw))
	deltaTheta[0] = thetaAxis2draw[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw[k+1]-thetaAxis2draw[k]
		k+=1

	thetaAxis2draw=thetaAxis2draw-deltaTheta/2
	thetaAxis2draw = np.radians(thetaAxis2draw)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw, rAxis2draw)

	# kernel2draw = scipy.ndimage.zoom(kernel2draw, 3)
	# rAxisMesh = scipy.ndimage.zoom(rAxisMesh, 3)
	# thetaAxisMesh = scipy.ndimage.zoom(thetaAxisMesh, 3)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, polar = True)

	# polarImg = ax.pcolor(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxis2draw, rAxis2draw, kernel2draw, cmap=cm.jet, norm=colors.LogNorm())
	# kernel2draw.max()*0.01
	# L = [kernel2draw.max()*0.001, kernel2draw.max()*0.01]

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw,L, cmap=cm.hot, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(False)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_rotado_45grados.png', dpi = 200, bbox_inches='tight')
	plt.show()