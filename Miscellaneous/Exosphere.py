import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter
from mpl_toolkits.mplot3d import Axes3D

x_min = -5
x_max = 5
y_min = -5
y_max = 5
z_min = -5
z_max = 5

dx = 0.001
dy = 0.001
dz = 0.001
xvalues = np.arange(x_min,x_max+1,dx)
yvalues = np.arange(y_min,y_max+1,dy)
zvalues = np.arange(z_min,z_max+1,dz)
xbox, ybox = np.meshgrid(xvalues,yvalues)
rad = np.sqrt(xbox*xbox + ybox*ybox)

c_1 = 4.2E-5
c_2 = 31

n_exo = np.zeros_like(rad)
mask = rad > 0.5
n_exo[mask] = c_1 * np.exp(c_2/rad[mask])

        
""" fig = plt.figure()
ax = plt.gca()
x_mid = int(np.shape(xbox)[0]/2)
levs = np.linspace(n_exo.min(),n_exo.max(),1000)
yz_plane = ax.contourf(ybox[x_mid,:,:],zbox[x_mid,:,:],n_exo[x_mid,:,:],cmap='plasma',levels=levs)
fig.colorbar(yz_plane,label=r'Density (cm$^{-3}$)')
ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
ax.set_xlabel(r'$y$ ($R_{U}$)')
ax.set_ylabel(r'$z$ ($R_{U}$)')
plt.title(r"Exosphere $y$-$z$ Projection") """

fig = plt.figure()
ax = plt.gca()
levs = np.linspace(n_exo.min(),n_exo.max(),1000)
yz_plane = ax.contourf(xbox,ybox,n_exo,cmap = 'plasma')
fig.colorbar(yz_plane,label=r'Density (cm$^{-3}$)')
ax.set_xlabel(r'$y$ ($R_{U}$)')
ax.set_ylabel(r'$z$ ($R_{U}$)')
plt.title(r"Exosphere $y$-$z$ Projection")

plt.show()


