import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.spatial import distance


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

r0 = 16
K = 0.6

x_min = -40
x_max = 40
dx = 1
y_min = -40
y_max = 40
dy = 1
z_min = -40
z_max = 40
dz = 1

X_grid, Y_grid, Z_grid = np.mgrid[x_min:x_max:dx, y_min:y_max:dy, z_min:z_max:dz]
rad = (X_grid**2 + Y_grid**2 + Z_grid**2)

m_theta =np.linspace(7/8*np.pi,0,80)
theta = np.linspace(-np.pi, 0, num=100)
exclude=np.pi
theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
theta = theta_exc
m_theta = theta
circle = np.linspace(0,2*np.pi,80)
circle_mpause, _x = np.meshgrid(circle,X_grid[:,0,0])
_r = r0*(2./(1.0+np.cos(m_theta)))**K
xmp = _r*np.cos(m_theta)
interpfunc = interpolate.InterpolatedUnivariateSpline(xmp, _r)
rmp = interpfunc(_x)
r_mp = _r
r_open = np.sqrt(rmp**2 - _x**2)
_y = r_open*np.sin(circle_mpause)
_z = r_open*np.cos(circle_mpause)        
x_mp = _x
y_mp = _y
z_mp = _z
ax.plot_surface(x_mp, y_mp, z_mp, color='red', edgecolor='none',alpha=0.2)

K = 0.88
r0 = 20
m_theta =np.linspace(3/4*np.pi,0,80)
circle = np.linspace(0,2*np.pi,80)
circle_mpause, _x = np.meshgrid(circle,X_grid[:,0,0])
_r = r0*(2./(1.0+np.cos(m_theta)))**K
xmp = _r*np.cos(m_theta)
interpfunc = interpolate.InterpolatedUnivariateSpline(xmp, _r)
r_bs = _r
rmp = interpfunc(_x)
r_open = np.sqrt(rmp**2 - _x**2)
_y = r_open*np.sin(circle_mpause)
_z = r_open*np.cos(circle_mpause)  
x_bs = _x
y_bs = _y
z_bs = _z
ax.plot_surface(x_bs, y_bs, z_bs, color='blue', edgecolor='none',alpha=0.2)

ax.set_xlabel(r'$x$ ($R_{U}$)')
ax.set_ylabel(r'$y$ ($R_{U}$)')
ax.set_zlabel(r'$z$ ($R_{U}$)')


###########################################################################################

# Flatten the surface data for magnetopause
x_mp_flat = x_mp.flatten()
y_mp_flat = y_mp.flatten()
z_mp_flat = z_mp.flatten()

# Flatten the surface data for bow shock
x_bs_flat = x_bs.flatten()
y_bs_flat = y_bs.flatten()
z_bs_flat = z_bs.flatten()

# Stack the magnetopause and bow shock points together
points_mp = np.vstack([x_mp_flat, y_mp_flat, z_mp_flat]).T
points_bs = np.vstack([x_bs_flat, y_bs_flat, z_bs_flat]).T

# Remove any NaN points before using griddata
points_mp = points_mp[np.all(np.isfinite(points_mp), axis=1)]
points_bs = points_bs[np.all(np.isfinite(points_bs), axis=1)]

mask_bs = np.isfinite(griddata(points_bs, np.ones(len(points_bs)), (X_grid, Y_grid, Z_grid), method='linear'))
mask_mp = np.isfinite(griddata(points_mp, np.ones(len(points_mp)), (X_grid, Y_grid, Z_grid), method='linear'))

values_bs = np.ones(len(points_bs))  # Placeholder values for bs points
interpolated_bs = griddata(
    points_bs, values_bs, (X_grid, Y_grid, Z_grid), method='linear'
)

values_mp = np.ones(len(points_mp))  # Placeholder values for bs points
interpolated_mp = griddata(
    points_mp, values_mp, (X_grid, Y_grid, Z_grid), method='linear'
)

density_grid = np.zeros_like(X_grid)
density_grid[(interpolated_bs >= rad) & (interpolated_mp <= rad)]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points_bs[:, 0], points_bs[:, 1], points_bs[:, 2], color='blue', alpha=0.1, label='points_bs')
ax.scatter(points_mp[:, 0], points_mp[:, 1], points_mp[:, 2], color='red', alpha=0.1, label='points_mp')
ax.set_xlabel(r'$x$ ($R_{U}$)')
ax.set_ylabel(r'$y$ ($R_{U}$)')
ax.set_zlabel(r'$z$ ($R_{U}$)')

print(interpolated_mp.max())
print(interpolated_mp.min())
print(np.shape(interpolated_mp))
print(type(interpolated_mp))

###########################################################################################

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#density_grid[np.isfinite(interpolated_bs) & np.isfinite(interpolated_mp)] = True
density_grid = np.zeros_like(rad)
""" density_grid[np.isfinite(interpolated_bs)] = True
density_grid[np.isfinite(interpolated_mp)] = False """
""" for i,j,k in range(0,np.shape(density_grid)[0]):
    if density_grid == True:
        density_grid[i,j,k] = 0.1 """

density_grid[np.isfinite(griddata(points_bs, np.ones(len(points_bs)), (X_grid, Y_grid, Z_grid), method='linear')) & 
             ~np.isfinite(griddata(points_mp, np.ones(len(points_mp)), (X_grid, Y_grid, Z_grid), method='linear'))] = True

indices = np.where(density_grid)
x_points = X_grid[indices]
y_points = Y_grid[indices]
z_points = Z_grid[indices]

ax.scatter(x_points,y_points,z_points,color='green',alpha=0.1)
ax.plot_surface(x_mp, y_mp, z_mp, color='red', edgecolor='none', alpha=0.2)
ax.plot_surface(x_bs, y_bs, z_bs, color='blue', edgecolor='none', alpha=0.2)

ax.set_xlabel(r'$x$ ($R_{U}$)')
ax.set_ylabel(r'$y$ ($R_{U}$)')
ax.set_zlabel(r'$z$ ($R_{U}$)')

###########################################################################################

valid = ~np.isnan(interpolated_mp) & ~np.isnan(interpolated_bs)
magnetosheath_region = (interpolated_mp <= rad) & (rad <= interpolated_bs)
density = np.zeros_like(rad)
density[magnetosheath_region] = 1
x_sheath = X_grid[magnetosheath_region]
y_sheath = Y_grid[magnetosheath_region]
z_sheath = Z_grid[magnetosheath_region]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_sheath,y_sheath,z_sheath,alpha=0.5)
ax.plot_surface(x_mp, y_mp, z_mp, color='red', edgecolor='none', alpha=0.2, label="Magnetopause")
ax.plot_surface(x_bs, y_bs, z_bs, color='blue', edgecolor='none', alpha=0.2, label="Bow Shock")
ax.set_xlabel(r'$x$ ($R_{U}$)')
ax.set_ylabel(r'$y$ ($R_{U}$)')
ax.set_zlabel(r'$z$ ($R_{U}$)')

plt.show()