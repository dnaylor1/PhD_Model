import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

############## MAGNETOPAUSE

r0 = 16
K = 0.6
#theta = np.linspace(-np.pi + ((np.pi) / 20), 0, num=100)  
theta = np.linspace(-np.pi, 0, num=100)
exclude=np.pi
theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
theta = theta_exc
phi = np.linspace(0, 2 * np.pi, num=100)  
r = r0 * (2 / (1 + np.cos(theta))) ** K  
theta, phi = np.meshgrid(theta, phi)

### surface plot
x = r * np.cos(theta)
y = r * np.sin(theta) * np.cos(phi)
z = r * np.sin(theta) * np.sin(phi)

### sets limits in x, y and z so the surface doesn't extend past the axes limits
x_limit, y_limit, z_limit = 80,80,80
mask = (np.abs(x) > x_limit) | (np.abs(y) > y_limit) | (np.abs(z) > z_limit)
x = np.where(mask, np.nan, x)
y = np.where(mask, np.nan, y)
z = np.where(mask, np.nan, z)

ax.plot_surface(x, y, z, cmap='plasma', edgecolor='none',alpha=0.2)

############## BOW SHOCK

r0 = 20
K = 0.88
r = r0 * (2 / (1 + np.cos(theta))) ** K 
x = r * np.cos(theta)
y = r * np.sin(theta) * np.cos(phi)
z = r * np.sin(theta) * np.sin(phi)

x_limit, y_limit, z_limit = 80,80,80
mask = (np.abs(x) > x_limit) | (np.abs(y) > y_limit) | (np.abs(z) > z_limit)
x = np.where(mask, np.nan, x)
y = np.where(mask, np.nan, y)
z = np.where(mask, np.nan, z)

ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none',alpha=0.2)

### plots the lines to check the surface matches - first plot in the x-y plane
#theta = np.linspace(-np.pi + ((np.pi) / 20), 0, num=100)  
theta = theta_exc
x = r * np.cos(theta)
y = np.sqrt(r**2 - x**2)
z = np.zeros_like(x)  

### sets the limits of the lines (can't use clip with surface as it distorts it)
x = np.clip(x, -50, 50)
y = np.clip(y, -50, 50)
z = np.clip(z, -50, 50)

""" ax.plot(x, y, z, color='blue')
ax.plot(x, -y, z, color='blue') """

### then plots in the x-z plane
x = r * np.cos(theta)
z = np.sqrt(r**2 - x**2)
y = np.zeros_like(x)  

x = np.clip(x, -50, 50)
y = np.clip(y, -50, 50)
z = np.clip(z, -50, 50)

""" ax.plot(x, y, z, color='green')
ax.plot(x, y, -z, color='green') """

""" x_min = -10
x_max = 10
y_min = -10
y_max = 10
z_min = -5
z_max = 5
dx = 1
dy = 1
dz = 1
xvalues = np.arange(x_min,x_max+1,dx)
yvalues = np.arange(y_min,y_max+1,dy)
zvalues = np.arange(z_min,z_max+1,dz)
xbox, ybox, zbox = np.meshgrid(xvalues,yvalues,zvalues)
rad = np.sqrt(xbox**2 + ybox**2 + zbox**2)

############## CHENG NEUTRALS 

neutrals_rad = [(2.9,9.9),(3.8,17),(4.7,29),(6.1,73),(6.9,150)]
neutrals_n = [0.06,0.11,0.12,0.02,0.002]
delta_Z = [0.38,0.68,1.1,2.3,3.6]
r0 = [5.1,7.5,10.4,17.0,22.8]
maps = ['Blues','Reds','Greens','Purples','Oranges']

Z1 = np.zeros_like(rad)

for (r_min, r_max), density, scale_height in zip(neutrals_rad, neutrals_n, delta_Z):
    mask_r = (rad>=r_min) & (rad<=r_max)
    mask_z = np.abs(ybox) < scale_height
    mask = mask_r & mask_z
    Z1[mask] += density

ax.scatter(xbox,ybox,zbox,c=Z1,cmap='plasma',alpha=0.1) """

ax.set_xlabel(r'$x$ ($R_{U}$)')
ax.set_ylabel(r'$y$ ($R_{U}$)')
ax.set_zlabel(r'$z$ ($R_{U}$)')
""" ax.set_xlim([-50, 50])
ax.set_ylim([-50, 50])
ax.set_zlim([-50, 50]) """

#ax.view_init(elev=20, azim=-140)
#ax.view_init(elev=75, azim=-61)
#ax.view_init(elev=35,azim=-120) 
ax.view_init(elev=20, azim=-130)
    # 20 -140 for back view
    # 75 -61 for top view
    # 35 -120 for MP/BS view
    #20 -130 for neutrals
#plt.savefig('3D Magnetopause Attempt 4 Neutrals Rotated', dpi=1200)
plt.savefig('3D Magnetopause Attempt 4 Bow Shock Full', dpi=1200)

""" # Define the spherical coordinates
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_sphere = np.outer(np.cos(u), np.sin(v))
y_sphere = np.outer(np.sin(u), np.sin(v))
z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the sphere at the origin with radius 1
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='black', alpha=1, edgecolor='none') """



############## Work backwards so x is meshgrid and get r from that? Maybe possible but not sure 

plt.show()
