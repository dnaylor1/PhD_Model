from Import import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Plotter import *
import json
from Magnetosheath import *

n_p = None #defaults when no solar wind variations used
v_sw = None
T_sw = None

###### SOLAR WIND VARIATIONS

""" ### NOVEMBER 7-8 2023 HIGH SW SPEED
n_p = 24.3e6
v_sw = 370e3
T_sw = 4.16e4 """


""" ### MAY 28 2008, 01:00-03:00 HIGH SW DENSITY
n_p = 1.856e6
v_sw = 690e3
T_sw = 2.54e5 """

###### SET GRID LIMITS AND RESOLUTION
grid_lims = [-40,40,-40,40,-40,40]
grid_res = [1,1,1]

###### MOONS
miranda = Miranda()
ariel = Ariel()
umbriel = Umbriel()
titania = Titania()
oberon = Oberon()
moons = [miranda,ariel,umbriel,titania,oberon]

###### INITIALISE SYSTEM AND PLOTTER CLASS
system = System(moons,B_eq=2e-5,grid_limits=grid_lims,grid_resolution=grid_res)
plotter = Plotter(grid_lims,system.X_grid,system.Y_grid,system.Z_grid)

###### SYSTEM NEUTRAL DENSITIES
total_density, moon_density, n_exo = system.calculate_total_density()
#plotter.plot_density(moon_density)
#plotter.plot_exo(n_exo)

###### MAGNETOPAUSE AND BOW SHOCK SURFACES
magnetopause = Surface(r0=16,K=0.6)
bow_shock = Surface(r0=20,K=0.88)
x_mp,y_mp,z_mp,r0_mag,k_mag = magnetopause.define_surface(system.X_grid,n_p,v_sw,"MP",grid_lims[0],grid_lims[1])
x_bs,y_bs,z_bs,r0_bow,k_bow = bow_shock.define_surface(system.X_grid,n_p,v_sw,"BS",grid_lims[2],grid_lims[3])
#plotter.plot_surfaces(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)

##### MAGNETOSHEATH
sheath = Magnetosheath(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,system.X_grid,system.Y_grid,system.Z_grid,system.rad)
sheath_density = sheath.sheath_surface()
#plotter.plot_sheath(sheath_density,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,r0_mag,k_mag,r0_bow,k_bow)

###### VOLUMETRIC EMISSION
ver = system.volumetric_emission(total_density,sheath_density,n_p,T_sw,v_sw)
plotter.plot_ver(ver,r0_mag,k_mag,r0_bow,k_bow)

plt.show()

