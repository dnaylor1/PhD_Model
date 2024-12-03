from Import import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
#from Moons import moons
from Plotter import *
import json
from Magnetosheath import *

""" ### MAY 28 2008, 01:00-03:00 HIGH SW DENSITY
n_p = 24.3E6
v_sw = 370E3
T = 4.16E4 """

n_p = None
v_sw = None
T_sw = None

### Solar wind variations
n_p = 24.3e6
v_sw = 370e3
T_sw = 4.16e4
###

miranda = Miranda()
ariel = Ariel()
umbriel = Umbriel()
titania = Titania()
oberon = Oberon()
moons = [miranda,ariel,umbriel,titania,oberon]
#system = System(moons,grid_limits=(-80,80,-80,80,-80,80),grid_resolution=(0.5,0.5,0.5))
#system = System(moons,grid_limits=(-80,80,-80,80,-80,80),grid_resolution=(2,2,2))
system = System(moons,2.3e-5,grid_limits=(-80,80,-80,80,-80,80))
#system = System(moons,grid_limits=(-40,20,-80,80,40,40),grid_resolution=(0.5,0.5,0.5))
#system = System(moons,grid_limits=(-40,40,-40,40,-40,40))
#system = System(moons,2.3e-5,grid_limits=(-10,10,-10,10,-10,10),grid_resolution=(0.05,0.05,0.05))
#system = System(moons)
magnetopause = Surface(r0=16,K=0.6)
bow_shock = Surface(r0=20,K=0.88)
plotter = Plotter(-80,80,-80,80,-80,80,system.X_grid,system.Y_grid,system.Z_grid)
#plotter = Plotter(-40,40,-40,40,-40,40,system.xbox,system.ybox,system.zbox,system.X_grid,system.Y_grid,system.Z_grid)
#plotter = Plotter(-10,10,-10,10,-10,10,system.X_grid,system.Y_grid,system.Z_grid)
total_density, moon_density, n_exo = system.calculate_total_density()
#plotter.plot_density(moon_density)
#plotter.plot_exo(n_exo)
x_mp,y_mp,z_mp,r0_mag,k_mag = magnetopause.define_surface(system.X_grid,n_p,v_sw,type="MP")
x_bs,y_bs,z_bs,r0_bow,k_bow = bow_shock.define_surface(system.X_grid,n_p,v_sw,type="BS")
#plotter.plot_surfaces(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)
#r_mp_grid = Surface.interpolate(x_mp,y_mp,z_mp)
#r_bs_grid = Surface.interpolate(x_bs,y_bs,z_bs)
#system.sheath(r_mp_grid,r_bs_grid,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)
sheath = Magnetosheath(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,system.X_grid,system.Y_grid,system.Z_grid,system.rad)
points_mp,points_bs,x_points,y_points,z_points,x_sheath,y_sheath,z_sheath,density_grid = sheath.sheath_surface()
#plotter.plot_sheath(points_bs,points_mp,x_points,y_points,z_points,x_sheath,y_sheath,z_sheath,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,density_grid)
#plotter.plot_sheath(density_grid,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,r0_mag,k_mag,r0_bow,k_bow)
ver = system.volumetric_emission(total_density,density_grid,n_p,T_sw,v_sw)
plotter.plot_ver(ver,r0_mag,k_mag,r0_bow,k_bow)


plt.show()

