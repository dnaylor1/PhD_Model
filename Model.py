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

miranda = Miranda()
ariel = Ariel()
umbriel = Umbriel()
titania = Titania()
oberon = Oberon()
moons =[miranda,ariel,umbriel,titania,oberon]
system = System(moons,grid_limits=(-80,80,-80,80,-80,80))
#system = System(moons)
magnetopause = Surface(r0=16,K=0.6)
bow_shock = Surface(r0=20,K=0.88)
plotter = Plotter(-80,80,-80,80,-80,80,system.xbox,system.ybox,system.zbox,system.X_grid,system.Y_grid,system.Z_grid)
#plotter = Plotter(-10,10,-10,10,-10,10,system.xbox,system.ybox,system.zbox,system.X_grid,system.Y_grid,system.Z_grid)
total_density = system.calculate_total_density()
#plotter.plot_density(total_density)
x_mp,y_mp,z_mp = magnetopause.define_surface_NEW(system.X_grid)
x_bs,y_bs,z_bs = bow_shock.define_surface_NEW(system.X_grid)
plotter.plot_surfaces(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)
#r_mp_grid = Surface.interpolate(x_mp,y_mp,z_mp)
#r_bs_grid = Surface.interpolate(x_bs,y_bs,z_bs)
#system.sheath(r_mp_grid,r_bs_grid,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)

plt.show()

