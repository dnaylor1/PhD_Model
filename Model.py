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
combd = None

###### SOLAR WIND VARIATIONS

### NOVEMBER 7-8 2023 HIGH SW SPEED
""" n_p = 24.3e6
v_sw = 370e3
T_sw = 4.16e4
 """

### MAY 28 2008, 01:00-03:00 HIGH SW DENSITY
""" n_p = 1.856e6
v_sw = 690e3
T_sw = 2.54e5 """

### Slow-Fast Distinction
combd = "Y"
v_slow = 400e3
v_fast = 800e3


###### SET GRID LIMITS AND RESOLUTION
#grid_lims = [-40,40,-40,40,-40,40]
grid_res = [1,1,1]
grid_lims=[-80,80,-80,80,-80,80]

###### MOONS
miranda = Miranda()
ariel = Ariel()
umbriel = Umbriel()
titania = Titania()
oberon = Oberon()
moons = [miranda,ariel,umbriel,titania,oberon]

###### INITIALISE SYSTEM AND PLOTTER CLASS
system = System(moons,B_eq=2.3e-5,grid_limits=grid_lims,grid_resolution=grid_res)
plotter = Plotter(grid_lims,system.X_grid,system.Y_grid,system.Z_grid)

###### SYSTEM NEUTRAL DENSITIES
total_density, moon_density, n_exo = system.calculate_total_density()
#plotter.plot_density(moon_density)
#plotter.plot_exo(n_exo)

###### MAGNETOPAUSE AND BOW SHOCK SURFACES - singular
if combd == None:
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

else:
    magnetopause = Surface(r0=16,K=0.6)
    bow_shock = Surface(r0=20,K=0.88)
    
    x_mp_s,y_mp_s,z_mp_s,r0_mag_s,k_mag_s = magnetopause.define_surface(system.X_grid,n_p,v_sw,"MP",grid_lims[0],grid_lims[1],combd="Y",v_type="S")
    x_bs_s,y_bs_s,z_bs_s,r0_bow_s,k_bow_s = bow_shock.define_surface(system.X_grid,n_p,v_sw,"BS",grid_lims[2],grid_lims[3],combd="Y",v_type="S")

    x_mp_f,y_mp_f,z_mp_f,r0_mag_f,k_mag_f = magnetopause.define_surface(system.X_grid,n_p,v_sw,"MP",grid_lims[0],grid_lims[1],combd="Y",v_type="F")
    x_bs_f,y_bs_f,z_bs_f,r0_bow_f,k_bow_f = bow_shock.define_surface(system.X_grid,n_p,v_sw,"BS",grid_lims[2],grid_lims[3],combd="Y",v_type="F")

    ##### MAGNETOSHEATH
    sheath_slow = Magnetosheath(x_mp_s,y_mp_s,z_mp_s,x_bs_s,y_bs_s,z_bs_s,system.X_grid,system.Y_grid,system.Z_grid,system.rad)
    sheath_density_slow = sheath_slow.sheath_surface()
    sheath_fast = Magnetosheath(x_mp_f,y_mp_f,z_mp_f,x_bs_f,y_bs_f,z_bs_f,system.X_grid,system.Y_grid,system.Z_grid,system.rad)
    sheath_density_fast = sheath_slow.sheath_surface()

    ###### VOLUMETRIC EMISSION
    ver_slow = system.volumetric_emission(total_density,sheath_density_slow,n_p,T_sw,v_sw)
    ver_fast = system.volumetric_emission(total_density,sheath_density_fast,n_p,T_sw,v_sw)
    
    
    
    
    
    
plt.show()

