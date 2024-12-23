from Import import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Plotter import *
from SXI import *
import json
from Magnetosheath import *
from PIL import Image
import io

#Solstice or equiniox ("S" or "E") #E default
config = "E"
config_combd = None

n_p = None #defaults when no solar wind variations used
v_sw = None
T_sw = None
j_s = None
r0_mp = None
K = None

###### SOLAR WIND VARIATIONS

""" ### Jasinski et al. (2024)
K = None
j_s = True """

""" ## 1985, DOY 290.16667
v_sw = 435.2e3     
n_p = 0.02139e6    
P_dyn = 0.0073806539 #nPa      
r0_mp = 20.319744 """


""" ## DOY: 350.58333
v_sw= 480.4e3    
n_p = 0.0030400000e6    
P_dyn = 0.0012781619       
r0_mp = 27.638616 """

### NOVEMBER 7-8 2023 HIGH SW SPEED
""" n_p = 24.3e6
v_sw = 370e3
T_sw = 4.16e4 """

### MAY 28 2008, 01:00-03:00 HIGH SW DENSITY
""" n_p = 1.856e6
v_sw = 690e3
T_sw = 2.54e5 """

### Slow-Fast Distinction - uncomment for combined figure
combd = None
v_slow = None
v_fast = None

""" combd = "Y"
v_slow = 400e3
v_fast = 800e3 """

###### SET GRID LIMITS AND RESOLUTION

""" grid_lims = [-40,40,-40,40,-40,40]
grid_res = [2,2,2] """

grid_res = [2,2,2]
grid_lims=[-80,80,-80,80,-80,80]

""" grid_res = [0.05,0.05,0.05]
grid_lims = [-10,10,-10,10,-10,10] """

###### MOONS
miranda = Miranda()
ariel = Ariel()
umbriel = Umbriel()
titania = Titania()
oberon = Oberon()
moons = [miranda,ariel,umbriel,titania,oberon]

###### INITIALISE SYSTEM AND PLOTTER CLASS
system = System(moons,B_eq=2.3e-5,grid_limits=grid_lims,grid_resolution=grid_res)
plotter = Plotter(grid_lims,system.X_grid,system.Y_grid,system.Z_grid,config)

###### SYSTEM NEUTRAL DENSITIES
total_density, moon_density, n_exo = system.calculate_total_density(config)
#plotter.plot_density(moon_density)
#plotter.plot_exo(n_exo)

###### MAGNETOPAUSE AND BOW SHOCK SURFACES - singular
if combd == None and config_combd == None:
    if j_s == True:
        magnetopause = Surface(r0=r0_mp,K=K)
        bow_shock = Surface(r0=r0_mp,K=K)
    else:
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
    ver = sheath.volumetric_emission(total_density,sheath_density,n_p,T_sw,v_sw,j_s=j_s)
    #plotter.plot_ver(ver,r0_mag,k_mag,r0_bow,k_bow)

    ###### FLUX
    SMILE = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,26.5,9.6)
    SMILE.SXI_distance()
    flux = SMILE.flux(ver)
    int_s,int_h = SMILE.integration_time(flux)
    #plotter.plot_all_ver_flux(ver,flux,int_s,int_h,x_pos,y_pos,z_pos,v_sw)
    #plotter.plot_flux_gif(ver,flux,int_s,int_h,config)
    #plotter.plot_ver_separate(ver)
    #plotter.plot_distance_counts(ver,dx=1,A_eff=9.6)

elif config_combd != None:
    if j_s == True:
        magnetopause = Surface(r0=r0_mp,K=K)
        bow_shock = Surface(r0=r0_mp,K=K)
    else:
        magnetopause = Surface(r0=16,K=0.6)
        bow_shock = Surface(r0=20,K=0.88)
    total_density_s, moon_density_s, n_exo_s = system.calculate_total_density(config="S")
    total_density_e, moon_density_e, n_exo_e = system.calculate_total_density(config="E")
    #plotter.plot_density_combined(moon_density_s,moon_density_e)
    x_mp,y_mp,z_mp,r0_mag,k_mag = magnetopause.define_surface(system.X_grid,n_p,v_sw,"MP",grid_lims[0],grid_lims[1])
    x_bs,y_bs,z_bs,r0_bow,k_bow = bow_shock.define_surface(system.X_grid,n_p,v_sw,"BS",grid_lims[2],grid_lims[3])
    sheath = Magnetosheath(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,system.X_grid,system.Y_grid,system.Z_grid,system.rad)
    sheath_density = sheath.sheath_surface()
    ver_equinox = sheath.volumetric_emission(total_density_e,sheath_density,n_p,T_sw,v_sw,j_s=j_s)
    ver_solstice = sheath.volumetric_emission(total_density_s,sheath_density,n_p,T_sw,v_sw,j_s=j_s)
    SMILE = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,26.5,9.6)
    SMILE.SXI_distance()
    #LEXI = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,9.1,44.18)
    #LEXI.SXI_distance()
    #future = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,53,19.2)
    #future.SXI_distance()
    flux_equinox = SMILE.flux(ver_equinox)
    int_s_equ,int_h_equ = SMILE.integration_time(flux_equinox)
    flux_solstice = SMILE.flux(ver_solstice)
    int_s_sol,int_h_sol = SMILE.integration_time(flux_solstice)
    plotter.plot_ver_flux_equsol2(ver_equinox,ver_solstice,flux_equinox,flux_solstice)

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
    sheath_density_fast = sheath_fast.sheath_surface()

    ###### VOLUMETRIC EMISSION
    ver_slow = sheath_slow.volumetric_emission(total_density,sheath_density_slow,n_p,T_sw,v_sw,v_slow)
    ver_fast = sheath_fast.volumetric_emission(total_density,sheath_density_fast,n_p,T_sw,v_sw,v_fast)
    SMILE_slow = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,FOV=26.5,Aeff=9.6)
    SMILE_fast = SXI(system.X_grid,system.Y_grid,system.Z_grid,grid_res,system.rad,FOV=26.5,Aeff=9.6)
    SMILE_slow.SXI_distance()
    SMILE_fast.SXI_distance()
    flux_slow = SMILE_slow.flux(ver_slow)
    flux_fast = SMILE_fast.flux(ver_fast)
    int_s_s, int_h_s = SMILE_slow.integration_time(flux_slow)
    int_s_f, int_h_f = SMILE_fast.integration_time(flux_fast)

plt.show()

