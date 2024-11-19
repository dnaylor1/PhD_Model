import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Moons import moons
import json

## pip freeze > file_name to set requirements (virtual environment)
## always use virtual environment - keeps versions the same so the code won't be broken using extensions
## open model folder directly rather than single file
## use launch json file to always launch model rather than the other files
## can use a config json file too to have all the numbers in
## keep main model file outside of the PhD_Work folder as it causes problems with git
## breakpoints- using F10 and F11 to step through code
## good practice to define data types for all variables
## also do the """ """ commenting for all methods

################## plot a 2D projection too to make sure scale heights are working

system = System(moons,grid_limits=(-80,80,-80,80,-80,80))
#system = System(moons)
magnetopause = Surface(r0=16,K=0.6)
bow_shock = Surface(r0=20,K=0.88)
total_density = system.calculate_total_density()
x_mp,y_mp,z_mp = magnetopause.define_surface()
x_bs,y_bs,z_bs = bow_shock.define_surface()
system.plot_surfaces(x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)
r_mp_grid = Surface.interpolate(x_mp,y_mp,z_mp)
r_bs_grid = Surface.interpolate(x_bs,y_bs,z_bs)
system.sheath(r_mp_grid,r_bs_grid,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs)

plt.show()


## SW solar wind variations - change flaring
