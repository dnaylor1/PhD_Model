from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import UnivariateSpline
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Miscellaneous.Moons import moons
from Plotter import *
import json
from scipy import constants


class SXI:
    def __init__(self,X_grid,Y_grid,Z_grid,ver,grid_res,rad):
        self.R_U = 25362*(10**5) #Uranus radius in cm
        self.X_grid = X_grid
        self.Y_grid = Y_grid
        self.Z_grid = Z_grid
        self.ver = ver
        self.dx, self.dy, self.dz = grid_res
        self.rad = rad
    
    def flux(self):
        X_grid = self.X_grid.astype(np.float64)
        Y_grid = self.Y_grid.astype(np.float64)
        Z_grid = self.Z_grid.astype(np.float64)
        R_U = np.float64(self.R_U)
        distance_SMILE_Particle = np.sqrt(((300*R_U)-(X_grid*R_U))**2 + (Y_grid*R_U)**2 + (Z_grid*R_U)**2)
        ver_scaled = self.ver * (1/(distance_SMILE_Particle)**2)
        dx = self.dx * R_U
        dy = self.dy * R_U
        dz = self.dz * R_U
        intensity_integral = ver_scaled * dx * dy * dz
        x_projection = np.sum(intensity_integral,axis=0)
        flux = x_projection
        return flux
    
    def integration_time(self,flux):
        counts_second = flux*9.6
        counts_hour = counts_second * 3600
        integration_time_s = 1/(counts_second.max())
        integration_time_h = 1/(counts_hour.max())
        integration_time_sec_3sf = f"{integration_time_s:.3g}" #to give integration time to 3sf on graph
        integration_time_hour_3sf = f"{integration_time_h:.3g}"
        return integration_time_sec_3sf, integration_time_hour_3sf

    

""" intensity_integral_slow = intensity_scaled_slow * dl * dl * (4*dl) ## assume 2Rs scale height  ## *dl*dl*(40dl) = dxdydz // units of cm^-2
x_projection_slow = np.sum(intensity_integral_slow,axis=1) ##axis = 0 for column, 1 for row.
num_stacks=202
stackedx_projection_slow = [x_projection_slow for x in range(num_stacks)]
verticalstack_slow = np.vstack(stackedx_projection_slow)
verticalstackarea_slow = verticalstack_slow*9.6 #area of detector
perhour_slow = verticalstackarea_slow*3600 #for time

############## FLUX AND INTEGRATION TIMES

R_U = 25362*(10**5) #R_U in cm
distance_SMILEToParticle = np.sqrt(((300*R_U) - (xbox*R_U))**2 + (ybox*R_U)**2) #kept at ~300 R_U (similar to Saturn (but R_U not R_S so significantly closer))
intensity_scaled = intensity_values * (1/(distance_SMILEToParticle)**2) #scaled to spacecraft distance
dl = dx*R_U ## dx = 1Rs, converted to cm

intensity_integrals = np.zeros_like(rad)
for (r_inner, r_outer), scale_height in zip(neutrals_rad, delta_Z):
    integral = intensity_scaled * dl * dl * (scale_height*R_U)
    intensity_integrals += integral """

""" #### avoiding overlapping regions so the integral isn't calculated for points twice.
intensity_integrals = np.zeros_like(rad)
covered = np.zeros_like(rad, dtype=bool)
for (r_inner, r_outer), scale_height in zip(neutrals_rad, delta_Z):
    torus_mask = (rad >= r_inner) & (rad <= r_outer) & (~covered)
    integral = torus_mask * intensity_scaled * dl * dl * (scale_height * R_U)
    intensity_integrals += integral
    covered |= torus_mask """

""" total_intensity_integral = intensity_integrals
print(np.shape(total_intensity_integral))
x_projection = np.sum(total_intensity_integral,axis=1) ##axis = 0 for column, 1 for row. So sums each row
num_stacks=402 #same number of grid cells in xbox and ybox (0.5*201 (1 for the 0 point))
stackedx_projection = [x_projection for x in range(num_stacks)] #stacks the x-projection such that is is the same shape as the meshgrid. 
verticalstack = np.vstack(stackedx_projection) #vertically stacks the stackedx_projection, necessary for contouring
verticalstackarea = verticalstack*9.6 #area of detector
perhour = verticalstackarea*3600 #for time

counts_sec = verticalstackarea #counts per second
counts_hour = perhour #counts per hour
flux_detected = verticalstack #flux

integration_time_sec = 1/(counts_sec.max()) #PEAK integration time here
integration_time_hour = 1/(counts_hour.max())
integration_time_sec_3sf = f"{integration_time_sec:.3g}" #to give integration time to 3sf on graph
integration_time_hour_3sf = f"{integration_time_hour:.3g}"

##################################################################################################################
###### Unscaled counts to give flux with distance plot

unscaled_intensity_integrals = np.zeros_like(rad)
for (r_inner, r_outer), scale_height in zip(neutrals_rad, delta_Z):
    integral = intensity_values * dl * dl * (scale_height*R_U)
    unscaled_intensity_integrals += integral
total_unscaled_intensity_integral = unscaled_intensity_integrals
unscaled_projection = np.sum(total_unscaled_intensity_integral,axis=1)
num_stacks=402
unscaled_stackedproj = [unscaled_projection for x in range(num_stacks)]
unscaled_vertical = np.vstack(unscaled_stackedproj)
unscaled_verticalarea = unscaled_vertical*9.6
unscaled_verticalhour = unscaled_verticalarea*3600
unscaled_counts = unscaled_verticalarea
total_unscaled = np.sum(unscaled_counts,axis=0)
unscaled_sum = np.sum(total_unscaled)
distance = np.linspace(200,1000,800)
flux_scaling = (unscaled_sum)*1/((distance*R_U)**2)
 """