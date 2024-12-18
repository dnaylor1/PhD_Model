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
    def __init__(self,X_grid,Y_grid,Z_grid,grid_res,rad,FOV,Aeff):
        self.R_U = 25362*(10**5) #Uranus radius in cm
        self.X_grid = X_grid
        self.Y_grid = Y_grid
        self.Z_grid = Z_grid
        self.dx, self.dy, self.dz = grid_res
        self.rad = rad
        self.FOV = FOV
        self.Aeff = Aeff

    def SXI_distance(self):
        self.distance = 80/(np.tan(self.FOV/2)) - 60
    
    def flux(self,ver):
        X_grid = self.X_grid.astype(np.float64)
        Y_grid = self.Y_grid.astype(np.float64)
        Z_grid = self.Z_grid.astype(np.float64)
        R_U = np.float64(self.R_U)
        distance_SMILE_Particle = np.sqrt(((self.distance*R_U)-(X_grid*R_U))**2 + (Y_grid*R_U)**2 + (Z_grid*R_U)**2)
        ver_scaled = ver * (1/(distance_SMILE_Particle)**2)
        dx = self.dx * R_U
        dy = self.dy * R_U
        dz = self.dz * R_U
        intensity_integral = ver_scaled * dx * dy * dz
        x_projection = np.sum(intensity_integral,axis=0)
        flux = x_projection
        return flux
    
    def integration_time(self,flux):
        counts_second = flux*self.Aeff
        counts_hour = counts_second * 3600
        integration_time_s = 1/(counts_second.max())
        integration_time_h = 1/(counts_hour.max())
        integration_time_sec_3sf = f"{integration_time_s:.3g}" #to give integration time to 3sf on graph
        integration_time_hour_3sf = f"{integration_time_h:.3g}"
        return integration_time_sec_3sf, integration_time_hour_3sf

    

