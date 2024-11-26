from Import import *
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

class Surface:
    def __init__(self, r0, K):
        self.r0 = r0
        self.K = K
    
    def define_surface(self):
        """
        Defines the surface for either the magnetopause or bow shock.
        """
        theta = np.linspace(-np.pi,0,162)
        exclude = np.pi
        theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
        theta = theta_exc
        phi = np.linspace(0,2*np.pi,161)
        r = self.r0 * (2 / (1 + np.cos(theta)))**self.K
        theta, phi = np.meshgrid(theta, phi)
        x = r*np.cos(theta)
        y = r*np.sin(theta)*np.cos(phi)
        z = r*np.sin(theta)*np.sin(phi)
        return x,y,z 
    
    def define_surface_NEW(self, X_grid):
        """
        Defines the surface for either the magnetopause or bow shock
        """
        theta = np.linspace(-np.pi,0,162)
        exclude = np.pi
        theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
        theta = theta_exc
        phi = np.linspace(0,2*np.pi,161)
        Phi, X = np.meshgrid(phi, X_grid[:,0,0])
        r_mp = self.r0 * (2 / (1 + np.cos(theta)))**self.K
        xmp = r_mp*np.cos(theta)
        interpfunc = interpolate.InterpolatedUnivariateSpline(xmp, r_mp)
        rmp = interpfunc(X)
        r_open = np.sqrt(rmp**2 - X**2)
        _y = r_open*np.sin(Phi)
        _z = r_open*np.cos(Phi)        
        x = X
        y = _y
        z = _z
        return x,y,z
    
    
    def sheath(x,y,z):
        x_flat = x.ravel()
        y_flat = y.ravel()
        z_flat = z.ravel()
        valid_points = ~np.isnan(x_flat) & ~np.isnan(y_flat) & ~np.isnan(z_flat) ## gets rid of nan values
        points = np.column_stack((x_flat[valid_points], y_flat[valid_points], z_flat[valid_points])) ##points is the coordinates of each point, input to interpolator to say where the surface is defined
        ## column stack used above because it combines the flattened surface arrays into a single 2D array where each row represents a 3D point. Required for interpolation.
        values_flat = np.sqrt(x**2 + y**2 + z**2).ravel() ##values is radial distance, output of interpolation, provide the MP/BS radius at each point.
        values = values_flat[valid_points] ##data cleaning: again filters out nan values.
        interp = LinearNDInterpolator(points,values)
        r_grid = interp(x,y,z)
        return r_grid
    