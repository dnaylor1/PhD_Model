from Import import *
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


class System: #class for the model as a whole including setting up the grid and plotting
    def __init__(self, moons, B_eq, grid_limits=(-10,10,-10,10,-10,10), grid_resolution=(1, 1, 1)): #sets these default values but they can be changed by passing through different values
        self.moons = moons #list of the Moons to iterate through later 
        self.dx, self.dy, self.dz = grid_resolution #sets dx,dy,dz so they can be accessed later on
        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = grid_limits #sets the grid limits so they can also be accessed later 
        #self.xbox, self.ybox, self.zbox, self.rad, self.Z1, self.X_grid, self.Y_grid, self. Z_grid = self.create_grid() #calls the create grid function to assign these variables. Doesn't need to pass anything due to using self.---
        self.rad, self.Z1, self.X_grid, self.Y_grid, self. Z_grid = self.create_grid()
        self.B_eq = B_eq

    def create_grid(self):
        """
        Creates the 3D meshgrid for the neutrals (MP/BS are currently on theta,phi meshgrid so created separately)

        Returns
            xbox (ndarray): x meshgrid
            ybox (ndarray): y meshgrid
            zbox (ndarray): z meshgrid
            rad (ndarray): radius of each point on the grid
            Z1 (ndarray): empty rad-like array to be updated with density contributions
        """
        X_grid, Y_grid, Z_grid = np.mgrid[self.x_min:self.x_max:self.dx, self.y_min:self.y_max:self.dy, self.z_min:self.z_max:self.dz]
        rad = np.sqrt(X_grid**2 + Y_grid**2 + Z_grid**2)
        Z1 = np.zeros_like(rad)
        return rad, Z1, X_grid, Y_grid, Z_grid

    def calculate_total_density(self,config):
        """
        Calculates the total density at each point, i.e. iterates through each torus to calculate total. Then calls the plot method.
        
        Returns
            total_density (ndarray): total density at each point (moon density + exospheric density).
            moon_density (ndarray): moon density at each point.
            n_exo (ndarray): exospheric density at each point.
        """
        self.moon_density = np.zeros_like(self.rad)
        self.total_density = np.zeros_like(self.rad)
        n_exo = self.exosphere()
        self.total_density += n_exo
        for moon in self.moons: #for each moon in the list of moons
            self.Z1 = moon.add_density(self.rad, self.X_grid, self.Y_grid, self.Z_grid, self.Z1, config) #calls the add density method in the moon class to calculate the density at each point due to this moon. Adds to the total density
            self.moon_density += self.Z1 #sets the total density to this variable due to have a better name (could be changed by renaming Z1)
            self.Z1 = np.zeros_like(self.rad) #resets Z1 so that when it is passed add_density again in the next loop, it doesn't pass a non-zero array.
        #return total_density #returns this to the bottom of the code, then passed back into the class when calling plot_density. Could be changed by using self.total_density maybe?
        #Plotter.plot_density(self, self.total_density) #once the total density at each point is calculated, the plot_density method is called. Saves having to have another line at the bottom of the code.
        self.total_density += self.moon_density
        return self.total_density, self.moon_density, n_exo
    
    def exosphere(self):
        """
        Calculates the density of the exosphere at each point, returns to the calculate total density method

        Returns:
            n_exo (ndarray): exospheric density of each point in the grid
        """
        c_1 = 4.2e-5
        c_2 = 31
        n_exo = np.zeros_like(self.rad)
        mask = self.rad > 1
        n_exo[mask] = c_1 * np.exp((c_2)/self.rad[mask])
        return n_exo

    


        
    




