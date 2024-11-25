import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Import import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Moons import moons
from Plotter import *
import json


class System: #class for the model as a whole including setting up the grid and plotting
    def __init__(self, moons, grid_limits=(-10,10,-10,10,-10,10), grid_resolution=(1, 1, 1)): #sets these default values but they can be changed by passing through different values
        self.moons = moons #list of the Moons to iterate through later 
        self.dx, self.dy, self.dz = grid_resolution #sets dx,dy,dz so they can be accessed later on
        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = grid_limits #sets the grid limits so they can also be accessed later 
        self.xbox, self.ybox, self.zbox, self.rad, self.Z1, self.X_grid, self.Y_grid, self. Z_grid = self.create_grid() #calls the create grid function to assign these variables. Doesn't need to pass anything due to using self.---

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
        #xvalues = np.arange(self.x_min, self.x_max + self.dx, self.dx)
        #yvalues = np.arange(self.y_min, self.y_max + self.dy, self.dy)
        #zvalues = np.arange(self.z_min, self.z_max + self.dz, self.dz)
        xvalues = np.arange(self.x_min, self.x_max, self.dx)
        yvalues = np.arange(self.y_min, self.y_max, self.dy)
        zvalues = np.arange(self.z_min, self.z_max, self.dz)
        X_grid, Y_grid, Z_grid = np.mgrid[self.x_min:self.x_max:self.dx, self.y_min:self.y_max:self.dy, self.z_min:self.z_max:self.dz]
        xbox, ybox, zbox = np.meshgrid(xvalues, yvalues, zvalues)
        rad = np.sqrt(xbox**2 + ybox**2 + zbox**2)
        Z1 = np.zeros_like(rad)
        return xbox, ybox, zbox, rad, Z1, X_grid, Y_grid, Z_grid

    def calculate_total_density(self):
        """
        Calculates the total density at each point, i.e. iterates through each torus to calculate total. Then calls the plot method.
        """
        self.total_density = 0
        for moon in self.moons: #for each moon in the list of moons
            self.Z1 = moon.add_density(self.rad, self.zbox, self.Z1) #calls the add density method in the moon class to calculate the density at each point due to this moon. Adds to the total density
            self.total_density += self.Z1 #sets the total density to this variable due to have a better name (could be changed by renaming Z1)
            self.Z1 = np.zeros_like(self.rad) #resets Z1 so that when it is passed add_density again in the next loop, it doesn't pass a non-zero array.
        #return total_density #returns this to the bottom of the code, then passed back into the class when calling plot_density. Could be changed by using self.total_density maybe?
        #Plotter.plot_density(self, self.total_density) #once the total density at each point is calculated, the plot_density method is called. Saves having to have another line at the bottom of the code.
        return self.total_density

"""     def sheath(self, r_mp_grid, r_bs_grid, x_mp, y_mp, z_mp, x_bs, y_bs, z_bs):
        valid = ~np.isnan(r_mp_grid) & ~np.isnan(r_bs_grid)
        magnetosheath_region = valid & (r_mp_grid <= self.rad) & (self.rad <= r_bs_grid)
        density = np.zeros_like(self.rad)  
        density[magnetosheath_region] = 1  
        ##extract the magnetosheath points 
        x_sheath = self.xbox[magnetosheath_region]
        y_sheath = self.ybox[magnetosheath_region]
        z_sheath = self.zbox[magnetosheath_region]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_sheath,y_sheath,z_sheath)
        ax.plot_surface(x_mp, y_mp, z_mp, cmap='plasma', edgecolor='none', alpha=0.2)
        ax.plot_surface(x_bs, y_bs, z_bs, cmap='viridis', edgecolor='none', alpha=0.2)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        ax.set_zlabel(r'$z$ ($R_{U}$)')
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_zlim(-80,80) """

