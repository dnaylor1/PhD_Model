import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class System: #class for the model as a whole including setting up the grid and plotting
    def __init__(self, moons, grid_limits=(-10,10,-10,10,-5,5), grid_resolution=(1, 1, 1)): #sets these default values but they can be changed by passing through different values
        self.moons = moons #list of the Moons to iterate through later 
        self.dx, self.dy, self.dz = grid_resolution #sets dx,dy,dz so they can be accessed later on
        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = grid_limits #sets the grid limits so they can also be accessed later 
        self.xbox, self.ybox, self.zbox, self.rad, self.Z1 = self.create_grid() #calls the create grid function to assign these variables. Doesn't need to pass anything due to using self.---

    def create_grid(self):
        xvalues = np.arange(self.x_min, self.x_max + self.dx, self.dx)
        yvalues = np.arange(self.y_min, self.y_max + self.dy, self.dy)
        zvalues = np.arange(self.z_min, self.z_max + self.dz, self.dz)
        xbox, ybox, zbox = np.meshgrid(xvalues, yvalues, zvalues)
        rad = np.sqrt(xbox**2 + ybox**2 + zbox**2)
        Z1 = np.zeros_like(rad)
        return xbox, ybox, zbox, rad, Z1

    def calculate_total_density(self):
        self.total_density = 0
        for moon in self.moons: #for each moon in the list of moons
            self.Z1 = moon.add_density(self.rad, self.zbox, self.Z1) #calls the add density method in the moon class to calculate the density at each point due to this moon. Adds to the total density
            self.total_density += self.Z1 #sets the total density to this variable due to have a better name (could be changed by renaming Z1)
            self.Z1 = np.zeros_like(self.rad) #resets Z1 so that when it is passed add_density again in the next loop, it doesn't pass a non-zero array.
        #return total_density #returns this to the bottom of the code, then passed back into the class when calling plot_density. Could be changed by using self.total_density maybe?
        self.plot_density() #once the total density at each point is calculated, the plot_density method is called. Saves having to have another line at the bottom of the code.

    def plot_density(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        im = ax.scatter(self.xbox, self.ybox, self.zbox, c=self.total_density, cmap='plasma', alpha=0.4)
        fig.colorbar(im, shrink=0.5)
        ax.set_xlabel(r"$x$ ($R_{U}$)")
        ax.set_ylabel(r"$y$ ($R_{U}$)")
        ax.set_zlabel(r"$z$ ($R_{U}$)")
        ax.set_xlim([self.x_min, self.x_max])
        ax.set_ylim([self.y_min, self.y_max])
        ax.set_zlim([self.z_min, self.z_max])
        plt.tight_layout()
        plt.show()
