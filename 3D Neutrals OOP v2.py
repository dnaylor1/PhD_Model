import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Moon: #moon class to give properties to all the moons in the system
    def __init__(self, r_min, r_max, density, delta_Z, ER_min_H, ER_min_O, ER_max_H, ER_max_O):
        self.r_min = r_min
        self.r_max = r_max
        self.density = density
        self.delta_Z = delta_Z
        self.ER_min_H = ER_min_H #Eviator and Richardson neutrals, minimum and maximum for both H and O based neutrals (important for cross sections)
        self.ER_min_O = ER_min_O
        self.ER_max_H = ER_max_H
        self.ER_max_O = ER_max_O

    def add_density(self, rad, zbox, Z1):
        mask_r = (rad >= self.r_min) & (rad <= self.r_max) #radial mask
        mask_z = np.abs(zbox) < self.delta_Z #z-direction mask
        mask = mask_r & mask_z
        Z1[mask] = self.density 
        return Z1

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

################## plot a 2D projection too to make sure scale heights are working

Miranda = Moon(2.9, 9.9, 0.06, 0.38, 0.59, 0.57, 13.8, 11.68)
Ariel = Moon(3.8, 17, 0.11, 0.68, 4.37, 3.97, 28.4, 38.74)
Umbriel = Moon(4.7, 29, 0.12, 1.1, 4.06, 4.11, 15.5, 24.67)
Titania = Moon(6.1, 73, 0.02, 2.3, 22.3, 14.5, 38.1, 34.4)
Oberon = Moon(6.9, 150, 0.002, 3.6, 31, 19.5, 39.2, 30.6)
moons = [Miranda, Ariel, Umbriel, Titania, Oberon]

#system = System(moons,grid_limits=(-40,40,-40,40,-4,4))
system = System(moons)
total_density = system.calculate_total_density()
#system.plot_density()
