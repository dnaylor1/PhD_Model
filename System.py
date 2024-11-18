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
        """
        Creates the 3D meshgrid for the neutrals (MP/BS are currently on theta,phi meshgrid so created separately)

        Returns
            xbox (ndarray): x meshgrid
            ybox (ndarray): y meshgrid
            zbox (ndarray): z meshgrid
            rad (ndarray): radius of each point on the grid
            Z1 (ndarray): empty rad-like array to be updated with density contributions
        """
        xvalues = np.arange(self.x_min, self.x_max + self.dx, self.dx)
        yvalues = np.arange(self.y_min, self.y_max + self.dy, self.dy)
        zvalues = np.arange(self.z_min, self.z_max + self.dz, self.dz)
        xbox, ybox, zbox = np.meshgrid(xvalues, yvalues, zvalues)
        rad = np.sqrt(xbox**2 + ybox**2 + zbox**2)
        Z1 = np.zeros_like(rad)
        return xbox, ybox, zbox, rad, Z1

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
        self.plot_density() #once the total density at each point is calculated, the plot_density method is called. Saves having to have another line at the bottom of the code.

    def plot_density(self):
        """
        Plots the total density at each point on a 3D grid.
        """
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
        ax.view_init(elev=10,azim=60)
        #ax.tick_params(axis='both', direction='out', top=True, bottom=True, left=True, right=True)
        ax.tick_params(axis='both', which='major', labelsize=8)
        #plt.savefig("3D Neutrals Small Labels",dpi=1200)

    def plot_surfaces(self, x_mp,y_mp,z_mp,x_bs,y_bs,z_bs):
        """
        Plots the magnetopause and bow shock surfaces.

        Parameters
            x_mp (ndarray): magnetopause x-coordinates
            y_mp (ndarray): magnetopause y-coordinates
            z_mp (ndarray): magnetopause z-coordinates
            x_bs (ndarray): bow shock x-coordinates
            y_bs (ndarray): bow shock y-coordinates
            z_bs (ndarray): bow shock z-coordinates
        """
        x_limit, y_limit, z_limit = 80,80,80
        mask = (np.abs(x_mp) > x_limit) | (np.abs(y_mp) > y_limit) | (np.abs(z_mp) > z_limit)
        mask2 = (np.abs(x_bs) > x_limit) | (np.abs(y_bs) > y_limit) | (np.abs(z_bs) > z_limit)
        x_mp = np.where(mask, np.nan, x_mp)
        y_mp = np.where(mask, np.nan, y_mp)
        z_mp = np.where(mask, np.nan, z_mp)
        x_bs = np.where(mask2, np.nan, x_bs)
        y_bs = np.where(mask2, np.nan, y_bs)
        z_bs = np.where(mask2, np.nan, z_bs)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.plot_surface(x_mp,y_mp,z_mp,cmap='plasma', edgecolor='none', alpha=0.2)
        #ax.plot_surface(x_bs,y_bs,z_bs,cmap='viridis', edgecolor='none',alpha=0.2)

        ax.plot_surface(x_mp,y_mp,z_mp,cmap='plasma', edgecolor='none', alpha=0.2)
        ax.plot_surface(x_bs,y_bs,z_bs,cmap='viridis', edgecolor='none',alpha=0.2)

        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        ax.set_zlabel(r'$z$ ($R_{U}$)')
        ax.view_init(elev=20, azim=-130)
        #plt.savefig("3D MP BS",dpi=1200)

