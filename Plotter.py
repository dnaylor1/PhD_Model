from Import import *
#import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *
from Surface import *
from Miscellaneous.Moons import moons
from Plotter import *
from Magnetosheath_New import *
import json

class Plotter:
    def __init__(self,x_min,x_max,y_min,y_max,z_min,z_max,xbox,ybox,zbox,X_grid,Y_grid,Z_grid):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max
        self.xbox, self.ybox, self.zbox = xbox,ybox,zbox
        self.X_grid, self.Y_grid, self.Z_grid = X_grid, Y_grid, Z_grid
    def plot_density(self,total_density):
        """
        Plots the total density at each point on a 3D grid.

        Parameters:
            total_density (ndarray): total density at each point in the grid from the moon tori and exosphere
        """
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        im = ax.scatter(self.xbox, self.ybox, self.zbox, c=total_density, cmap='plasma', alpha=0.4)
        fig.colorbar(im, shrink=0.5)
        ax.set_xlabel(r"$x$ ($R_{U}$)")
        ax.set_ylabel(r"$y$ ($R_{U}$)")
        ax.set_zlabel(r"$z$ ($R_{U}$)")
        #ax.set_xlim([self.x_min, self.x_max])
        #ax.set_ylim([self.y_min, self.y_max])
        #ax.set_zlim([self.z_min, self.z_max])
        ax.set_xlim(-10,10)
        ax.set_ylim(-10,10)
        ax.set_zlim(-10,10)
        plt.tight_layout()
        ax.view_init(elev=10,azim=60)
        #ax.tick_params(axis='both', direction='out', top=True, bottom=True, left=True, right=True)
        plt.title("3D Neutral Tori")
        ax.tick_params(axis='both', which='major', labelsize=8)
        #plt.savefig("3D Neutrals Small Labels",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        xmid = int(np.shape(self.Y_grid)[0]/2)
        yzplane = ax.contourf(self.Y_grid[xmid,:,:],self.Z_grid[xmid,:,:],total_density[xmid,:,:],cmap='plasma')
        fig.colorbar(yzplane,shrink=0.5,aspect=5)
        plt.title('Neutral Tori in y-z plane')
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("Neutrals y-z Projection",dpi=1200)

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
        plt.title("3D Magnetopause and Bow Shock Surfaces")
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        ax.set_zlabel(r'$z$ ($R_{U}$)')
        ax.view_init(elev=20, azim=-130)
        #plt.savefig("3D MP BS",dpi=1200)


    def plot_sheath(self,density_grid):
        """
        Plot the magnetosheath region

        Parameters:
            density_grid (ndarray): needs renaming*. Grid that defines whether points are inside or outside the magnetosheath. 
                                    Density set to 1 in the sheath and 0 elsewhere.
        """
        z_slice = 0
        # Indices where the grid intersects the slicing plane
        slice_indices = np.abs(self.Z_grid - z_slice) < 1e-2 
        ## grid points may not align exactly with the slicing plane due to interpolation so this makes sure they are close enough
        # Data points for the slicing plane
        x_slice = self.X_grid[slice_indices]
        y_slice = self.Y_grid[slice_indices]
        density_slice = density_grid[slice_indices]

        fig, ax = plt.subplots()
        scatter = ax.scatter(x_slice, y_slice, c=density_slice, cmap='Blues')

        from Surface import Surface
        bs = Surface(20,0.88)
        mp = Surface(16,0.6)
        x_bow, y_bow = bs.surf_2D()
        x_mag, y_mag = mp.surf_2D()
        
        ax.plot(x_bow,y_bow,color='red')
        ax.plot(x_mag,y_mag,color='blue')
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        fig.colorbar(scatter, label='Density')
        #plt.title(f'Slice through region at z = {z_slice}')
        plt.title(r"Magnetosheath $x$-$y$ Projection")
        plt.savefig("High res sheath x-y",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        x_mid = int(np.shape(self.X_grid)[0]/2)
        yz_plane = ax.contourf(self.Y_grid[x_mid,:,:],self.Z_grid[x_mid,:,:],density_grid[x_mid,:,:],cmap='Blues')
        fig.colorbar(yz_plane)
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        plt.title(r"Magnetosheath $y$-$z$ Projection")
        #plt.savefig("Magnetosheath y-z Projection Low Res",dpi=1200)
        plt.savefig("High res sheath y-z",dpi=1200)

        