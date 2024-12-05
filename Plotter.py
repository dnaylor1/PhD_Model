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
from Magnetosheath import *
import json

class Plotter:
    def __init__(self,grid_limits=[-10,10,-10,10,-10,10],X_grid=None,Y_grid=None,Z_grid=None):
        self.x_min = grid_limits[0]
        self.x_max = grid_limits[1]
        self.y_min = grid_limits[2]
        self.y_max = grid_limits[3]
        self.z_min = grid_limits[4]
        self.z_max = grid_limits[5]
        #self.xbox, self.ybox, self.zbox = xbox,ybox,zbox
        self.X_grid, self.Y_grid, self.Z_grid = X_grid, Y_grid, Z_grid
    def plot_density(self,moon_density):
        """
        Plots the total moon-sourced density at each point on a 3D grid.

        Parameters:
            moon_density (ndarray): total moon-sourced density at each point in the grid from the moon tori and exosphere
        """
        #fig = plt.figure()
        #ax = fig.add_subplot(projection='3d')
        #im = ax.scatter(self.X_grid,self.Y_grid,self.Z_grid,c=moon_density,cmap='plasma',alpha=0.4)
        #fig.colorbar(im, shrink=0.5)
        #ax.set_xlabel(r"$x$ ($R_{U}$)")
        #ax.set_ylabel(r"$y$ ($R_{U}$)")
        #ax.set_zlabel(r"$z$ ($R_{U}$)")
        #ax.set_xlim(-10,10)
        #ax.set_ylim(-10,10)
        #ax.set_zlim(-10,10)
        #plt.tight_layout()
        #ax.view_init(elev=10,azim=60)
        #plt.title("3D Neutral Tori")
        #ax.tick_params(axis='both', which='major', labelsize=8)
        #plt.savefig("3D Neutrals Small Labels",dpi=1200)

        ##solstice: pole pointed at Sun
        ##equinox: equator pointed at Sun

        fig = plt.figure()
        ax = plt.gca()
        xmid = int(np.shape(self.Y_grid)[0]/2)
        yzplane = ax.contourf(self.Y_grid[xmid,:,:],self.Z_grid[xmid,:,:],moon_density[xmid,:,:],cmap='plasma')
        fig.colorbar(yzplane,label=r"Density (cm$^{-3}$)")
        plt.title(r'Equinox: Neutral Tori in $y$-$z$ plane')
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("neutrals_solstice_yz",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        zmid = int(np.shape(self.Z_grid)[0]/2)
        xyplane = ax.contourf(self.X_grid[:,:,zmid],self.Y_grid[:,:,zmid],moon_density[:,:,zmid],cmap='plasma')
        fig.colorbar(xyplane,label=r"Density (cm$^{-3}$)")
        plt.title(r'Equinox: Neutral Tori in $x$-$y$ plane')
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("neutrals_solstice_xy",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        ymid = int(np.shape(self.Y_grid)[0]/2)
        xzplane = ax.contourf(self.X_grid[:,ymid,:],self.Z_grid[:,ymid,:],moon_density[:,ymid,:],cmap='plasma')
        fig.colorbar(xzplane,label=r"Density (cm$^{-3}$)")
        plt.title(r'Equinox: Neutral Tori in $x$-$z$ plane')
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("neutrals_solstice_xz",dpi=1200)

    def plot_exo(self,n_exo):
        """
        Plots the exosphere density

        Parameters:
            n_exo (ndarray): exosphere density at each point
        """
        fig = plt.figure()
        ax = plt.gca()
        x_mid = int(np.shape(self.X_grid)[0]/2)
        levs = np.linspace(n_exo.min(),n_exo.max(),1000)
        yz_plane = ax.contourf(self.Y_grid[x_mid,:,:],self.Z_grid[x_mid,:,:],n_exo[x_mid,:,:],cmap='plasma',levels=levs)
        fig.colorbar(yz_plane,label=r'Density (cm$^{-3}$)')
        ax.set_xlim(self.y_min,self.y_max)
        ax.set_ylim(self.z_min,self.z_max)
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        plt.title(r"Exosphere $y$-$z$ Projection")
        #plt.savefig("Exosphere y-z",dpi=1200)

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


    def plot_sheath(self,density_grid,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,r0_mag,k_mag,r0_bow,k_bow):
        """
        Plot the magnetosheath region

        Parameters:
            density_grid (ndarray): needs renaming*. Grid that defines whether points are inside or outside the magnetosheath. 
                                    Density set to 1 in the sheath and 0 elsewhere.
            x_mp,y_mp,z_mp,x_bs,y_bs,z_bs (ndarray): MP and BS coordinates used to plot the 3D surface in 2D
            r0_mag,k_mag,r0_bow,k_mag (float): standoff distance and flaring parameter for the magnetopause and bow shock to retrieve the 2D surfaces for plotting.
        """

        fig = plt.figure()
        ax = plt.gca()
        x_mid = int(np.shape(self.X_grid)[0]/2)
        #x_mid = -20
        yz_plane = ax.contourf(self.Y_grid[x_mid,:,:],self.Z_grid[x_mid,:,:],density_grid[x_mid,:,:],cmap='Blues')
        fig.colorbar(yz_plane,label=r'Density (cm$^{-3}$)')
        ax.set_xlim(self.y_min,self.y_max)
        ax.set_ylim(self.z_min,self.z_max)
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        plt.title(r"Magnetosheath $y$-$z$ Projection")
        #plt.savefig("Magnetosheath y-z Projection Low Res",dpi=1200)
        #plt.savefig("High res sheath y-z",dpi=1200)
        #plt.savefig("Magnetosheath y-z Value",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        xy_plane = ax.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],density_grid[:,:,z_pos],cmap='Greens')
        ax.plot(x_mp,y_mp,color='red',alpha=0.1)
        ax.plot(x_bs,y_bs,color='blue',alpha=0.1)
        fig.colorbar(xy_plane,label=r'Density (cm$^{-3}$)')
        ax.set_xlim(self.x_min,self.x_max)
        ax.set_ylim(self.y_min,self.y_max)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        plt.title(r"Magnetosheath $x$-$y$ Projection")       
        #plt.savefig("Magnetosheath Surfaces Plotted",dpi=1200)
        #plt.savefig("Magnetosheath x-y Value",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        from Surface import Surface
        bs = Surface(r0_bow,k_bow)
        mp = Surface(r0_mag,k_mag)
        x_bow, y_bow = bs.surf_2D(self.x_max,self.x_min)
        x_mag, y_mag = mp.surf_2D(self.x_max,self.x_min)        
        ax.plot(x_bow,y_bow,color='red')
        ax.plot(x_mag,y_mag,color='green')
        ax.set_xlim(self.x_min,30)
        ax.set_ylim(self.y_min,self.y_max)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        xy2 = ax.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],density_grid[:,:,z_pos],cmap='Blues')
        fig.colorbar(xy2, label=r'Density (cm$^{-3}$)')
        #plt.title(f'Slice through region at z = {z_slice}')
        plt.title(r"Magnetosheath $x$-$y$ Projection")
        #plt.savefig("Magnetosheath x-y Reduced Grid",dpi=1200)

    def plot_ver(self,ver,r0_mag,k_mag,r0_bow,k_bow):
        """
        Plots the volumetric emission

        Parameters:
            ver (ndarray): volumetric emission in photon cm^-3 s^-1
            r0_mag,k_mag (float): magnetopause standoff distance and flaring parameter
            r0_bow,k_bow (float): bow shock standoff distance and flaring parameter
        """
        #fig = plt.figure()
        #ax = fig.add_subplot(111,projection='3d')
        #ax.set_xlim([self.x_min, self.x_max])
        #ax.set_ylim([self.y_min, self.y_max])
        #ax.set_zlim([self.z_min, self.z_max])
        #levs = np.linspace(ver.min(),ver.max(),100)
        #ver_cont = plt.contourf(self.X_grid, self.Y_grid, self.Z_grid, ver, cmap='YlOrRd',levels = levs)
        #ver_scatter = plt.scatter(self.X_grid,self.Y_grid,self.Z_grid,c=ver,cmap='YlOrRd')
        #fig.colorbar(ver_scatter, label=r'Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)')
        #ax.set_xlabel(r'$x$ ($R_{U}$)')
        #ax.set_ylabel(r'$y$ ($R_{U}$)')
        #ax.set_zlabel(r'$z$ ($R_{U}$)')
        #plt.savefig("3D VER not working",dpi=1200)

        fig = plt.figure(figsize=(6.4,5.8)) #default: width = 6.4inches, height = 4.8inches
        ax = plt.gca()
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        xy_plane = ax.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver[:,:,z_pos],cmap='YlOrRd')
        from Surface import Surface
        bs = Surface(r0_bow,k_bow)
        mp = Surface(r0_mag,k_mag)
        x_bow, y_bow = bs.surf_2D(self.x_max,self.x_min)
        x_mag, y_mag = mp.surf_2D(self.x_max,self.x_min)        
        ax.plot(x_bow,y_bow,color='red')
        ax.plot(x_mag,y_mag,color='blue')
        fig.colorbar(xy_plane,label=r'Emission Rate (photon cm$^{-3}$ s$^{-1}$)')
        ax.set_xlim(self.x_min,40)
        ax.set_ylim(self.y_min,self.y_max)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        #plt.title(r"VER $x$-$y$ Projection: $v_{\mathrm{SW}}=400$ km s$^{-1}$")
        plt.title(r"VER: $v_{\mathrm{SW}}$ = 690 km s$^{-1}$, $n_{\mathrm{SW,1AU}}=1.856$ cm$^{-3}$",fontsize=12)
        #plt.title("Volumetric Emission: Reduced Grid Limits")
        ver_max = ver.max()
        ver_mean = ver.mean()
        ver_max_3 = f"{ver_max:.3g}"
        ver_mean_3 = f"{ver_mean:.3g}"
        plt.gcf().text(0.05, 0.05, f"Max VER: {ver_max_3}", ha='left', fontsize=12)
        plt.gcf().text(0.05, 0.01, f"Mean VER: {ver_mean_3}", ha='left', fontsize=12)   
        #plt.savefig("VER_x-y_high_vsw_v3",dpi=1200)






        