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
from matplotlib.gridspec import GridSpec
from PIL import Image
import io

class Plotter:
    def __init__(self,grid_limits=[-10,10,-10,10,-10,10],X_grid=None,Y_grid=None,Z_grid=None,config=None):
        self.x_min = grid_limits[0]
        self.x_max = grid_limits[1]
        self.y_min = grid_limits[2]
        self.y_max = grid_limits[3]
        self.z_min = grid_limits[4]
        self.z_max = grid_limits[5]
        #self.xbox, self.ybox, self.zbox = xbox,ybox,zbox
        self.X_grid, self.Y_grid, self.Z_grid = X_grid, Y_grid, Z_grid
        self.config = config
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

        #fig = plt.figure()
        #ax = plt.gca()
        #xmid = int(np.shape(self.Y_grid)[0]/2)
        #yzplane = ax.contourf(self.Y_grid[xmid,:,:],self.Z_grid[xmid,:,:],moon_density[xmid,:,:],cmap='plasma')
        #fig.colorbar(yzplane,label=r"Density (cm$^{-3}$)")
        #plt.title(r'Equinox: Neutral Tori in $y$-$z$ plane')
        #ax.set_xlabel(r'$y$ ($R_{U}$)')
        #ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("neutrals_equinox_yz_reduced",dpi=1200)

        #fig = plt.figure()
        #ax = plt.gca()
        #zmid = int(np.shape(self.Z_grid)[0]/2)
        #xyplane = ax.contourf(self.X_grid[:,:,zmid],self.Y_grid[:,:,zmid],moon_density[:,:,zmid],cmap='plasma')
        #fig.colorbar(xyplane,label=r"Density (cm$^{-3}$)")
        #plt.title(r'Equinox: Neutral Tori in $x$-$y$ plane')
        #ax.set_xlabel(r'$x$ ($R_{U}$)')
        #ax.set_ylabel(r'$y$ ($R_{U}$)')
        #plt.savefig("neutrals_equinox_xy_reduced",dpi=1200)

        #fig = plt.figure()
        #ax = plt.gca()
        #ymid = int(np.shape(self.Y_grid)[0]/2)
        #xzplane = ax.contourf(self.X_grid[:,ymid,:],self.Z_grid[:,ymid,:],moon_density[:,ymid,:],cmap='plasma')
        #fig.colorbar(xzplane,label=r"Density (cm$^{-3}$)")
        #plt.title(r'Equinox: Neutral Tori in $x$-$z$ plane')
        #ax.set_xlabel(r'$x$ ($R_{U}$)')
        #ax.set_ylabel(r'$z$ ($R_{U}$)')
        #plt.savefig("neutrals_equinox_xz_reduced",dpi=1200)

        #fig,ax = plt.subplots(1,3,figsize=(18,10),gridspec_kw={'width_ratios': [1, 1, 1]})
        fig = plt.figure(figsize=(14,5))
        gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 0.1], wspace=0.5)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        x_pos = int(np.shape(self.Y_grid)[0]/2)
        y_pos = int(np.shape(self.Y_grid)[0]/2)
        yzplane = ax0.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],moon_density[x_pos,:,:],cmap='plasma')
        xyplane = ax1.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],moon_density[:,:,z_pos],cmap='plasma')
        xzplane = ax2.contourf(self.X_grid[:,y_pos,:],self.Z_grid[:,y_pos,:],moon_density[:,y_pos,:],cmap='plasma')
        #cbar = fig.colorbar(xzplane, ax=[ax0, ax1, ax2], location='right', shrink=0.65, pad=0.02)
        #cbar.set_label(r'Neutral Density (cm$^{-3}$)')
        cax1 = fig.add_axes([0.85, 0.22, 0.01, 0.55])  # Manually define position of colorbar
        cbar1 = fig.colorbar(yzplane, cax=cax1,label=r"Density (cm$^{-3}$)",shrink=0.3)
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax0.set_xlabel(r'$y$ ($R_{U}$)')
        ax0.set_ylabel(r'$z$ ($R_{U}$)')
        ax1.set_xlabel(r'$x$ ($R_{U}$)')
        ax1.set_ylabel(r'$y$ ($R_{U}$)')
        ax2.set_xlabel(r'$x$ ($R_{U}$)')
        ax2.set_ylabel(r'$z$ ($R_{U}$)')
        ax0.set_title(r"(a) $y$-$z$ plane")
        ax1.set_title(r"(b) $x$-$y$ plane")
        ax2.set_title(r"(c) $x$-$z$ plane")
        if self.config == "E":
            fig.suptitle("Moon-Sourced Neutral Tori: Equinox",fontsize=13,y=0.9)
        elif self.config == "S":
            fig.suptitle("Moon-Sourced Neutral Tori: Solstice",fontsize=13,y=0.9)
        #plt.tight_layout()
        #plt.savefig("neutrals_equinox_reduced",dpi=1200)

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

        fig = plt.figure(figsize=(17,5))
        gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 0.1], wspace=0.5)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        x_pos = int(np.shape(self.Y_grid)[0]/2)
        y_pos = int(np.shape(self.Y_grid)[0]/2)
        yzplane = ax0.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],ver[x_pos,:,:],cmap='YlOrRd')
        xyplane = ax1.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver[:,:,z_pos],cmap='YlOrRd')
        xzplane = ax2.contourf(self.X_grid[:,y_pos,:],self.Z_grid[:,y_pos,:],ver[:,y_pos,:],cmap='YlOrRd')
        if self.config == "N":
            ax2.plot(x_bow,y_bow,color='red')
            ax2.plot(x_mag,y_mag,color='blue')
            ax1.plot(x_bow,y_bow,color='red')
            ax1.plot(x_mag,y_mag,color='blue')
        if self.config == "E":
            ax1.plot(x_bow,y_bow,color='red')
            ax1.plot(x_mag,y_mag,color='blue')
            plt.suptitle(r"Volumetric Emission Rate: Equinox, $v_{\mathrm{SW}} = 400$ km s$^{-1}$")
        if self.config == "S":
            ax2.plot(x_bow,y_bow,color='red')
            ax2.plot(x_mag,y_mag,color='blue')
            plt.suptitle(r"Volumetric Emission Rate: Solstice, $v_{\mathrm{SW}} = 400$ km s$^{-1}$")
        cbar = fig.colorbar(xzplane, ax=[ax0, ax1, ax2], location='right', shrink=0.65, pad=0.02)
        cbar.set_label(r'Emission Rate (photon cm$^{-3}$ s$^{-1}$)')
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax1.set_xlim(self.x_min,self.x_max)
        ax1.set_ylim(self.y_min,self.y_max)
        ax2.set_xlim(self.x_min,self.x_max)
        ax2.set_ylim(self.z_min,self.z_max)
        ax0.set_xlabel(r'$y$ ($R_{U}$)')
        ax0.set_ylabel(r'$z$ ($R_{U}$)')
        ax1.set_xlabel(r'$x$ ($R_{U}$)')
        ax1.set_ylabel(r'$y$ ($R_{U}$)')
        ax2.set_xlabel(r'$x$ ($R_{U}$)')
        ax2.set_ylabel(r'$z$ ($R_{U}$)')
        ax0.set_title(r"(a) $y$-$z$ plane")
        ax1.set_title(r"(b) $x$-$y$ plane")
        ax2.set_title(r"(c) $x$-$z$ plane")
        if self.config == "S":
            plt.suptitle(r"Volumetric Emission Rate: Solstice, $v_{\mathrm{SW}} = 400$ km s$^{-1}$")
        if self.config == "E":
            plt.suptitle(r"Volumetric Emission Rate: Equinox, $v_{\mathrm{SW}} = 400$ km s$^{-1}$")
        ver_max_s = ver.max()
        ver_mean_s = ver.mean()
        ver_max_3_s = f"{ver_max_s:.3g}"
        ver_mean_3_s = f"{ver_mean_s:.3g}"
        plt.gcf().text(0.01, 0.05, f"Max VER: {ver_max_3_s}", ha='left', fontsize=12)
        plt.gcf().text(0.01, 0.01, f"Mean VER: {ver_mean_3_s}", ha='left', fontsize=12)   
        #plt.savefig("VER_equinox_allplanes",dpi=1200)
        #plt.savefig("VER_allplanes_no_dZ",dpi=1200)

    def plot_ver_combined(self,ver_slow,ver_fast,r_mag_f,k_mag_f,r_bow_f,k_bow_f):

        fig,ax = plt.subplots(1,2,figsize=(10,5))
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        from Model import magnetopause, bow_shock
        from Surface import Surface
        mp = Surface(magnetopause.r0,magnetopause.K)
        bs = Surface(bow_shock.r0,bow_shock.K)
        xy_plane_slow = ax[0].contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver_slow[:,:,z_pos],cmap='YlOrRd')
        x_bow_s, y_bow_s = bs.surf_2D(self.x_max,self.x_min)
        x_mag_s, y_mag_s = mp.surf_2D(self.x_max,self.x_min)        
        ax[0].plot(x_bow_s,y_bow_s,color='red')
        ax[0].plot(x_mag_s,y_mag_s,color='blue')
        plt.colorbar(xy_plane_slow,label=r'Emission Rate (photon cm$^{-3}$ s$^{-1}$)',ax=ax[0])
        ax[0].set_xlim(self.x_min,40)
        ax[0].set_ylim(self.y_min,self.y_max)
        ax[0].set_xlabel(r'$x$ ($R_{U}$)')
        ax[0].set_ylabel(r'$y$ ($R_{U}$)')
        ax[0].set_title(r"(a) Slow Wind: $v_{\mathrm{SW}}=400$ km s$^{-1}$")
        xy_plane_fast = ax[1].contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver_fast[:,:,z_pos],cmap='YlOrRd')
        x_bow_f, y_bow_f = bs.surf_2D(self.x_max,self.x_min,r_bow_f,k_bow_f)
        x_mag_f, y_mag_f = mp.surf_2D(self.x_max,self.x_min,r_mag_f,k_mag_f)        
        ax[1].plot(x_bow_f,y_bow_f,color='red')
        ax[1].plot(x_mag_f,y_mag_f,color='blue')
        plt.colorbar(xy_plane_fast,label=r'Emission Rate (photon cm$^{-3}$ s$^{-1}$)',ax=ax[1])
        ax[1].set_xlim(self.x_min,40)
        ax[1].set_ylim(self.y_min,self.y_max)
        ax[1].set_xlabel(r'$x$ ($R_{U}$)')
        ax[1].set_ylabel(r'$y$ ($R_{U}$)')
        ax[1].set_title(r"(b) Fast Wind: $v_{\mathrm{SW}}=800$ km s$^{-1}$")
        ver_max_s = ver_slow.max()
        ver_mean_s = ver_slow.mean()
        ver_max_3_s = f"{ver_max_s:.3g}"
        ver_mean_3_s = f"{ver_mean_s:.3g}"
        ver_max_f = ver_fast.max()
        ver_mean_f = ver_fast.mean()
        ver_max_3_f = f"{ver_max_f:.3g}"
        ver_mean_3_f = f"{ver_mean_f:.3g}"
        plt.gcf().text(0.01, 0.05, f"Max Slow VER: {ver_max_3_s}", ha='left', fontsize=12)
        plt.gcf().text(0.01, 0.01, f"Mean Slow VER: {ver_mean_3_s}", ha='left', fontsize=12)   
        plt.gcf().text(0.99, 0.05, f"Max Fast VER: {ver_max_3_f}", ha='right', fontsize=12)
        plt.gcf().text(0.99, 0.01, f"Mean Fast VER: {ver_mean_3_f}", ha='right', fontsize=12)
        plt.tight_layout()
        #plt.savefig("VER_combined",dpi=1200)

    def plot_combined_ver_flux(self,ver_slow,ver_fast,flux_slow,flux_fast,int_s_s,int_h_s,int_s_f,int_h_f):
        fig = plt.figure(figsize=(15,10))
        gs = GridSpec(2, 4, width_ratios=[1, 1, 1, 1], wspace=0.6)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])
        ax3 = fig.add_subplot(gs[0, 3])
        ax4 = fig.add_subplot(gs[1, 0])
        ax5 = fig.add_subplot(gs[1, 1])
        ax6 = fig.add_subplot(gs[1, 2])
        ax7 = fig.add_subplot(gs[1, 3])
        z_pos = int(np.shape(self.Z_grid)[0]/2)
        x_pos = int(np.shape(self.Y_grid)[0]/2)
        y_pos = int(np.shape(self.Y_grid)[0]/2)
        yzplane_v_s = ax0.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],ver_slow[x_pos,:,:],cmap='YlOrRd',levels=20)
        xyplane_v_s = ax1.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver_slow[:,:,z_pos],cmap='YlOrRd',levels=20)
        xzplane_v_s = ax2.contourf(self.X_grid[:,y_pos,:],self.Z_grid[:,y_pos,:],ver_slow[:,y_pos,:],cmap='YlOrRd',levels=20)
        yzplane_f_s = ax3.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],flux_slow,cmap='plasma',levels=20)  
        yzplane_v_f = ax4.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],ver_fast[x_pos,:,:],cmap='YlOrRd',levels=20)
        xyplane_v_f = ax5.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver_fast[:,:,z_pos],cmap='YlOrRd',levels=20)
        xzplane_v_f = ax6.contourf(self.X_grid[:,y_pos,:],self.Z_grid[:,y_pos,:],ver_fast[:,y_pos,:],cmap='YlOrRd',levels=20)
        yzplane_f_f = ax7.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],flux_fast,cmap='plasma',levels=20)  
        if self.config == "E":
            plt.suptitle(r"VER and Flux: Equinox, $v_{\mathrm{SW, Slow}}=400$ km s$^{-1}$, $v_{\mathrm{SW, Fast}}=800$ km s$^{-1}$, $x$,$y$,$z$ Slice Position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            #plt.suptitle(r"VER & Flux: Equinox, $v_{\mathrm{SW}} = 690$ km s$^{-1}$, $n_{\mathrm{SW,1 AU}}=1.86$ cm$^{-3}$",x=0.5,y=0.9)
        if self.config == "S":
            plt.suptitle(r"VER and Flux: Solstice, $v_{\mathrm{SW, Slow}}=400$ km s$^{-1}$, $v_{\mathrm{SW, Fast}}=800$ km s$^{-1}$, $x$,$y$,$z$ Slice Position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            #plt.suptitle(r"VER & Flux: Solstice, $v_{\mathrm{SW}} = 690$ km s$^{-1}$, $n_{\mathrm{SW,1 AU}}=1.86$ cm$^{-3}$",x=0.5,y=0.9)
        cax1 = fig.add_axes([0.93, 0.57, 0.01, 0.275])  # Manually define position of colorbar
        cbar1 = fig.colorbar(yzplane_f_s, cax=cax1,label=r"Flux (photon cm$^{-2}$ s$^{-1}$)",shrink=0.3)
        cax4 = fig.add_axes([0.93, 0.14, 0.01, 0.275])  # Manually define position of colorbar
        cbar4 = fig.colorbar(yzplane_f_f, cax=cax4,label=r"Flux (photon cm$^{-2}$ s$^{-1}$)",shrink=0.3)
        cax2 = fig.add_axes([0.04, 0.57, 0.01, 0.275])  # Manually define position of colorbar #left, bottom, width, height
        cax3 = fig.add_axes([0.04, 0.14, 0.01, 0.275])
        if self.config == "E":
            cbar2 = fig.colorbar(xzplane_v_s, cax=cax2,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
            cbar3 = fig.colorbar(xzplane_v_f, cax=cax3,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
        if self.config == "S":
            cbar2 = fig.colorbar(yzplane_v_s, cax=cax2,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
            cbar3 = fig.colorbar(yzplane_v_f, cax=cax3,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
        cbar2.ax.yaxis.set_label_position('left')  # Move the label to the left
        cbar3.ax.yaxis.set_label_position('left')
        #cbar2.ax.set_ylim(cbar2.get_ticks()[0], cbar2.get_ticks()[0] + 0.5)  # Shrink the range if necessary
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax3.set_aspect('equal')
        ax0.set_xlim(self.x_min,self.x_max)
        ax0.set_ylim(self.y_min,self.y_max)
        ax1.set_xlim(self.x_min,self.x_max)
        ax1.set_ylim(self.z_min,self.z_max)
        ax2.set_xlim(self.y_min,self.y_max)
        ax2.set_ylim(self.z_min,self.z_max)
        ax0.set_xlabel(r'$y$ ($R_{U}$)')
        ax0.set_ylabel(r'$z$ ($R_{U}$)')
        ax1.set_xlabel(r'$x$ ($R_{U}$)')
        ax1.set_ylabel(r'$y$ ($R_{U}$)')
        ax2.set_xlabel(r'$x$ ($R_{U}$)')
        ax2.set_ylabel(r'$z$ ($R_{U}$)')
        ax3.set_xlabel(r'$y$ ($R_{U}$)')
        ax3.set_ylabel(r'$z$ ($R_{U}$)')
        ax4.set_aspect('equal')
        ax5.set_aspect('equal')
        ax6.set_aspect('equal')
        ax7.set_aspect('equal')
        ax4.set_xlim(self.x_min,self.x_max)
        ax4.set_ylim(self.y_min,self.y_max)
        ax5.set_xlim(self.x_min,self.x_max)
        ax5.set_ylim(self.z_min,self.z_max)
        ax6.set_xlim(self.y_min,self.y_max)
        ax6.set_ylim(self.z_min,self.z_max)
        ax3.set_xlabel(r'$y$ ($R_{U}$)')
        ax3.set_ylabel(r'$z$ ($R_{U}$)')
        ax4.set_xlabel(r'$x$ ($R_{U}$)')
        ax4.set_ylabel(r'$y$ ($R_{U}$)')
        ax5.set_xlabel(r'$x$ ($R_{U}$)')
        ax5.set_ylabel(r'$z$ ($R_{U}$)')
        ax6.set_xlabel(r'$y$ ($R_{U}$)')
        ax6.set_ylabel(r'$z$ ($R_{U}$)')
        ax0.set_title(r"(a) $y$-$z$ plane (Slow)")
        ax1.set_title(r"(b) $x$-$y$ plane (Slow)")
        ax2.set_title(r"(c) $x$-$z$ plane (Slow)")
        ax3.set_title(r"(d) Flux (Slow)")
        ax4.set_title(r"(e) $y$-$z$ plane (Fast)")
        ax5.set_title(r"(f) $x$-$y$ plane (Fast)")
        ax6.set_title(r"(g) $x$-$z$ plane (Fast)")
        ax7.set_title(r"(h) Flux (Fast)")
        plt.gcf().text(0.05, 0.53, f"Max VER (Slow): {f"{ver_slow.max():.3g}"}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)
        plt.gcf().text(0.05, 0.50, f"Mean VER (Slow): {ver_slow.mean():.3g}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)   
        plt.gcf().text(0.05, 0.08, f"Max VER (Fast): {f"{ver_fast.max():.3g}"}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)
        plt.gcf().text(0.05, 0.05, f"Mean VER (Fast): {ver_fast.mean():.3g}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)   
        plt.gcf().text(0.75, 0.53, f"Integration time (Slow) (s): {int_s_s} s", ha='left', fontsize=12)
        plt.gcf().text(0.75, 0.50, f"Integration time (Slow) (h): {int_h_s} h", ha='left', fontsize=12) 
        plt.gcf().text(0.75, 0.08, f"Integration time (Fast) (s): {int_s_f} s", ha='left', fontsize=12)
        plt.gcf().text(0.75, 0.05, f"Integration time (Fast) (h): {int_h_f} h", ha='left', fontsize=12) 
        #return fig
        #if self.config == "S":
            #plt.savefig("ver_flux_combined_solstice",dpi=1200)  
        #elif self.config == "E":
            #plt.savefig("ver_flux_combined_equinox",dpi=1200)



    def plot_flux(self,flux):
        """
        Plots the flux projection obtained from summing VER along the x-axis

        Parameters
            flux (ndarray): flux detected by SXI
        """
        fig = plt.figure()
        ax = plt.gca()
        levs = np.linspace(flux.min(),flux.max(),20)
        xmid = int(np.shape(self.Y_grid)[0]/2)
        yzplane = ax.contourf(self.Y_grid[xmid,:,:],self.Z_grid[xmid,:,:],flux,cmap='plasma',levels=levs)
        fig.colorbar(yzplane,label=r"Flux (photon cm$^{-2}$ s$^{-1}$)")
        ax.set_xlim(self.y_min,self.y_max)
        ax.set_ylim(self.z_min,self.z_max)
        ax.set_xlabel(r"$y$ ($R_{U}$)")
        ax.set_ylabel(r"$z$ ($R_{U}$)")
        #plt.title(r"Flux: Equinox, SMILE-like SXI 300 $R_{U}$ Upstream",fontsize=11)
        if self.config == "S":
            plt.title(r"Flux: Solstice, SMILE-like SXI 300 $R_{U}$ Upstream",fontsize=11)
        if self.config == "E":
            plt.title(r"Flux: Equinox, SMILE-like SXI 300 $R_{U}$ Upstream",fontsize=11)
        #plt.savefig("flux_image_no_dZ",dpi=1200)
        #plt.savefig("flux_solstice",dpi=1200)

    def plot_flux_ver(self,ver,flux,r0_mag,k_mag,r0_bow,k_bow):
        fig, ax = plt.subplots(1,2,figsize=(10,5))
        #plt.rcParams.update({'font.size': 13})
        contour_levels1 = np.linspace(ver.min(), ver.max(), 10)
        if self.config == "S":
            plt.suptitle(r'VER & Flux: Solstice Configuration',fontsize=13)
            pos = int(np.shape(self.Y_grid)[0]/2)
            xz_plane = ax[0].contourf(self.X_grid[:,pos,:],self.Z_grid[:,pos,:],ver[:,pos,:],cmap='YlOrRd',levels=contour_levels1)
            plt.colorbar(xz_plane, label='Volumetric Emission (photon cm'r"$^{-3}$"' s'r"$^{-1}$"')', ax=ax[0])
            print("S")
        elif self.config == "E":
            plt.suptitle(r'VER & Flux: Equinox Configuration',fontsize=13)
            pos = int(np.shape(self.Z_grid)[0]/2)
            xy_plane = ax[0].contourf(self.X_grid[:,:,pos],self.Y_grid[:,:,pos],ver[:,:,pos],cmap='YlOrRd',levels=contour_levels1)
            plt.colorbar(xy_plane, label='Volumetric Emission (photon cm'r"$^{-3}$"' s'r"$^{-1}$"')', ax=ax[0])
            print("E")
        ax[0].set_xlim([self.x_min, self.x_max])
        ax[0].set_ylim([self.y_min, self.y_max])
        #from Model import magnetopause, bow_shock ## can't use this as it reruns the whole program :(
        from Surface import Surface
        bs = Surface(r0_bow,k_bow)
        mp = Surface(r0_mag,k_mag)
        x_bow, y_bow = bs.surf_2D(self.x_max,self.x_min)
        x_mag, y_mag = mp.surf_2D(self.x_max,self.x_min)
        x_bow, y_bow = bs.surf_2D(self.x_max,self.x_min)
        x_mag, y_mag = mp.surf_2D(self.x_max,self.x_min)   
        ax[0].plot(x_mag, y_mag, color='blue')
        ax[0].plot(x_bow, y_bow, color='red')
        ax[0].set_xlabel(r"$x$ ($R_{U}$)",fontsize=13)
        ax[0].set_ylabel(r"$y$ ($R_{U}$)",fontsize=13)
        ax[0].set_title('(a) Volumetric Emission',fontsize=13)
        ax[0].grid(False)
        ax[0].tick_params(axis='both', direction='out', top=True, bottom=True, left=True, right=True,labelsize=13)
        x_pos = int(np.shape(self.Y_grid)[0]/2)
        levs = np.linspace(flux.min(),flux.max(),10)
        yz_plane = ax[1].contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],flux,cmap='plasma',levels=levs)
        plt.colorbar(yz_plane, label='Detected Flux (photon cm'r"$^{-2}$"' s'r"$^{-1}$"')', ax=ax[1])
        ax[1].set_xlim([self.y_min, self.y_max])
        ax[1].set_ylim([self.z_min, self.z_max])
        ax[1].set_title('(b) SMILE-like SXI Image',fontsize=13)
        #axes[1].set_title('(b)')
        ax[1].set_ylabel(r'$z$ ($R_{U}$)',fontsize=13)
        ax[1].set_xlabel(r"$y$ ($R_{U}$)",fontsize=13)
        ax[1].grid(False)
        ax[1].tick_params(axis='both', direction='out', top=True, bottom=True, left=True, right=True,labelsize=13)
        #plt.gcf().text(0.48, 0.05, f"Peak Integration time (s): {integration_time_sec_3sf}", ha='center', fontsize=12)
        #plt.gcf().text(0.48, 0.01, f"Peak Integration time (h): {integration_time_hour_3sf}", ha='center', fontsize=12)
        #plt.subplots_adjust(bottom=0.2)
        plt.tight_layout()
        #plt.savefig('ver_flux_equinox',dpi=1200)


    def plot_all_ver_flux(self,ver,flux,int_s,int_h,x_pos,y_pos,z_pos,v_sw):

        fig = plt.figure(figsize=(15,5))
        gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 1], wspace=0.6)
        #from Surface import Surface
        #bs = Surface(r0_bow,k_bow)
        #mp = Surface(r0_mag,k_mag)
        #x_bow, y_bow = bs.surf_2D(self.x_max,self.x_min)
        #x_mag, y_mag = mp.surf_2D(self.x_max,self.x_min)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])
        ax3 = fig.add_subplot(gs[0, 3])
        #z_pos = int(np.shape(self.Z_grid)[0]/2)
        #x_pos = int(np.shape(self.Y_grid)[0]/2)
        #y_pos = int(np.shape(self.Y_grid)[0]/2)
        #print(x_pos,y_pos,z_pos)
        #x_pos, y_pos, z_pos = 40,40,40
        levs1 = np.linspace(flux.min(),flux.max(),20)
        levs2 = np.linspace(ver.min(),ver.max(),20)
        yzplane_v = ax0.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],ver[x_pos,:,:],cmap='YlOrRd',levels=levs2)
        xyplane_v = ax1.contourf(self.X_grid[:,:,z_pos],self.Y_grid[:,:,z_pos],ver[:,:,z_pos],cmap='YlOrRd',levels=levs2)
        xzplane_v = ax2.contourf(self.X_grid[:,y_pos,:],self.Z_grid[:,y_pos,:],ver[:,y_pos,:],cmap='YlOrRd',levels=levs2)
        yzplane_f = ax3.contourf(self.Y_grid[x_pos,:,:],self.Z_grid[x_pos,:,:],flux,cmap='plasma',levels=levs1)  
        if self.config == "E":
            if v_sw != None:
                plt.suptitle(r"VER and Flux: Equinox, $v_{\mathrm{SW}}=$"+f"{v_sw/1000:.3g}"+r" km s$^{-1}$, $x$,$y$,$z$ position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            else:
                plt.suptitle(r"VER and Flux: Equinox, $v_{\mathrm{SW}}=$400 km s$^{-1}$, $x$,$y$,$z$ position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            #plt.suptitle(r"VER & Flux: Equinox, $v_{\mathrm{SW}} = 690$ km s$^{-1}$, $n_{\mathrm{SW,1 AU}}=1.86$ cm$^{-3}$",x=0.5,y=0.9)
        if self.config == "S":
            if v_sw != None:
                plt.suptitle(r"VER and Flux: Solstice, $v_{\mathrm{SW}}=$"+f"{v_sw/1000:.3g}"+r" km s$^{-1}$, $x$,$y$,$z$ position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            else:
                plt.suptitle(r"VER and Flux: Solstice, $v_{\mathrm{SW}}=$400 km s$^{-1}$, $x$,$y$,$z$ position = "+f"{x_pos},{y_pos},{z_pos}",x=0.5,y=0.9)
            #plt.suptitle(r"VER & Flux: Solstice, $v_{\mathrm{SW}} = 690$ km s$^{-1}$, $n_{\mathrm{SW,1 AU}}=1.86$ cm$^{-3}$",x=0.5,y=0.9)
        cax1 = fig.add_axes([0.91, 0.22, 0.01, 0.55])  # Manually define position of colorbar
        cbar1 = fig.colorbar(yzplane_f, cax=cax1,label=r"Flux (photon cm$^{-2}$ s$^{-1}$)",shrink=0.3)
        cax2 = fig.add_axes([0.04, 0.22, 0.01, 0.55])  # Manually define position of colorbar
        if self.config == "E":
            cbar2 = fig.colorbar(xzplane_v, cax=cax2,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
        if self.config == "S":
            cbar2 = fig.colorbar(yzplane_v, cax=cax2,label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)",shrink=0.3)
        cbar2.ax.yaxis.set_label_position('left')  # Move the label to the left
        #cbar2.ax.set_ylim(cbar2.get_ticks()[0], cbar2.get_ticks()[0] + 0.5)  # Shrink the range if necessary
        ax0.set_aspect('equal')
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')
        ax3.set_aspect('equal')
        ax0.set_xlim(self.x_min,self.x_max)
        ax0.set_ylim(self.y_min,self.y_max)
        ax1.set_xlim(self.x_min,self.x_max)
        ax1.set_ylim(self.z_min,self.z_max)
        ax2.set_xlim(self.y_min,self.y_max)
        ax2.set_ylim(self.z_min,self.z_max)
        ax0.set_xlabel(r'$y$ ($R_{U}$)')
        ax0.set_ylabel(r'$z$ ($R_{U}$)')
        ax1.set_xlabel(r'$x$ ($R_{U}$)')
        ax1.set_ylabel(r'$y$ ($R_{U}$)')
        ax2.set_xlabel(r'$x$ ($R_{U}$)')
        ax2.set_ylabel(r'$z$ ($R_{U}$)')
        ax3.set_xlabel(r'$y$ ($R_{U}$)')
        ax3.set_ylabel(r'$z$ ($R_{U}$)')
        ax0.set_title(r"(a) $y$-$z$ plane")
        ax1.set_title(r"(b) $x$-$y$ plane")
        ax2.set_title(r"(c) $x$-$z$ plane")
        ax3.set_title(r"(d) Flux ($y$-$z$ plane)")
        plt.gcf().text(0.05, 0.10, f"Max VER: {f"{ver.max():.3g}"}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)
        plt.gcf().text(0.05, 0.05, f"Mean VER: {ver.mean():.3g}"+r" photon cm$^{-3}$ s$^{-1}$", ha='left', fontsize=12)   
        plt.gcf().text(0.75, 0.10, f"Integration time (s): {int_s}s", ha='left', fontsize=12)
        plt.gcf().text(0.75, 0.05, f"Integration time (h): {int_h}h", ha='left', fontsize=12) 
        #plt.savefig("ver_flux_equinox_j2",dpi=1200)
        #plt.savefig("ver_flux_high_vsw_solstice",dpi=1200)

    def plot_flux_gif(self,ver,flux,int_s,int_h,config):
        x_pos,y_pos,z_pos = 10,10,10
        flux_figs = []
        for x in range(7):
            fig = self.plot_all_ver_flux(ver,flux,int_s,int_h,x_pos,y_pos,z_pos)
            flux_figs.append(fig)
            x_pos += 10
            y_pos += 10
            z_pos += 10

        # Step 2: Save each figure to an in-memory buffer
        images = []
        for fig in flux_figs:
            buf = io.BytesIO()
            fig.savefig(buf, format='png')  # Save to buffer
            buf.seek(0)  # Move to the start of the buffer
            images.append(Image.open(buf))
            plt.close(fig)  # Close the figure to free memory

        # Step 3: Create the GIF
        if config == "S":
            gif_path = "ver_flux_gif_solstice.gif"
        elif config == "E":
            gif_path = "ver_flux_gif_equinox.gif"
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            duration=2000,  # Duration of each frame in milliseconds
            loop=0          # Infinite loop
        )






