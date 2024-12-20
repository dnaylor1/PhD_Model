import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.spatial import distance
from scipy import constants

class Magnetosheath:
    def __init__(self,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,X_grid,Y_grid,Z_grid,rad):
        self.x_mp = x_mp
        self.y_mp = y_mp
        self.z_mp = z_mp
        self.x_bs = x_bs
        self.y_bs = y_bs
        self.z_bs = z_bs
        self.X_grid = X_grid
        self.Y_grid = Y_grid
        self.Z_grid = Z_grid
        self.rad = rad

    def sheath_surface(self):
        """
        Defines the magnetosheath surface, the area between the magnetopause and the bow shock

        Returns:
            density_grid (ndarray): array of magnetosheath points with density set to 0.1 cm^-3 inside the sheath.
        """
        # Flatten the surface data for magnetopause
        x_mp_flat = self.x_mp.flatten()
        y_mp_flat = self.y_mp.flatten()
        z_mp_flat = self.z_mp.flatten()
        # Flatten the surface data for bow shock
        x_bs_flat = self.x_bs.flatten()
        y_bs_flat = self.y_bs.flatten()
        z_bs_flat = self.z_bs.flatten()
        # Stack the magnetopause and bow shock points together
        points_mp = np.vstack([x_mp_flat, y_mp_flat, z_mp_flat]).T
        points_bs = np.vstack([x_bs_flat, y_bs_flat, z_bs_flat]).T
        # Remove any NaN points before using griddata
        points_mp = points_mp[np.all(np.isfinite(points_mp), axis=1)]
        points_bs = points_bs[np.all(np.isfinite(points_bs), axis=1)]
        # Placeholder values
        values_mp = np.ones(len(points_mp))
        values_bs = np.ones(len(points_bs)) 
        # Interpolate the surfaces
        interpolated_bs = griddata(
            points_bs, values_bs, (self.X_grid, self.Y_grid, self.Z_grid), method='linear'
        )
 
        interpolated_mp = griddata(
            points_mp, values_mp, (self.X_grid, self.Y_grid, self.Z_grid), method='linear'
        )

        ###########################################################################################

        density_grid = np.zeros_like(self.rad)

        #density_grid[np.isfinite(griddata(points_bs, np.ones(len(points_bs)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear')) & 
                    #~np.isfinite(griddata(points_mp, np.ones(len(points_mp)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear'))] = True
        
        sheath_region = (
            np.isfinite(griddata(points_bs, np.ones(len(points_bs)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear')) &
            ~np.isfinite(griddata(points_mp, np.ones(len(points_mp)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear'))
        )

        density_grid[sheath_region] = 0.1

        indices = np.where(density_grid)
        x_points = self.X_grid[indices]
        y_points = self.Y_grid[indices]
        z_points = self.Z_grid[indices]

        ###########################################################################################

        magnetosheath_region = (interpolated_mp <= self.rad) & (self.rad <= interpolated_bs)
        density = np.zeros_like(self.rad)
        density[magnetosheath_region] = 1
        x_sheath = self.X_grid[magnetosheath_region]
        y_sheath = self.Y_grid[magnetosheath_region]
        z_sheath = self.Z_grid[magnetosheath_region]

        return density_grid
    
    def volumetric_emission(self,n_n,n_q,n_p=None,T_sw=None,v_sw=None,v_sf=None,j_s=None):
        """
        Calculates volumetric emission of soft x-rays in the magnetosheath

        Parameters:
            n_n (ndarray): neutral density array
            n_q (ndarray): sheath ion density array
            n_p (float): solar wind density in the case of solar wind variations
            T_sw (float): solar wind temperature in the case of variations
            v_sw (float): solar wind speed in the case of variations

        Returns:
            ver (ndarray): volumetric emission at each point in the magnetosheath
        """
        if j_s == True:
            n_sw = n_p * 1e-6
            T_sheath = 5.45e4
            n_q = n_q * (1/n_q.max()) * n_sw
        if n_p != None and j_s == None:
            n_scaled = n_p/((19.2)**2)
            n_sw = n_scaled * 1e-6
            T_scaled = T_sw/((19.2)**0.5) #Richardson paper for temp scaling
            T_sheath = T_scaled
            v_bulk = v_sw
            n_q = n_q * (1/n_q.max()) * n_sw
        if v_sf != None:
            v_bulk = v_sf
            v_sw = v_bulk
            T_sheath = 5.45e4
        else:
            T_sheath = 5.45e4
            v_bulk = 450e3
            v_sw = v_bulk

        abundance_slow = 1.48E-5 #from Whittaker and Sembay (2016)
        abundance_fast = 6.69E-6
        if v_sw > 500e3:
            n_n = n_n * abundance_fast #magnetosheath ion density = 0.1*abundance
        else:
            n_n = n_n * abundance_slow
        v_therm = np.sqrt((3*constants.Boltzmann*T_sheath)/constants.m_p) 
        v_rel = (np.sqrt(v_bulk**2 + v_therm**2))*(1e2)
        sigma_sqn_slow = (1/3)*((34+10+11+1.3+0.79+1.3+0.06)*(1e-16)) + (2/3)*(12e-15)
        sigma_sqn_fast = (1/3)*((32+9.9+11+1.2+1.2+0.68+0.02)*(1e-16)) + (2/3)*(12e-15)
        if v_sw > 500e3:
            sigma_sqn = sigma_sqn_fast
        else:
            sigma_sqn = sigma_sqn_slow
        n_q = n_q * 4
        ver = n_n * n_q * v_rel * sigma_sqn * 1/(4*np.pi)
        return ver