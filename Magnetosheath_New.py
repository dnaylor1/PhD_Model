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

        #density_grid[np.isfinite(interpolated_bs) & np.isfinite(interpolated_mp)] = True
        density_grid = np.zeros_like(self.rad)
        """ density_grid[np.isfinite(interpolated_bs)] = True
        density_grid[np.isfinite(interpolated_mp)] = False """
        """ for i,j,k in range(0,np.shape(density_grid)[0]):
            if density_grid == True:
                density_grid[i,j,k] = 0.1 """

        density_grid[np.isfinite(griddata(points_bs, np.ones(len(points_bs)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear')) & 
                    ~np.isfinite(griddata(points_mp, np.ones(len(points_mp)), (self.X_grid, self.Y_grid, self.Z_grid), method='linear'))] = True

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

        return points_mp,points_bs,x_points,y_points,z_points,x_sheath,y_sheath,z_sheath,density_grid