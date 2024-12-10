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

    def volumetric_emission(self,n_n,n_q,n_p=None,T_sw=None,v_sw=None,v_sf=None):
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
        if n_p != None:
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
            v_bulk = 400e3
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
        ver = n_n * n_q * v_rel * sigma_sqn * 1/(4*np.pi)
        return ver
    


        
    




