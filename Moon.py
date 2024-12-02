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

    def add_density(self, rad, Z_grid, Z1):
        """
        Creates the neutral tori and finds the density contribution of that torus at each point in the model
        
        Parameters
            rad (ndarray): radius of each point
            zbox (ndarray): z meshgrid
            Z1 (ndarray): density contribution of current torus at each point, initialised as zero everywhere

        Returns
            Z1 (ndarray): updated density contribution array
        """
        mask_r = (rad >= self.r_min) & (rad <= self.r_max) #radial mask
        mask_z = np.abs(Z_grid) < self.delta_Z #z-direction mask
        mask = mask_r & mask_z
        #Z1[mask] = self.density
        Z1[mask] = self.ER_max_H+self.ER_max_O 
        return Z1
    
class Miranda(Moon):
    def __init__(self):
        super().__init__(r_min = 2.9,
                         r_max = 9.9,
                         density = 0.06,
                         delta_Z = 0.38,
                         ER_min_H = 0.59,
                         ER_min_O = 0.57,
                         ER_max_H = 13.8,
                         ER_max_O = 11.68)
        
class Ariel(Moon):
    def __init__(self):
        super().__init__(r_min = 3.8,
                         r_max = 17,
                         density = 0.11,
                         delta_Z = 0.68,
                         ER_min_H = 4.37,
                         ER_min_O = 3.97,
                         ER_max_H = 28.4,
                         ER_max_O = 38.74)
        
class Umbriel(Moon):
        def __init__(self):
             super().__init__(r_min = 4.7,
                              r_max = 29,
                              density = 0.12,
                              delta_Z = 1.1,
                              ER_min_H = 4.06,
                              ER_min_O = 4.11,
                              ER_max_H = 15.5,
                              ER_max_O = 24.67)

class Titania(Moon):
        def __init__(self):
             super().__init__(r_min = 6.1,
                              r_max = 73,
                              density = 0.02,
                              delta_Z = 2.3,
                              ER_min_H = 22.3,
                              ER_min_O = 14.5,
                              ER_max_H = 38.1,
                              ER_max_O = 34.4)
             
class Oberon(Moon):
        def __init__(self):
             super().__init__(r_min = 6.9,
                              r_max = 150,
                              density = 0.002,
                              delta_Z = 3.6,
                              ER_min_H = 31,
                              ER_min_O = 19.5,
                              ER_max_H = 39.2,
                              ER_max_O = 30.6)

