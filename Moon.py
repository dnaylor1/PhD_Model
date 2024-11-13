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