import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Surface:
    def __init__(self, r0, K):
        self.r0 = r0
        self.K = K
    
    def define_surface(self):
        """
        Defines the surface for either the magnetopause or bow shock.
        """
        theta = np.linspace(-np.pi,0,100)
        exclude = np.pi
        theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
        theta = theta_exc
        phi = np.linspace(0,2*np.pi,100)
        r = self.r0 * (2 / (1 + np.cos(theta)))**self.K
        theta, phi = np.meshgrid(theta, phi)
        x = r*np.cos(theta)
        y = r*np.sin(theta)*np.cos(phi)
        z = r*np.sin(theta)*np.sin(phi)
        return x,y,z
    