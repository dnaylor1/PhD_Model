from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import UnivariateSpline
from scipy import interpolate
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

class Surface:
    def __init__(self, r0, K):
        self.r0 = r0
        self.K = K

    
    def define_surface(self, X_grid, n_p=None, v_sw=None, type=None, x_min=10, x_max=10, combd=None, v_type = None):
        """
        Defines the surface for either the magnetopause or bow shock
        
        Parameters:
            X_grid (ndarray): x grid to create a meshgrid

        Returns:
            x,y,z (Ndarray): arrays for the x,y,z coordinates of the surface
            r0, k (float): returns updated standoff distances and flaring parameters if solar wind variations are used
        """
        if n_p != None and self.K != None:
            surf_type = type
            
            n_scaled = n_p/((19.2)**2)
            #n_scaled = 0.1e6
            #T_scaled = T_sw/((19.2)**0.5) #Richardson paper for temp scaling

            R_standoff = ((2*(2.3e-5)**2)/(constants.mu_0*constants.m_p*n_scaled*(v_sw)**2))**(1/6)

            from Model import magnetopause
            diff = (R_standoff-magnetopause.r0)/magnetopause.r0

            if surf_type == "MP":
                r0 = R_standoff
            elif surf_type == "BS":
                r0 = self.r0 + (diff*self.r0)
            else:
                print("Invalid surface type, must be MP (magnetopause) or BS (bow shock)")
                exit()
            K = self.K + (diff*0.25*self.K)
            #K = self.K
        elif combd == "Y":
            if v_type == "S":
                    r0 = self.r0
                    K = self.K
                #if surf_type == "MP":
                    #r0 = self.r0
                    #K = self.K
                #elif surf_type == "BS":
                    #r0 = self.r0
                    #K = self.K
                #else:
                    #print("Invalid surface type")
                    #exit()
            elif v_type == "F":
                r0 = self.r0/1.26
                diff = (r0-self.r0)/self.r0
                K = self.K + (self.K*0.25*diff)                
                #if surf_type == "MP":
                 #   r0 = self.r0/1.26
                  #  diff = (r0-self.r0)/self.r0
                   # K = self.K + (self.K*0.25*diff)
                #elif surf_type == "BS":
                 #   from Model import magnetopause
                  #  r0 = self.r0/1.26
                   # K = self.K/(1.26*0.25)
                #else:
                #    print("Invalid surface type")
                #    exit()
            else:
                print("Invalid v_type")
                exit()

        elif self.K == None:
            s_type = type
            if s_type == "MP":
                r0 = self.r0
                diff = (r0-16)/16
                K = 0.6 + (diff*0.25*0.6)
            elif s_type == "BS":
                diff = (self.r0-16)/16
                r0 = 20 + (20*diff)
                K = 0.88 + (diff*0.25*0.6)
        else:
            r0 = self.r0
            K = self.K
        
        #print(type)
        #print(self.r0)
        #print(R_standoff)
        #print(diff)
        #print(r0)
        #print(self.K,K)

        theta_num = np.abs(x_max)+np.abs(x_min)+1
        theta = np.linspace(-np.pi,0,theta_num)
        exclude = np.pi
        theta_exc = theta[~np.isclose(np.abs(theta),exclude)]
        theta = theta_exc
        phi = np.linspace(0,2*np.pi,161)
        Phi, X = np.meshgrid(phi, X_grid[:,0,0])
        r_mp = r0 * (2 / (1 + np.cos(theta)))**K
        xmp = r_mp*np.cos(theta)
        interpfunc = interpolate.InterpolatedUnivariateSpline(xmp, r_mp)
        rmp = interpfunc(X)
        r_open = np.sqrt(rmp**2 - X**2)
        _y = r_open*np.sin(Phi)
        _z = r_open*np.cos(Phi)        
        x = X
        y = _y
        z = _z
        return x,y,z,r0,K
    
    def surf_2D(self,x_max,x_min,r0=None,K=None):
        """
        Defines the 2D magnetopause and bow shock surfaces for visualisation

        Parameters:
            x_max, x_min (float): grid limits in the x plane, used to define theta

        Returns:
            x, y (ndarray): x and y points of the magnetopause and bow shock
        """
        theta_num = np.abs(x_max)+np.abs(x_min)+1
        theta = np.linspace(-np.pi, np.pi, num=theta_num)
        exclude=np.pi
        theta_ex = theta[~np.isclose(np.abs(theta),exclude)]
        theta = theta_ex

        if r0 == None:
            r0 = self.r0
            K = self.K

        R = r0 * (2/(1+np.cos(theta)))**K
        x = R * np.cos(theta)
        y = R * np.sin(theta)

        return x,y

    