import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from Moon import *
from System import *

## pip freeze > file_name

################## plot a 2D projection too to make sure scale heights are working

Miranda = Moon(2.9, 9.9, 0.06, 0.38, 0.59, 0.57, 13.8, 11.68)
Ariel = Moon(3.8, 17, 0.11, 0.68, 4.37, 3.97, 28.4, 38.74)
Umbriel = Moon(4.7, 29, 0.12, 1.1, 4.06, 4.11, 15.5, 24.67)
Titania = Moon(6.1, 73, 0.02, 2.3, 22.3, 14.5, 38.1, 34.4)
Oberon = Moon(6.9, 150, 0.002, 3.6, 31, 19.5, 39.2, 30.6)
moons = [Miranda, Ariel, Umbriel, Titania, Oberon]

#system = System(moons,grid_limits=(-40,40,-40,40,-4,4))
system = System(moons)
total_density = system.calculate_total_density()
#system.plot_density()
