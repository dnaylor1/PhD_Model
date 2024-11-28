"""     def sheath(self, r_mp_grid, r_bs_grid, x_mp, y_mp, z_mp, x_bs, y_bs, z_bs):
        valid = ~np.isnan(r_mp_grid) & ~np.isnan(r_bs_grid)
        magnetosheath_region = valid & (r_mp_grid <= self.rad) & (self.rad <= r_bs_grid)
        density = np.zeros_like(self.rad)  
        density[magnetosheath_region] = 1  
        ##extract the magnetosheath points 
        x_sheath = self.xbox[magnetosheath_region]
        y_sheath = self.ybox[magnetosheath_region]
        z_sheath = self.zbox[magnetosheath_region]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_sheath,y_sheath,z_sheath)
        ax.plot_surface(x_mp, y_mp, z_mp, cmap='plasma', edgecolor='none', alpha=0.2)
        ax.plot_surface(x_bs, y_bs, z_bs, cmap='viridis', edgecolor='none', alpha=0.2)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        ax.set_zlabel(r'$z$ ($R_{U}$)')
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_zlim(-80,80) """ ### from system class