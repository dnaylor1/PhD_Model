    def plot_sheath(self,points_bs,points_mp,x_points,y_points,z_points,x_sheath,y_sheath,z_sheath,x_mp,y_mp,z_mp,x_bs,y_bs,z_bs,density_grid):
        """
        Plot the magnetosheath region

        Parameters:
        
        """
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(points_bs[:, 0], points_bs[:, 1], points_bs[:, 2], color='blue', alpha=0.1, label='points_bs')
        #ax.scatter(points_mp[:, 0], points_mp[:, 1], points_mp[:, 2], color='red', alpha=0.1, label='points_mp')
        #ax.set_xlabel(r'$x$ ($R_{U}$)')
        #ax.set_ylabel(r'$y$ ($R_{U}$)')
        #ax.set_zlabel(r'$z$ ($R_{U}$)')

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_points,y_points,z_points,color='green',alpha=0.05)
        #ax.plot_surface(x_mp, y_mp, z_mp, color='red', edgecolor='none', alpha=0.05)
        #ax.plot_surface(x_bs, y_bs, z_bs, color='blue', edgecolor='none', alpha=0.05)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        ax.set_zlabel(r'$z$ ($R_{U}$)')
        #plt.savefig("Sheath maybe working",dpi=1200)

        z_slice = 0
        # Indices where the grid intersects the slicing plane
        slice_indices = np.abs(self.Z_grid - z_slice) < 1e-2 
        ## grid points may not align exactly with the slicing plane due to interpolation so this makes sure they are close enough
        # Data points for the slicing plane
        x_slice = self.X_grid[slice_indices]
        y_slice = self.Y_grid[slice_indices]
        density_slice = density_grid[slice_indices]

        fig, ax = plt.subplots()
        scatter = ax.scatter(x_slice, y_slice, c=density_slice, cmap='Blues', alpha=0.7)

        from Surface import Surface
        bs = Surface(20,0.88)
        mp = Surface(16,0.6)
        x_bow, y_bow = bs.surf_2D()
        x_mag, y_mag = mp.surf_2D()
        
        ax.plot(x_bow,y_bow,color='red')
        ax.plot(x_mag,y_mag,color='blue')
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$y$ ($R_{U}$)')
        fig.colorbar(scatter, label='Density')
        #plt.title(f'Slice through region at z = {z_slice}')
        plt.title(r"Magnetosheath $x$-$y$ Projection")
        plt.savefig("High res sheath x-y",dpi=1200)

        y_slice = 0
        # Indices where the grid intersects the slicing plane
        slice_indices_xz = np.abs(self.Y_grid - y_slice) < 1e-2 
        ## grid points may not align exactly with the slicing plane due to interpolation so this makes sure they are close enough
        # Data points for the slicing plane
        x_slice_xz = self.X_grid[slice_indices_xz]
        z_slice_xz = self.Z_grid[slice_indices_xz]
        density_slice_xz = density_grid[slice_indices_xz]
        fig = plt.figure()
        ax = plt.gca()
        scatter = ax.scatter(x_slice_xz, z_slice_xz, c=density_slice_xz, cmap='Blues')
        #ax.plot(x_bow,y_bow,color='red')
        #ax.plot(x_mag,y_mag,color='blue')
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_xlabel(r'$x$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        fig.colorbar(scatter, label='Density')
        #plt.title(f'Slice through region at z = {z_slice}')
        plt.title(r"Magnetosheath $x$-$z$ Projection")
        #plt.savefig("GGG x-y",dpi=1200)

        print(f"x_slice_xz range: [{np.min(x_slice_xz)}, {np.max(x_slice_xz)}]")
        print(f"z_slice_xz range: [{np.min(z_slice_xz)}, {np.max(z_slice_xz)}]")

        print(f"x_slice_xy range: [{np.min(x_slice)}, {np.max(x_slice)}]")
        print(f"y_slice_xy range: [{np.min(y_slice)}, {np.max(y_slice)}]")

        fig = plt.figure()
        ax = plt.gca()
        y_mid = int(np.shape(self.Y_grid)[0]/2)
        xz_plane = ax.contourf(self.X_grid[:,y_mid,:],self.Z_grid[:,y_mid,:],density_grid[:,y_mid,:],cmap='plasma')
        fig.colorbar(xz_plane)

        fig = plt.figure()
        ax = plt.gca()
        x_mid = int(np.shape(self.X_grid)[0]/2)
        yz_plane = ax.contourf(self.Y_grid[x_mid,:,:],self.Z_grid[x_mid,:,:],density_grid[x_mid,:,:],cmap='Blues')
        fig.colorbar(yz_plane)
        ax.set_xlim(-80,80)
        ax.set_ylim(-80,80)
        ax.set_xlabel(r'$y$ ($R_{U}$)')
        ax.set_ylabel(r'$z$ ($R_{U}$)')
        plt.title(r"Magnetosheath $y$-$z$ Projection")
        #plt.savefig("Magnetosheath y-z Projection Low Res",dpi=1200)
        plt.savefig("High res sheath y-z",dpi=1200)

        fig = plt.figure()
        ax = plt.gca()
        z_mid = int(np.shape(self.Z_grid)[0]/2)
        xy_plane = ax.contourf(self.X_grid[:,:,z_mid],self.Y_grid[:,:,z_mid],density_grid[:,:,z_mid],cmap='plasma')
        fig.colorbar(xy_plane)
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(x_sheath,y_sheath,z_sheath,alpha=0.5)
        #ax.plot_surface(x_mp, y_mp, z_mp, color='red', edgecolor='none', alpha=0.2, label="Magnetopause")
        #ax.plot_surface(x_bs, y_bs, z_bs, color='blue', edgecolor='none', alpha=0.2, label="Bow Shock")
        #ax.set_xlabel(r'$x$ ($R_{U}$)')
        #ax.set_ylabel(r'$y$ ($R_{U}$)')
        #ax.set_zlabel(r'$z$ ($R_{U}$)')