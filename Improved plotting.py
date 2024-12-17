import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

def plot_combined_ver_flux(self, ver_slow, ver_fast, flux_slow, flux_fast, int_s_s, int_h_s, int_s_f, int_h_f):
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 4, width_ratios=[1, 1, 1, 1], wspace=0.6)

    # Create axes for each subplot
    axes = [fig.add_subplot(gs[i // 4, i % 4]) for i in range(8)]
    z_pos = np.shape(self.Z_grid)[0] // 2
    x_pos = np.shape(self.Y_grid)[0] // 2
    y_pos = np.shape(self.Y_grid)[0] // 2

    # Plot slow VER and flux
    axes[0].contourf(self.Y_grid[x_pos, :, :], self.Z_grid[x_pos, :, :], ver_slow[x_pos, :, :], cmap='YlOrRd', levels=20)
    axes[1].contourf(self.X_grid[:, :, z_pos], self.Y_grid[:, :, z_pos], ver_slow[:, :, z_pos], cmap='YlOrRd', levels=20)
    axes[2].contourf(self.X_grid[:, y_pos, :], self.Z_grid[:, y_pos, :], ver_slow[:, y_pos, :], cmap='YlOrRd', levels=20)
    axes[3].contourf(self.Y_grid[x_pos, :, :], self.Z_grid[x_pos, :, :], flux_slow, cmap='plasma', levels=20)

    # Plot fast VER and flux
    axes[4].contourf(self.Y_grid[x_pos, :, :], self.Z_grid[x_pos, :, :], ver_fast[x_pos, :, :], cmap='YlOrRd', levels=20)
    axes[5].contourf(self.X_grid[:, :, z_pos], self.Y_grid[:, :, z_pos], ver_fast[:, :, z_pos], cmap='YlOrRd', levels=20)
    axes[6].contourf(self.X_grid[:, y_pos, :], self.Z_grid[:, y_pos, :], ver_fast[:, y_pos, :], cmap='YlOrRd', levels=20)
    axes[7].contourf(self.Y_grid[x_pos, :, :], self.Z_grid[x_pos, :, :], flux_fast, cmap='plasma', levels=20)

    # Set titles
    titles = [
        r"(a) $y$-$z$ plane (Slow)", r"(b) $x$-$y$ plane (Slow)", r"(c) $x$-$z$ plane (Slow)", r"(d) Flux (Slow)",
        r"(e) $y$-$z$ plane (Fast)", r"(f) $x$-$y$ plane (Fast)", r"(g) $x$-$z$ plane (Fast)", r"(h) Flux (Fast)"
    ]
    for ax, title in zip(axes, titles):
        ax.set_title(title)

    # Configure axes
    for ax, (xlabel, ylabel) in zip(axes, [
        (r'$y$ ($R_{U}$)', r'$z$ ($R_{U}$)'), (r'$x$ ($R_{U}$)', r'$y$ ($R_{U}$)'),
        (r'$x$ ($R_{U}$)', r'$z$ ($R_{U}$)'), (r'$y$ ($R_{U}$)', r'$z$ ($R_{U}$)'),
        (r'$y$ ($R_{U}$)', r'$z$ ($R_{U}$)'), (r'$x$ ($R_{U}$)', r'$y$ ($R_{U}$)'),
        (r'$x$ ($R_{U}$)', r'$z$ ($R_{U}$)'), (r'$y$ ($R_{U}$)', r'$z$ ($R_{U}$)')
    ]):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_aspect('equal')
        ax.set_xlim(self.x_min, self.x_max)
        ax.set_ylim(self.y_min, self.y_max)

    # Add colorbars
    cbar_axes = [
        fig.add_axes([0.93, 0.57, 0.01, 0.275]), fig.add_axes([0.04, 0.57, 0.01, 0.275]),
        fig.add_axes([0.93, 0.14, 0.01, 0.275]), fig.add_axes([0.04, 0.14, 0.01, 0.275])
    ]
    fig.colorbar(axes[3].collections[0], cax=cbar_axes[0], label=r"Flux (photon cm$^{-2}$ s$^{-1}$)")
    fig.colorbar(axes[7].collections[0], cax=cbar_axes[2], label=r"Flux (photon cm$^{-2}$ s$^{-1}$)")

    if self.config == "E":
        fig.colorbar(axes[2].collections[0], cax=cbar_axes[1], label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)")
        fig.colorbar(axes[6].collections[0], cax=cbar_axes[3], label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)")
    else:
        fig.colorbar(axes[0].collections[0], cax=cbar_axes[1], label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)")
        fig.colorbar(axes[4].collections[0], cax=cbar_axes[3], label=r"Volumetric Emission (photon cm$^{-3}$ s$^{-1}$)")

    for cbar_ax in cbar_axes[1:]:
        cbar_ax.yaxis.set_label_position('left')

    # Add global title
    title_prefix = "Equinox" if self.config == "E" else "Solstice"
    plt.suptitle(
        fr"VER and Flux: {title_prefix}, $v_{{\mathrm{{SW, Slow}}}}=400$ km s$^{{-1}}$, $v_{{\mathrm{{SW, Fast}}}}=800$ km s$^{{-1}}$, $x$, $y$, $z$ Slice Position = {x_pos},{y_pos},{z_pos}",
        x=0.5, y=0.9
    )

    # Add annotations
    plt.gcf().text(0.05, 0.53, fr"Max VER (Slow): {ver_slow.max():.3g} photon cm$^{{-3}}$ s$^{{-1}}$", ha='left', fontsize=12)
    plt.gcf().text(0.05, 0.50, fr"Mean VER (Slow): {ver_slow.mean():.3g} photon cm$^{{-3}}$ s$^{{-1}}$", ha='left', fontsize=12)
    plt.gcf().text(0.05, 0.08, fr"Max VER (Fast): {ver_fast.max():.3g} photon cm$^{{-3}}$ s$^{{-1}}$", ha='left', fontsize=12)
    plt.gcf().text(0.05, 0.05, fr"Mean VER (Fast): {ver_fast.mean():.3g} photon cm$^{{-3}}$ s$^{{-1}}$", ha='left', fontsize=12)
    plt.gcf().text(0.75, 0.53, f"Integration time (Slow) (s): {int_s_s} s", ha='left', fontsize=12)
    plt.gcf().text(0.75, 0.50, f"Integration time (Slow) (h): {int_h_s} h", ha='left', fontsize=12)
    plt.gcf().text(0.75, 0.08, f"Integration time (Fast) (s): {int_s_f} s", ha='left', fontsize=12)
    plt.gcf().text(0.75, 0.05, f"Integration time (Fast) (h): {int_h_f} h", ha='left', fontsize=12)

    # Save figure
    filename = "ver_flux_combined_equinox" if self.config == "E" else "ver_flux_combined_solstice"
    plt.savefig(filename, dpi=1200)