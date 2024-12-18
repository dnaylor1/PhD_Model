# Soft X-Ray Emission from Uranus's Magnetosheath

This models simulates charge exchange-driven soft x-ray emission within Uranus's magnetosheath. It considers a simple bullet-shaped magnetopause and bow shock and a much simplified configuration to the reality at Uranus. It is the first step in developing an detailed model of emission to determine the viability of sending an SXI along with the planned Uranus flagship mission.

## Model

This is the main file for the project and includes settings for seasonal configuration, solar wind variations and which figures are to be plotted. Variations are set to None by default, in which the Voyager 2, Toth et al. (2004) conditions are used.

## Plotter

This class contains all the code to plot the different figures for the project, including neutral densities, magnetopause and bow shock surfaces, volumetric emission rate (in all planes), flux, and a combination of VER and flux. It generates the figures which are then displayed upon the model returning to Model file and reaching plt.show() at the end of the code.

## System

This file contains the System class, containing methods that are related to the system as a whole. For example, there are methods to set up the grid and to calculate the total neutral density in the system (given the exopheric density and each moon's density at each point). 

## Surface

This contains the Surface class, which has methods to set up the magnetopause and bow shock surfaces, both in 3D (for the general model) and 2D (if the surfaces need to be visualised on a 2D slice).

## Magnetosheath

This contains the Magnetosheath class which simply defines the sheath region by interpolating the magnetopause and bow shock surfaces and looking for the region inside the bow shock, but outside the magnetopause. It also calculates volumetric emission from the region, given by
$$P = \sum_{n} n_{n}n_{q}v_{\mathrm{rel}}\sigma{sqn}b_{sqj} / 4\pi$$
Where $n_{n}$ is neutral density (with $\sum_{n}$ representing the summation over the different neutral species), $n_{q}$ is the magnetosheath ion density, $v_{\mathrm{rel}}$ is the relative velocity between the ions and neutrals, $\sigma_{sqn}$ is the interaction cross section, and $b_{sqj}$ is the branching ratio of the interaction (the probablility for soft x-ray emission to occur). The cross sections and branching ratios are combined into one $\sigma$ term using the work of Bodewits et al. (2007). Volumetric emission is given in photon cm$^{-3}$ s$^{-1}$.

## SXI

This contains the SXI class which can be instantiated as SMILE, LEXI etc. This then uses the configurations of the specific SXI to determine how far from the magnetosheath the SXI must be to image the system on approach. It takes VER and calculates flux and integration times from it.

## Moon

This file contains the Moon parent class and subclasses for each Moon of the system (Miranda, Ariel, Umbriel, Titania, Oberon), defining its distance from Uranus, its density (Cheng, ER Max, ER Min) and scale height, and the radius of the tori. It has a method to calculate the density contribution of each moon.