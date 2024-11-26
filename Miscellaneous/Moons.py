from Moon import *
import json

#r_min, r_max, density, delta_Z, ER_min_H, ER_min_O, ER_max_H, ER_max_O

file = open('config.json')
moons = json.load(file)

def moon_params(moon_name):
    """
    Retrieves each moon's parameters from the config file and returns a dictionary of params.

    Parameters
        moon_name(str): name of the moon for which parameters need retrieving. Will match a string in the config file

    Returns
        params(dict): dictionary of parameters and the name of those parameters, e.g. r_max = "22" etc.
    """
    for i in moons["Moons"]:
        if i["moon"] == moon_name:
            params = i["params"]
            return params
    return None
    
Miranda_params = moon_params("Miranda")
Miranda = Moon(**Miranda_params)
Ariel_params = moon_params("Ariel")
Ariel = Moon(**Ariel_params)
Umbriel_params = moon_params("Umbriel")
Umbriel = Moon(**Umbriel_params)
Titania_params = moon_params("Titania")
Titania = Moon(**Titania_params)
Oberon_params = moon_params("Oberon")
Oberon = Moon(**Oberon_params)

moons = [Miranda, Ariel, Umbriel, Titania, Oberon]

#print("Moon parameters:", ", ".join(f"{key}={value}" for key, value in Miranda_params.items()))

""" Miranda = Moon(2.9, 9.9, 0.06, 0.38, 0.59, 0.57, 13.8, 11.68)
Ariel = Moon(3.8, 17, 0.11, 0.68, 4.37, 3.97, 28.4, 38.74)
Umbriel = Moon(4.7, 29, 0.12, 1.1, 4.06, 4.11, 15.5, 24.67)
Titania = Moon(6.1, 73, 0.02, 2.3, 22.3, 14.5, 38.1, 34.4)
Oberon = Moon(6.9, 150, 0.002, 3.6, 31, 19.5, 39.2, 30.6) """
