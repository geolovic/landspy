# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import DEM, Flow, Network, BNetwork, extract_points, Channel
import matplotlib.pyplot as plt
import numpy as np
from profiler import TProfile

# Test longitudinal profile
# Load profiles (TProfile)
prof_path = "/Users/vicen/Desktop/genil_profiler/scripts/files/graph_profiles.npy"
perfiles = np.load(prof_path, allow_pickle=True)
perfil = perfiles[0]

# Load channels (Channel)
net = Network("/Users/vicen/Desktop/genil_profiler/gisdata/genil_net.dat")
net.calculate_gradients(20, "slp")
heads = extract_points("/Users/vicen/Desktop/genil_profiler/gisdata/heads.shp", "id")
canal = net.get_channel(heads[0])

fig, ax = plt.subplots()

# di = perfil.get_l()
# zi = perfil.get_z()
# ax.plot(di, zi, c="b")

# di2 = canal.get_d()
# zi2 = canal.get_z()
# ax.plot(di2, zi2, c="r")


slopes = perfil.get_slope()
areas = perfil.get_area()
ax.plot(areas, slopes, "b.", mew=0.5)

slopes = canal.get_slope()
areas = canal.get_a(cells=False)
ax.plot(areas, slopes, "r.", mew=0.5)

ax.set_xlabel("Area $m^2$")
ax.set_ylabel("Slope (reg. points: {0})".format(perfil.slope_reg_points))
ax.set_xscale("log")
ax.set_yscale("log")
