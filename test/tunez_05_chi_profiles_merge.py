# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import TProfile, BNetwork, canal_2_profile
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import ogr, osr
import os



# GET MEDJERDA PROFILES
basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks"
bnet = BNetwork(basedir + "/bnet_medjerda.net")
bnet.calculate_chi(0.3)
canales = bnet.get_channels(25)

med_perfiles = []
rids = [0, 1, 2, 9, 10] # Right IDs for perfiles
rids.extend(range(14, 34))

for n, canal in enumerate(canales):
    perfil = canal_2_profile(canal, bnet, rids[n])
    med_perfiles.append(perfil)
    

# GET MAIN PROFILES
basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks"
nbasins = 33
main_perfiles = []
for n in range(1, nbasins +1):
    path = basedir + "/bnet_{:02}.net".format(n)
    bnet = BNetwork(path)
    bnet.calculate_chi(0.3)
    canal = bnet.get_main_channel()
    main_perfiles.append(canal_2_profile(canal, bnet, n))
   
out_perfiles = []

# GET OUTER PROFILES
for n in [3, 4, 5, 6, 7, 8, 11, 12, 13]:
    perfil = main_perfiles[n-1]
    perfil.calculate_chi()
    out_perfiles.append(main_perfiles[n-1])
    
    
# SAVE PROFILES ARRAYS
np.save(basedir + "/perfiles.npy", main_perfiles)
np.save(basedir + "/perfiles_med.npy", med_perfiles)
np.save(basedir + "/perfiles_out.npy", out_perfiles)

fig, ax = plt.subplots()
ramp = plt.get_cmap("viridis")
med = med_perfiles[0]
for n, perfil in enumerate(med_perfiles[1:]):
    color = ramp(n/len(med_perfiles[1:]))
    ax.plot(perfil.get_chi(), perfil.get_z(), label=str(perfil.rid), c=color, lw=0.9)

ramp = plt.get_cmap("cool")
for n, perfil in enumerate(out_perfiles):
    color =  color = ramp(n/len(out_perfiles))
    ax.plot(perfil.get_chi(relative=True), perfil.get_z(relative=True), label=str(perfil.rid), c=color, lw=0.9)
    
ax.plot(med.get_chi(), med.get_z(), label=str(med.rid), c="k", lw=1.2)
ax.set_xlim(xmin=0)
ax.set_ylim(ymin=0)
ax.set_xlabel("Chi [m]")
ax.set_ylabel("Elevation [m]")
