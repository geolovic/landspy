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

basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks/out_files"

perfiles = np.load(basedir + "/perfiles.npy")

# GET MEDJERDA PROFILES
med_perfiles = perfiles[1:27]
out_perfiles = perfiles[27:]
med = perfiles[0]

fig, ax = plt.subplots()
ramp = plt.get_cmap("viridis")
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