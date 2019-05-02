# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


basedir = "/Users/vicen/Desktop/tunez/gisdata/BNetworks/out_files"
perfiles = np.load(basedir + "/perfiles.npy")

# GET MEDJERDA PROFILES
med_perfiles = perfiles[1:25]
out_perfiles = perfiles[25:]
med = perfiles[0]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

width = 15
height = 10
sep = 0

x1 = 10
x2 = x1 + width
y1 = 10
y2 = y1 + height

ramp = plt.get_cmap("viridis")
for n, perfil in enumerate(med_perfiles):
    color = ramp(n/len(med_perfiles[1:]))
    data = np.array([[x1, x1, x2, x2], [y1, y2, y2, y1]]).T
    polygon = mpatches.Polygon(data, facecolor=color, linewidth=0.7, edgecolor="k")
    ax1.add_patch(polygon)
    y1 = y2 + sep
    y2 = y2 + height + sep
    ax2.plot(perfil.get_chi(), perfil.get_z(), label=str(perfil.rid), c=color, lw=0.9)    

ax1.set_xlim(0, 60)
ax1.set_ylim(0, 500)



x1 = 40
x2 = x1 + width
y1 = 10
y2 = y1 + height

ramp = plt.get_cmap("cool")
for n, perfil in enumerate(out_perfiles):
    color =  color = ramp(n/len(out_perfiles))
    ax2.plot(perfil.get_chi(relative=True), perfil.get_z(relative=True), label=str(perfil.rid), c=color, lw=0.9)
    data = np.array([[x1, x1, x2, x2], [y1, y2, y2, y1]]).T
    polygon = mpatches.Polygon(data, facecolor=color, linewidth=0.7, edgecolor="k")
    ax1.add_patch(polygon)
    y1 = y2 + sep
    y2 = y2 + height + sep


ax2.plot(med.get_chi(), med.get_z(), label=str(med.rid), c="k", lw=1.2)
ax2.set_xlim(xmin=0)
ax2.set_ylim(ymin=0)
ax2.set_xlabel("Chi [m]")
ax2.set_ylabel("Elevation [m]")