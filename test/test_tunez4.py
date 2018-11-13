# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Network, BNetwork, Grid, TProfile, canal_2_profile
import numpy as np
import ogr
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# Number of basins
nbasins = 3

# Get network and basins
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"
bnet = BNetwork(basedir + "/chi_analysis/bnet_05.net")
##bnet.calculate_chi(0.3)
#canales = bnet.get_channels()
#perfiles = []
#
#fig, ax = plt.subplots()
#canal = bnet.get_main_channel()
#
#ax.plot(canal[:, 3], canal[:, 7])
#
#perfil = canal_2_profile(canal, bnet)
#perfil.calculate_ksn(5)
#
#ax.plot(perfil.get_l(), perfil.get_ksn())
#canales = bnet.get_channels()
canal = bnet.get_main_channel()
ksn = canal[:, 7]
di = canal[:, 3]
length = di[0] - di[-1]
di -= di[-1]
di = length - di

perfil = canal_2_profile(canal, bnet)
perfil.calculate_ksn(25)
ksn2 = perfil.get_ksn()
di2 = perfil.get_l()

fig, ax = plt.subplots()
ax.plot(di, ksn)
ax.plot(di2, ksn2)