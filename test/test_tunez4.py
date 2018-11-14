# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import DEM, Network, BNetwork, Grid, TProfile, canal_2_profile, ProfilerApp
import numpy as np
import ogr
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# Number of basins
nbasins = 33

# Get network and basins
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata/chi_analysis"
bnet = BNetwork(basedir + "/bnet_04.net")

heads = bnet.calculate_gradients(5)
print(heads)