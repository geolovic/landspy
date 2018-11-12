# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:36:31 2018

@author: Usuario
"""

from topopy import Network, BNetwork, Grid
import numpy as np
import ogr
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Get network and basins
basedir = "C:/Users/Usuario/Desktop/tunez/gisdata"
bnet = BNetwork(basedir + "/chi_analysis/bnet_05.net")

canales = bnet.get_channels()


