#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:15:55 2020

@author: vicen
"""

from topopy import Grid, DEM, Flow, Network, BNetwork
import ogr
import numpy as np
import matplotlib.pyplot as plt

dem_path = "../data/in/jebja30.tif"

dem = DEM(dem_path)
fd = Flow(dem)
basins = fd.get_drainage_basins()
net = Network(fd)
canales = net.get_streams()

basins.save("C:/Users/Usuario/Desktop/temp/cuencas.tif")
canales.save("C:/Users/Usuario/Desktop/temp/canales.tif")







