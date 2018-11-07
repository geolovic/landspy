#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:53:42 2018

@author: vicen
"""

from topopy import Network
import numpy as np
import matplotlib.pyplot as plt
import ogr, osr
import os

dem = DEM("../data/in/morocco.tif")
fd = Flow(dem)
net = Network(dem, fd, 1000)


basin = fd.get_drainage_basins([579223, 504380], asgrid=False)
