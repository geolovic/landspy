#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:10:26 2019

@author: vicen
"""

from topopy import DEM, Flow, Network
import numpy as np

dem = DEM("data/in/morocco.tif")
flw = Flow(dem)

net = Network(flw, 1000)

net.export_2_points("data/out/prueba_morocco.txt")