# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 11:45:01 2021

@author: Usuario
"""

from topopy import Network, rivers_to_channels, Channel
import numpy as np
import os


path = "../test/data/in/jebja_channels.shp"
net = Network("../test/data/in/jebja30_net.dat")
field_n = "segid"

canales = rivers_to_channels(path, net, field_n)
canales = np.array(canales)

canal = canales[0]


np.save("canales.npy", canales)