# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:06:28 2019

@author: Usuario
"""
from topopy import DEM, Flow, Network, BNetwork


fld = Flow("../data/in/morocco_fd.tif")
net = Network(fld, gradients=True)

basins = fld.get_drainage_basins()
bnet = BNetwork(net, basin=basins, bid=4)
