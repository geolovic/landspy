#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08 October, 2018
Testing suite for Network get_streams functions
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 08 October, 2018
"""

from topopy import Flow, DEM, Network
infolder = "data/in"
outfolder = "data/out"

file = "jebja30"

flw_path = infolder +  "/{0}_fd.tif".format(file)
dem_path = infolder +  "/{0}.tif".format(file)
dem = DEM(dem_path)
fd = Flow(flw_path)
thr = int(fd.get_ncells() * 0.01)
# Simplemente probamos que no hay fallos al crear el objeto net
net = Network(dem, fd, thr)
net.save(outfolder + "/{0}_net.txt".format(file))
net.load(outfolder + "/{0}_net.txt".format(file))
