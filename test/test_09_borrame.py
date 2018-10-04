#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 02 October, 2018
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 04 October, 2018
"""

from topopy import Flow, DEM, Network
infolder = "data/in"
outfolder = "data/out"

dem = DEM("data/in/jebja30.tif")
fd = Flow("data/in/jebja30_fd.tif")
thr = int(fd.get_ncells() * 0.01)
# Simplemente probamos que no hay fallos al crear el objeto net
net = Network(dem, fd, 1000)
streams = net.get_stream_segments()
streams.plot()

#net.export_2_points(outfolder + "data/out/jebja30.txt")
