# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import DEM, Flow, Network, BNetwork

infolder = "/Users/vicen/Desktop/genil_profiler/gisdata"

dem = DEM(infolder + "/genil10.tif")
fd = Flow(dem)
fd.save(infolder + "/genil10_fd.tif")
fac = fd.get_flow_accumulation()
fac.save(infolder + "/genil10_fac.tif")

