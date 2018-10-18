# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

from topopy import DEM, Flow, Network

infolder = "data/in"
outfolder = "data/out"

files = ["small25", "morocco", "tunez", "jebja30"]
for file in files:
    flw_path = infolder +  "/{0}_fd.tif".format(file)
    dem_path = infolder +  "/{0}.tif".format(file)
    dem = DEM(dem_path)
    fd = Flow(flw_path)
    thr = int(fd.get_ncells() * 0.001)
    # Simplemente probamos que no hay fallos al crear el objeto net
    net = Network(dem, fd, thr)
    net.save(infolder + "/{0}_network.net".format(file))