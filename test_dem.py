# -*- coding: utf-8 -*-
"""
Created on Mon May 28 13:54:33 2018

@author: Usuario
"""

from topopy import DEM, Flow, Network

dem = DEM("D:/GIS/01_PROYECTOS/Marruecos/gisdata/nmorocco30.tif")
flow = Flow(dem)
del(dem)
flow.save_gtiff("D:/GIS/01_PROYECTOS/Marruecos/gisdata/fd_nmorocco30.tif")
fac = flow.get_flow_accumulation()
fac.save("D:/GIS/01_PROYECTOS/Marruecos/gisdata/fac_nmorocco30.tif")
del(fac)
basins = flow.get_drainage_basins()
basins.save("D:/GIS/01_PROYECTOS/Marruecos/gisdata/bss_nmorocco30.tif")
del(basins)
net = Network(flow, 2500)
streams = net.get_streams()
streams.save("D:/GIS/01_PROYECTOS/Marruecos/gisdata/str_nmorocco30.tif")
del(streams)
del(flow)
del(net)
