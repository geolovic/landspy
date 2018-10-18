# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""
import numpy as np
from topopy import DEM, Flow, Network
import matplotlib.pyplot as plt
import gdal, ogr, osr

dem = DEM("../data/in/jebja30.tif")
fd = Flow("../data/in/jebja30_fd.tif")
net = Network(dem, fd, 2000)
path = "../data/out/jebja_streams.shp"

net.export_to_shapefile(path)