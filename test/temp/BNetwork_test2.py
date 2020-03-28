#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 march, 2020 (COVID19 quarantine)
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 28 march, 2020 (COVID19 quarantine)
"""

from topopy import Network, Grid, Basin, Flow, DEM, BNetwork
import numpy as np

files = ["small25", "jebja30", "tunez"]

for file in files:
    file = "jebja30"
    # Cargamos DEM, Flow, Network
    fd = Flow("../data/in/{}_fd.tif".format(file))
    net = Network("../data/in/{}_network.net".format(file))
    dem = DEM("../data/in/{}.tif".format(file))
    
    # Generamos todas las cuencas
    cuencas = fd.get_drainage_basins(min_area = 0.0025)
    cuencas.save("../data/out/{}_testbasins.tif".format(file))
    
    # Generamos 20 puntos aleatorios
    xmin, xmax, ymin, ymax = net.get_extent()
    xi = np.random.randint(xmin, xmax, 50)
    yi = np.random.randint(ymin, ymax, 50)
    heads = np.array((xi, yi)).T
    
    np.savetxt("../data/out/{}_raw_heads.txt".format(file), heads, delimiter=";",
               header="X;Y", comments="", encoding="utf-8")
    
    # Cogemos 5 cuencas aleatorias
    bids = np.random.choice(np.unique(cuencas.read_array())[1:], 5)

    for bid in bids:
        bnet = BNetwork(net, cuencas, heads, bid)
        canales = bnet.get_streams()
        canales.save("../data/out/{}_bnet{}.tif".format(file, bid))
        row, col = bnet.ind_2_cell(bnet._heads)
        x, y = bnet.cell_2_xy(row, col)
        arr = np.array((x, y)).T
        np.savetxt("../data/out/{}_{}_heads.txt".format(file, bid), arr, delimiter=";",
               header="X;Y", comments="", encoding="utf-8")
                     
    
    
