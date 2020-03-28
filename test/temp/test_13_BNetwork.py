#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 march, 2020
Testing suite for Network export_to_shp function
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 26 march, 2020 (COVID19 quarantine)
"""

import unittest
from topopy import Network, BNetwork, Flow, DEM
import numpy as np
import matplotlib.pyplot as plt
infolder = "data/in"
outfolder = "data/out"

def draw_net(net, heads=None, ax = None):
    if ax is None:
        fig, ax = plt.subplots()
    row, col = net.ind_2_cell(net._ix)
    x, y = net.cell_2_xy(row, col)
    ax.plot(x, y, "b.")
    
    if heads is not None:
        row, col = net.ind_2_cell(heads)
        x, y = net.cell_2_xy(row, col)
        ax.plot(x, y, "ro")



class BNetwork_test(unittest.TestCase):
    
    def test_create_BNetwork(self):
        files = ["small25", "jebja30", "tunez"]
        for file in files:
            # Cargamos flow, dem and network 
            flw_path = "{}/{}_fd.tif".format(infolder, file)
            dem_path = "{}/{}.tif".format(infolder, file)
            net_path = "{}/{}_network.net".format(infolder, file)
            dem = DEM(dem_path)
            fd = Flow(flw_path)
            net = Network(net_path)
            
            # Generamos 50 puntos aleatorios (cabeceras)
            xi = np.random.randint(dem.get_extent()[0], dem.get_extent()[1], 50)
            yi = np.random.randint(dem.get_extent()[2], dem.get_extent()[3], 50)
            raw_heads = np.array((xi, yi)).T
            # Guardamos puntos (para comprobaciones)
            np.savetxt("{}/{}_rdn-heads".format(outfolder, file), raw_heads, comments="",
                       header="X;Y", delimiter=";", encoding="utf-8")
            
            # Obtenemos 5 cuencas aleatorias
            cuencas = fd.get_drainage_basins(min_area = 0.0025)
            bids = np.unique(cuencas.read_array())
            if 0 in bids:
                bids = bids[1:]
            bids = np.random.choice(bids, 5)
            
            
            # Generamos BNetwork
            for bid in bids:
                bnet = BNetwork(net, cuencas, raw_heads, bid)
                bnet_path = "{}/{}_Bnetwork.net".format(outfolder, file)
                bnet.save(bnet_path)
                print(file, bid, bnet._heads)
                draw_net(bnet, bnet._heads)
                
            
            

if __name__ == "__main__":
    unittest.main()