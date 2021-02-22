#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 february, 2021
Testing suite for BNetwork class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 09 february, 2021
"""

import unittest
import os
import numpy as np
from topopy import Flow, Basin, Network, BNetwork, DEM
from topopy.network import NetworkError
infolder = "data/in"
outfolder = "data/out"


class BNetworkClassTest(unittest.TestCase):
    
    # Indices de las cabeceras que deben de salir (para comprobar)
    results = {"small25":dict([(1, 16171), (2, 9354),  (3,1463)]),
               "jebja30":dict([(1, 151755), (2, 44786), (3, 48709), (4, 3819)]),
               "tunez":dict([(1, 77552), (2, 30013), (3, 7247)])}
    
    def test_BNetwork_class00(self):
        """
        Test00 Crea BNetwork para cuencas de prueba a partir de un Grid de cuencas
        Sin utilizar cabeceras 
        """
        files = ["small25", "jebja30", "tunez"]
        for file in files:
            # Cargamos DEM, Flow, Network
            fd = Flow("{}/{}_fd.tif".format(infolder, file))
            net = Network("{}/{}_net.dat".format(infolder, file))
            
            # Cargamos outlets y generamos cuencas
            outlets = np.loadtxt("{}/{}_bnet_outlets.txt".format(infolder, file), delimiter=";")
            outlets = net.snap_points(outlets)
            cuencas = fd.get_drainage_basins(outlets)
            
            for bid in np.unique(cuencas.read_array()):
                if bid == 0:
                    continue
                bnet = BNetwork(net, cuencas, None, bid)
                self.assertEqual(int(bnet._heads[0]), self.results[file][bid])
                            
    def test_BNetwork_class01(self):
        """
        Test00 Crea BNetwork  para cuencas de prueba a partir de un objeto Basin
        Sin utilizar cabeceras
        """
        files = ["small25", "jebja30", "tunez"]
        for file in files:
            # Cargamos DEM, Flow, Network
            fd = Flow("{}/{}_fd.tif".format(infolder, file))
            dem = DEM("{}/{}.tif".format(infolder, file))
            net = Network("{}/{}_net.dat".format(infolder, file))
            
            # Cargamos outlets y generamos cuencas
            outlets = np.loadtxt("{}/{}_bnet_outlets.txt".format(infolder, file), delimiter=";")
            outlets = net.snap_points(outlets)
            cuencas = fd.get_drainage_basins(outlets)
            
            for bid in np.unique(cuencas.read_array()):
                if bid == 0:
                    continue
                basin = Basin(dem, cuencas, bid)
                bnet = BNetwork(net, basin)
                # Este test solo verifica que se realice sin fallos y que
                # el objeto bnet tiene una única cabecera
                bnet = BNetwork(net, cuencas, None, bid)
                self.assertEqual(int(bnet._heads[0]), self.results[file][bid])    
       
    def test_BNetwork_class03(self):
        """
        Test que prueba cabeceras en cuenca 1 con small25
        474260.9;4114339.6;3
        474856.9;4114711.1;2           
        """
        # Cargamos DEM, Flow, Network
        fd = Flow("{}/{}_fd.tif".format(infolder, "small25"))
        net = Network("{}/{}_net.dat".format(infolder, "small25"))
        
        # Cargamos outlets, heads y generamos cuencas
        outlets = np.loadtxt("{}/{}_bnet_outlets.txt".format(infolder, "small25"), delimiter=";")
        heads = np.loadtxt("{}/{}_bnet_heads.txt".format(infolder, "small25"), delimiter=";")
        outlets = net.snap_points(outlets)
        cuencas = fd.get_drainage_basins(outlets)
        
        bid = 1
        bnet = BNetwork(net, cuencas, heads, bid)
        self.assertEqual(np.array_equal(bnet._heads, np.array([13494, 16171])), True)

    def test_BNetwork_class04(self):
        """
        Test que prueba cabeceras en cuenca 1 con small25 (sin utilizar id field)
        474260.9;4114339.6;3
        474856.9;4114711.1;2           
        """
        # Cargamos DEM, Flow, Network
        fd = Flow("{}/{}_fd.tif".format(infolder, "small25"))
        net = Network("{}/{}_net.dat".format(infolder, "small25"))
        
        # Cargamos outlets, heads y generamos cuencas
        outlets = np.loadtxt("{}/{}_bnet_outlets.txt".format(infolder, "small25"), delimiter=";")
        heads = np.loadtxt("{}/{}_bnet_heads.txt".format(infolder, "small25"), delimiter=";")
        # Remove the id column
        heads = heads[:,:-1]
        outlets = net.snap_points(outlets)
        cuencas = fd.get_drainage_basins(outlets)
        
        bid = 1
        bnet = BNetwork(net, cuencas, heads, bid)
        self.assertEqual(np.array_equal(bnet._heads, np.array([16171, 13494])), True)
               
    def test_BNetwork_class05(self):
        """
        Test de creado masivo de cuencas con cabeceras aleatorias
        """

        files = ["small25", "jebja30", "tunez"]
        for file in files:
            # Cargamos DEM, Flow, Network
            fd = Flow("{}/{}_fd.tif".format(infolder, file))
            net = Network("{}/{}_net.dat".format(infolder, file))
            dem = DEM("{}/{}.tif".format(infolder, file))
            
            # Generamos todas las cuencas
            cuencas = fd.get_drainage_basins(min_area = 0.0025)
            
            # Generamos 50 puntos aleatorios dentro de la extensión del objeto Network
            # Estos 50 puntos se usaran como cabeceras
            xmin, xmax, ymin, ymax = net.get_extent()
            xi = np.random.randint(xmin, xmax, 50)
            yi = np.random.randint(ymin, ymax, 50)
            heads = np.array((xi, yi)).T
            
            # Cogemos 5 cuencas aleatorias
            bids = np.random.choice(np.unique(cuencas.read_array())[1:], 5)
            for bid in bids:
                try:
                    if np.random.randint(100) < 70:
                        bnet = BNetwork(net, cuencas, heads, bid)
                    else:
                        basin = Basin(dem, cuencas, bid)
                        bnet = BNetwork(net, basin, heads)
                except NetworkError:
                    print("Network of {} file inside the basin {} has not enough pixels".format(file, bid))
                    continue
                
                # Salvamos BNetwork y volvemos a cargar para comprobar que se cargan-guardan bien
                bnet_path = "{}/{}_{}_bnet.dat".format(outfolder, file, bid)
                bnet.save(bnet_path)
                bnet2 = BNetwork(bnet_path)
                computed = np.array_equal(bnet._ix, bnet2._ix)
                self.assertEqual(computed, True)
                # borramos archivo
                os.remove(bnet_path)


if __name__ == "__main__":
    unittest.main()