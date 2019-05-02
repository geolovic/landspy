#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 25 april, 2019
Testing suite for Network class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 25 april, 2019
"""

import unittest
import numpy as np
from topopy import Grid, Flow, Network
infolder = "data/in"
outfolder = "data/out"

class NetworkClassTest(unittest.TestCase):
    
    def test_create_network_zx_ax(self):
        """
        Test the creation of self._zx and self._ax arrays of Network object
        """
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            # Cargamos dem, flow direction y flow accumulation
            dem_path = infolder + "/{0}.tif".format(file)
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            fac_path = infolder + "/{0}_fac.tif".format(file)
            dem = Grid(dem_path)
            flw = Flow(flw_path)
            fac = Grid(fac_path)
            
            # Creamos Network object with default threshold
            net = Network(flw)
            
            # Seleccionamos 25 puntos al azar
            ptos = np.random.randint(0, net._ix.size, 25)
            
            # Tomamos los id de las celdas, sus alturas y el área
            cell_id = net._ix[ptos]
            cell_zx = net._zx[ptos]
            cell_ax = net._ax[ptos]
            
            # Cogemos los valores de altura y flow accumulation del raster
            row, col = dem.ind_2_cell(cell_id)
            zx = dem.get_value(row, col).astype(np.float)
            ax = fac.get_value(row, col) * fac.get_cellsize()[0] * fac.get_cellsize()[1]  * -1
                      
            res = (np.array_equal(zx, cell_zx), np.array_equal(ax, cell_ax))
            self.assertEqual(res, (True, True))
        
    def test_save_load_without_gradients(self):
        files = ["small25", "morocco", "tunez", "jebja30"]
        for file in files:
            # Cargamos dem, flow direction y flow accumulation
            flw_path = infolder +  "/{0}_fd.tif".format(file)
            net_path = outfolder +  "/{0}.net".format(file)
            flw = Flow(flw_path)
            
            net = Network(flw)
            net.save(net_path)
            net2 = Network(net_path)
            
            res = (np.array_equal(net._ix, net2._ix), np.array_equal(net._ixc, net2._ixc))
            self.assertEqual(res, (True, True))
            
            props_01 = (net._ncells, net._proj, net._size, net._cellsize, net._geot)
            props_02 = (net2._ncells, net2._proj, net2._size, net2._cellsize, net2._geot)
            self.assertEqual(props_01, props_02)


if __name__ == "__main__":
    unittest.main()  
