#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 february, 2021
Testing suite for BNetwork class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@last_modified: 15 october, 2025
"""

import unittest
import numpy as np
from landspy import Network, Channel, shp_to_channels
from osgeo import ogr

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"


class ChannelSaveTest(unittest.TestCase):

    def test_save_load(self):
        """
        Este test obtiene canales de un shapefile, crea knickpoints y regressiones, los salva y los
        vuelve a cargar
        """
        # Cargamos objetos y obtenemos canales
        net = Network("{}/jebja30_net.dat".format(outfolder))
        shp = "{}/jebja_channels.shp".format(infolder)
        canales = shp_to_channels(shp, net, "segid")
        
        for canal in canales:
            # Creamos 5 knickpoints aleatorios
            for n in range(5):
                index = np.random.randint(0, canal._ix.size)
                tipo = np.random.choice([0, 1, 2, 3])
                canal.addKP(index, tipo)
            # Creamos una regression
            canal.addRegression(0, canal._ix.size-1)
            path = outfolder + "/jebja_chan{}.dat".format(canal.getOid())
            
            # Guardamos canal creado
            canal.save(path)
            
            # Cargamos canal creado
            auxch = Channel(path)
            
            # Verificamos algunas propiedades para ver si se guardó-cargó bien
            self.assertEqual(np.array_equal(canal._ix, auxch._ix), True)
            self.assertEqual(np.array_equal(canal._ax, auxch._ax), True)
            self.assertEqual(np.array_equal(canal._dx, auxch._dx), True)
            self.assertEqual(np.array_equal(canal._zx, auxch._zx), True)
            self.assertEqual(np.array_equal(canal._kp, auxch._kp), True)

if __name__ == "__main__":
    unittest.main()