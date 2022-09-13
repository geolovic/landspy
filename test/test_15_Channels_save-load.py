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
import numpy as np
from topopy import Network, Channel, shp_to_channels
import ogr
infolder = "data/in"
outfolder = "data/out"


class ChannelSaveTest(unittest.TestCase):

    def test_save_load(self):
        """
        Este test obtiene canales de un shapefile, crea knickpoints y regressiones, los salva y los
        vuelve a cargar
        """
        # Cargamos objetos y obtenemos canales
        net = Network("data/in/jebja30_net.dat")
        shp = "data/in/jebja_channels.shp"
        canales = shp_to_channels(shp, net, "segid")
        
        for canal in canales:
            # Creamos 5 knickpoints aleatorios
            for n in range(5):
                index = np.random.randint(0, canal._ix.size)
                tipo = np.random.choice([0, 1, 2, 3])
                canal.add_kp(index, tipo)
            # Creamos una regression
            canal.add_regression(0, canal._ix.size-1)
            path = outfolder + "/jebja_chan{}.dat".format(canal.get_oid())
            
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