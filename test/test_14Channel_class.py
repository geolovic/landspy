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
from topopy import Network, Channel
import ogr
infolder = "data/in"
outfolder = "data/out"


class ChannelClassTest(unittest.TestCase):

    
    def test_get_channel(self, plot_canales=False):
        """
        Test Obtiene shapefile con canales, selecciona 5 grupos de 10 canales
        aleaorios y crea objeto Channel (con puntos inicial y final de línea)
        """
        files = ["small25", "jebja30", "tunez"]
        for file in files:
            
            net = net = Network("{}/{}_net.dat".format(infolder, file))
            shp_path = outfolder + "/canales_{}.shp".format(file)
            net.export_to_shp(shp_path, False)
            
            dataset = ogr.Open(shp_path)
            layer = dataset.GetLayer(0)
            canales = []
            for n in range(5):
                for n in np.random.randint(0, layer.GetFeatureCount(), 10):
                    feat = layer.GetFeature(n)
                    geom = feat.GetGeometryRef()
                    head = geom.GetPoint(0)
                    mouth = geom.GetPoint(geom.GetPointCount() - 1)
                    canal = net.get_channel(head, mouth)
                    canales.append(canal)
                    # Verificamos que canal se ha creado bien
                    self.assertEqual(isinstance(canal, Channel), True)
           
    def test_get_channel2(self):
        """
        Test crea los canales hasta outlet para todas las cabeceras (sin mouth)
        """
        files = ["small25", "jebja30", "tunez"]

        for file in files:
            net = Network("data/in/{}_net.dat".format(file))
            heads = net.get_stream_poi(kind="heads", coords="XY")
            for head in heads:
                canal = net.get_channel(head)
                # Verificamos que canal se ha creado bien
                self.assertEqual(isinstance(canal, Channel), True)

if __name__ == "__main__":
    unittest.main()