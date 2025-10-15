#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 february, 2021
Testing suite for BNetwork class
@author: J. Vicente Perez
@email: geolovic@hotmail.com
@date: 09 february, 2021
@last_modified: 15 october, 2025
"""

import unittest
import numpy as np
from landspy  import Flow, Basin, Network, BNetwork, DEM
from landspy.network import NetworkError
from osgeo import ogr, osr

import sys, os
# Forzar el directorio actual al del archivo
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
infolder = "data/in"
outfolder = "data/out"

def canales_to_shapefile(path, canales):
    # Creamos shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataset = driver.CreateDataSource(path)
    sp = osr.SpatialReference()
    proj = canales[0].getCRS()
    sp.ImportFromWkt(proj)
    
    # Creamos layer
    layer = dataset.CreateLayer("canales", sp, ogr.wkbLineString)
    
    # Add fields
    campos = ["oid", "name", "flowto", "thetaref", "chi0", "slp_np", "ksn_np"]
    tipos = [0, 4, 0, 2, 2, 0, 0]
    
    for n in range(len(tipos)):
        layer.CreateField(ogr.FieldDefn(campos[n], tipos[n]))
    
    # Add channels to shapefile
    for canal in canales:
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetField("oid", canal.getOid())
        feat.SetField("name", canal.getName())
        feat.SetField("flowto", canal.getFlow())
        feat.SetField("thetaref", canal._thetaref)
        feat.SetField("chi0", canal._chi0)
        feat.SetField("slp_np", canal._slp_np)
        feat.SetField("ksn_np", canal._ksn_np)
        geom = ogr.Geometry(ogr.wkbLineString)
        xy = canal.getXY()
        for row in xy:
            geom.AddPoint(row[0], row[1])
            
        feat.SetGeometry(geom)
        layer.CreateFeature(feat)
        

class BNetworkGetMainChannelTest(unittest.TestCase):

    def test_BNetwork_get_main_channel(self):
        """
        Test00 Crea Bnetwork y extrae main channel
        Los canales principales se guardarán en la carpeta "data/out" como archivos de puntos y pueden
        verificarse con QGIS. 
        """
        files = ["jebja30", "tunez"]
        for file in files:
            # Cargamos DEM, Flow, Network
            fd = Flow("{}/{}_fd.tif".format(outfolder, file))
            net = Network("{}/{}_net.dat".format(outfolder, file))
            
            # Cargamos outlets y generamos cuencas
            outlets = np.loadtxt("{}/{}_bnet_outlets.txt".format(infolder, file), delimiter=";")
            outlets = net.snapPoints(outlets)
            cuencas = fd.drainageBasins(outlets)
            cuencas.save("{}/{}_tempbasins.tif".format(outfolder, file))
            
            heads = np.loadtxt("{}/{}_bnet_heads.txt".format(infolder, file), delimiter=";")
            outlets = net.snapPoints(outlets)
            
            for bid in np.unique(cuencas.readArray()):
                if bid == 0:
                    continue
                bnet = BNetwork(net, cuencas, heads, bid)
                # Extraemos canal principal
                mainc = bnet.getChannels(1)[0]
                # Guardamos canal principal
                xy = mainc.getXY()
                outf = open("{}/{}_tempbasin{}.txt".format(outfolder, file, bid), "w")
                outf.write("X;Y\n")
                for row in xy:
                    outf.write("{};{}\n".format(row[0], row[1]))
                outf.close()
   

if __name__ == "__main__":
    unittest.main()