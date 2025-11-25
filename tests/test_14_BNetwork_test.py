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
from matplotlib import pyplot as plt

from landspy  import Flow, Basin, Network, BNetwork, DEM, Grid
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

    def test_BNetwork(self):
        """
        Test00 Crea Bnetwork y testea varias de sus funciones principales:
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

            heads = np.loadtxt("{}/{}_bnet_heads.txt".format(infolder, file), delimiter=";")
            outlets = net.snapPoints(outlets)
            
            for bid in np.unique(cuencas.readArray()):
                if bid == 0:
                    continue
                # Creamos objeto BNetwork
                bnet = BNetwork(net, cuencas, heads, bid)

               # Test canal principal
                mainc = bnet.getChannels(1)[0]
                xy = mainc.getXY()
                outf = open("{}/{}_tempbasin{}.txt".format(outfolder, file, bid), "w")
                outf.write("X;Y\n")
                for row in xy:
                    outf.write("{};{}\n".format(row[0], row[1]))
                outf.close()

                # ChiPlot (relativo y absoluto)
                bnet.calculateChi()
                fig = plt.figure()
                ax1 = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                bnet.chiPlot(ax1, relative=True)
                bnet.chiPlot(ax2)
                ax1.set_title("Chi plot (z-relative)")
                ax2.set_title("Chi plot (z-absolute)")
                fig.tight_layout()
                fig.savefig("{}/{}_chiPlot_{}.png".format(outfolder, file, bid))
                plt.close(fig)

                # Sensitivity Analysis
                bestheta, r2, fig = bnet.chiSensitivityAnalysis(0.2, 0.65, 0.01, True)
                fig.savefig("{}/{}_chiSensitivity_{}.png".format(outfolder, file, bid))
                plt.close(fig)

                # Print ChiPlots for sensitivity analysis
                star = bestheta - 0.02 * 4
                stop = bestheta + 0.02 * 5
                fig = plt.figure()
                for n in range(9):
                    ax = plt.subplot(3, 3, n + 1)
                    bnet.calculateChi(star)
                    bnet.chiPlot(ax, relative=True)
                    star += 0.02
                    ax.set_title("m/n = {0:.2f}".format(star))
                fig.tight_layout()
                fig.savefig("{}/{}_chiSensitivityPlot_{}.png".format(outfolder, file, bid))
                plt.close(fig)



if __name__ == "__main__":
    unittest.main()