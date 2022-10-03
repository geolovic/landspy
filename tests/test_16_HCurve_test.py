#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 February, 2020
Testing suite for HCurve class
@author: J. Vicente Perez
@email: geolovic@gmail.com
@date: 01 October 2022
@last_modified: 02 October 2022
"""

import unittest
import numpy as np
import sys, os
# Add to the path code folder and data folder
sys.path.append("../src/")
from landspy import DEM, Basin, Grid, Flow, HCurve
infolder = "data/in"
outfolder = "data/out"

class HCurveTest(unittest.TestCase):

    def setUp(self):
        # Create the datasets (if do not exists)
        self.files = ["jebja30", "tunez"]
        for file in self.files:
            flw_path = "{}/{}_fd.tif".format(outfolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem_path = "{}/{}.tif".format(infolder, file)
            dem = DEM(dem_path)
            if not os.path.exists(flw_path):
                fd = Flow(dem)
                fd.save(flw_path)
            else:
                fd = Flow(flw_path)

            if not os.path.exists(basin_path):
                basins = fd.drainageBasins(min_area=0.0025)
                basins.save(basin_path)
            else:
                basins = Grid(basin_path)

    def test_HCurve_01(self):
        # Create hypsometric curves from a DEM and a basin Grid
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor y una aleatoria
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]
            idrdn = np.random.choice(bids)
            
            bids = [idmax, idrdn]
            lbls = ["max", "random"]

            for n, bid in enumerate(bids):
                hcurve = HCurve(dem, basins, bid)
                # Check is HCurve instance
                self.assertIsInstance(hcurve, HCurve)
                # Check all moments
                self.assertNotEqual(hcurve.getHI(), 0)
                self.assertNotEqual(hcurve.getHI2(), 0)
                self.assertNotEqual(hcurve.getKurtosis(), 0)
                self.assertNotEqual(hcurve.getSkewness(), 0)
                self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
                self.assertNotEqual(hcurve.getDensitySkewness(), 0)


    def test_HCurve_01b(self):
        # Create hipsometric curves from a DEM and a path to a basins tiff file
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor y una aleatoria
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]
            idrdn = np.random.choice(bids)

            bids = [idmax, idrdn]
            lbls = ["max", "random"]

            for n, bid in enumerate(bids):
                hcurve = HCurve(dem, basin_path, bid)
                # Check is HCurve instance
                self.assertIsInstance(hcurve, HCurve)
                # Check all moments
                self.assertNotEqual(hcurve.getHI(), 0)
                self.assertNotEqual(hcurve.getHI2(), 0)
                self.assertNotEqual(hcurve.getKurtosis(), 0)
                self.assertNotEqual(hcurve.getSkewness(), 0)
                self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
                self.assertNotEqual(hcurve.getDensitySkewness(), 0)

    def test_HCurve_02(self):
        # Create hipsometric curves from a Basin
        for file in self.files:
            dem_path = "{}/{}.tif".format(infolder, file)
            basin_path = "{}/{}_basins.tif".format(outfolder, file)
            dem = DEM(dem_path)
            basins = Grid(basin_path)

            # Obtenemos la cuenca mayor
            bids, counts = np.unique(basins.readArray(), return_counts=True)
            if 0 in bids:
                bids = bids[1:]
                counts = counts[1:]

            idmax = bids[np.argmax(counts)]

            # Creamos objeto Cuenca (Basin)
            basin = Basin(dem, basins, idmax)

            # Creamos HCurve
            hcurve = HCurve(basin)
            # Check is HCurve instance
            self.assertIsInstance(hcurve, HCurve)
            # Check all moments
            self.assertNotEqual(hcurve.getHI(), 0)
            self.assertNotEqual(hcurve.getHI2(), 0)
            self.assertNotEqual(hcurve.getKurtosis(), 0)
            self.assertNotEqual(hcurve.getSkewness(), 0)
            self.assertNotEqual(hcurve.getDensityKurtosis(), 0)
            self.assertNotEqual(hcurve.getDensitySkewness(), 0)


if __name__ == "__main__":
    unittest.main()